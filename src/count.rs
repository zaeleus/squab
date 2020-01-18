use std::{
    collections::{HashMap, HashSet},
    convert::TryFrom,
    io,
};

use interval_tree::IntervalTree;
use noodles_bam as bam;
use noodles_gff as gff;

use crate::{CigarToIntervals, Entry, Features, PairPosition, RecordPairs, StrandSpecification};

#[derive(Clone)]
pub struct Filter {
    min_mapq: u8,
    with_secondary_records: bool,
    with_supplementary_records: bool,
    with_nonunique_records: bool,
}

impl Filter {
    pub fn new(
        min_mapq: u8,
        with_secondary_records: bool,
        with_supplementary_records: bool,
        with_nonunique_records: bool,
    ) -> Filter {
        Self {
            min_mapq,
            with_secondary_records,
            with_supplementary_records,
            with_nonunique_records,
        }
    }

    pub fn filter(&self, ctx: &mut Context, record: &bam::Record) -> io::Result<bool> {
        let flag = record.flag();

        if flag.is_unmapped() {
            ctx.unmapped += 1;
            return Ok(true);
        }

        if (!self.with_secondary_records && flag.is_secondary())
            || (!self.with_supplementary_records && flag.is_supplementary())
        {
            return Ok(true);
        }

        if !self.with_nonunique_records && is_nonunique_record(&record)? {
            ctx.nonunique += 1;
            return Ok(true);
        }

        if record.mapq() < self.min_mapq {
            ctx.low_quality += 1;
            return Ok(true);
        }

        Ok(false)
    }

    pub fn filter_pair(
        &self,
        ctx: &mut Context,
        r1: &bam::Record,
        r2: &bam::Record,
    ) -> io::Result<bool> {
        let f1 = r1.flag();
        let f2 = r2.flag();

        if f1.is_unmapped() && f2.is_unmapped() {
            ctx.unmapped += 1;
            return Ok(true);
        }

        if (!self.with_secondary_records && (f1.is_secondary() || f2.is_secondary()))
            || (!self.with_supplementary_records
                && (f1.is_supplementary() || f2.is_supplementary()))
        {
            return Ok(true);
        }

        if !self.with_nonunique_records && (is_nonunique_record(&r1)? || is_nonunique_record(&r2)?)
        {
            ctx.nonunique += 1;
            return Ok(true);
        }

        if r1.mapq() < self.min_mapq || r2.mapq() < self.min_mapq {
            ctx.low_quality += 1;
            return Ok(true);
        }

        Ok(false)
    }
}

#[derive(Default)]
pub struct Context {
    pub counts: HashMap<String, u64>,
    pub no_feature: u64,
    pub ambiguous: u64,
    pub low_quality: u64,
    pub unmapped: u64,
    pub nonunique: u64,
}

impl Context {
    pub fn add(&mut self, other: &Context) {
        for (name, count) in other.counts.iter() {
            let entry = self.counts.entry(name.to_string()).or_insert(0);
            *entry += count;
        }

        self.no_feature += other.no_feature;
        self.ambiguous += other.ambiguous;
        self.low_quality += other.low_quality;
        self.unmapped += other.unmapped;
        self.nonunique += other.nonunique;
    }
}

pub fn count_single_end_records<I>(
    records: I,
    features: &Features,
    references: &[bam::Reference],
    filter: &Filter,
    strand_specification: StrandSpecification,
) -> io::Result<Context>
where
    I: Iterator<Item = io::Result<bam::Record>>,
{
    let mut ctx = Context::default();

    for result in records {
        let record = result?;

        count_single_end_record(
            &mut ctx,
            features,
            references,
            filter,
            strand_specification,
            &record,
        )?;
    }

    Ok(ctx)
}

pub fn count_single_end_record(
    ctx: &mut Context,
    features: &Features,
    references: &[bam::Reference],
    filter: &Filter,
    strand_specification: StrandSpecification,
    record: &bam::Record,
) -> io::Result<()> {
    if filter.filter(ctx, record)? {
        return Ok(());
    }

    let cigar = record.cigar();
    let start = record.pos() as u64;
    let flag = record.flag();

    let reverse = match strand_specification {
        StrandSpecification::Reverse => true,
        _ => false,
    };

    let intervals = CigarToIntervals::new(&cigar, start, flag, reverse);

    let ref_id = record.ref_id();
    let tree = match get_tree(ctx, features, references, ref_id)? {
        Some(t) => t,
        None => return Ok(()),
    };

    let set = find(tree, intervals, strand_specification);

    update_intersections(ctx, set);

    Ok(())
}

pub fn count_paired_end_records<I>(
    records: I,
    features: &Features,
    references: &[bam::Reference],
    filter: &Filter,
    strand_specification: StrandSpecification,
) -> io::Result<(Context, RecordPairs<I>)>
where
    I: Iterator<Item = io::Result<bam::Record>>,
{
    let mut ctx = Context::default();

    let primary_only = !filter.with_secondary_records && !filter.with_supplementary_records;
    let mut pairs = RecordPairs::new(records, primary_only);

    for pair in &mut pairs {
        let (r1, r2) = pair?;

        if filter.filter_pair(&mut ctx, &r1, &r2)? {
            continue;
        }

        let cigar = r1.cigar();
        let start = r1.pos() as u64;
        let f1 = r1.flag();

        let reverse = match strand_specification {
            StrandSpecification::Reverse => true,
            _ => false,
        };

        let intervals = CigarToIntervals::new(&cigar, start, f1, reverse);

        let ref_id = r1.ref_id();
        let tree = match get_tree(&mut ctx, features, references, ref_id)? {
            Some(t) => t,
            None => continue,
        };

        let mut set = find(tree, intervals, strand_specification);

        let cigar = r2.cigar();
        let start = r2.pos() as u64;
        let f2 = r2.flag();

        let reverse = match strand_specification {
            StrandSpecification::Reverse => false,
            _ => true,
        };

        let intervals = CigarToIntervals::new(&cigar, start, f2, reverse);

        let ref_id = r2.ref_id();
        let tree = match get_tree(&mut ctx, features, references, ref_id)? {
            Some(t) => t,
            None => continue,
        };

        let set2 = find(tree, intervals, strand_specification);

        set.extend(set2.into_iter());

        update_intersections(&mut ctx, set);
    }

    Ok((ctx, pairs))
}

pub fn count_paired_end_record_singletons<I>(
    records: I,
    features: &Features,
    references: &[bam::Reference],
    filter: &Filter,
    strand_specification: StrandSpecification,
) -> io::Result<Context>
where
    I: Iterator<Item = io::Result<bam::Record>>,
{
    let mut ctx = Context::default();

    for result in records {
        let record = result?;

        if filter.filter(&mut ctx, &record)? {
            continue;
        }

        let cigar = record.cigar();
        let start = record.pos() as u64;

        let reverse = match PairPosition::try_from(&record) {
            Ok(PairPosition::First) => false,
            Ok(PairPosition::Second) => true,
            Err(_) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "record is neither read 1 nor 2",
                ))
            }
        };

        let reverse = match strand_specification {
            StrandSpecification::Reverse => !reverse,
            _ => reverse,
        };

        let flag = record.flag();
        let intervals = CigarToIntervals::new(&cigar, start, flag, reverse);

        let ref_id = record.ref_id();
        let tree = match get_tree(&mut ctx, features, references, ref_id)? {
            Some(t) => t,
            None => continue,
        };

        let set = find(tree, intervals, strand_specification);

        update_intersections(&mut ctx, set);
    }

    Ok(ctx)
}

fn find(
    tree: &IntervalTree<u64, Entry>,
    intervals: CigarToIntervals,
    strand_specification: StrandSpecification,
) -> HashSet<String> {
    let mut set = HashSet::new();

    for (interval, is_reverse) in intervals {
        for entry in tree.find(interval.clone()) {
            let gene_name = &entry.value.0;
            let strand = &entry.value.1;

            match strand_specification {
                StrandSpecification::None => {
                    set.insert(gene_name.to_string());
                }
                StrandSpecification::Forward | StrandSpecification::Reverse => {
                    if (strand == &gff::Strand::Reverse && is_reverse)
                        || (strand == &gff::Strand::Forward && !is_reverse)
                    {
                        set.insert(gene_name.to_string());
                    }
                }
            }
        }
    }

    set
}

fn is_nonunique_record(record: &bam::Record) -> io::Result<bool> {
    let data = record.data();

    for result in data.fields() {
        let field = result?;

        if field.tag() == "NH" {
            if let bam::data::Value::Int8(n) = field.value() {
                return Ok(*n > 1);
            }
        }
    }

    Ok(false)
}

fn get_reference<'a>(
    references: &'a [bam::Reference],
    ref_id: i32,
) -> io::Result<&'a bam::Reference> {
    if ref_id < 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("expected ref id >= 0, got {}", ref_id),
        ));
    }

    references.get(ref_id as usize).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("expected ref id < {}, got {}", references.len(), ref_id),
        )
    })
}

fn update_intersections(ctx: &mut Context, intersections: HashSet<String>) {
    if intersections.is_empty() {
        ctx.no_feature += 1;
    } else if intersections.len() == 1 {
        for name in intersections {
            let count = ctx.counts.entry(name).or_insert(0);
            *count += 1;
        }
    } else if intersections.len() > 1 {
        ctx.ambiguous += 1;
    }
}

pub fn get_tree<'t>(
    ctx: &mut Context,
    features: &'t Features,
    references: &[bam::Reference],
    ref_id: i32,
) -> io::Result<Option<&'t IntervalTree<u64, Entry>>> {
    let reference = get_reference(&references, ref_id)?;
    let name = reference.name();

    match features.get(name) {
        Some(t) => Ok(Some(t)),
        None => {
            ctx.no_feature += 1;
            Ok(None)
        }
    }
}

#[cfg(test)]
mod tests {
    use noodles_bam as bam;

    use super::*;

    fn build_references() -> Vec<bam::Reference> {
        return vec![
            bam::Reference::new(String::from("chr1"), 7),
            bam::Reference::new(String::from("chr2"), 12),
            bam::Reference::new(String::from("chr3"), 148),
        ];
    }

    #[test]
    fn test_get_reference() {
        let references = build_references();

        let reference = get_reference(&references, 1).unwrap();
        assert_eq!(reference.name(), "chr2");
        assert_eq!(reference.len(), 12);

        let reference = get_reference(&references, -2);
        assert!(reference.is_err());

        let reference = get_reference(&references, 5);
        assert!(reference.is_err());
    }
}
