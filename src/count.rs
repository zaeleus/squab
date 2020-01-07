use std::collections::{HashMap, HashSet};
use std::io;

use interval_tree::IntervalTree;
use noodles_bam::{self as bam, Record, Reference};
use noodles_gff as gff;

use crate::{CigarToIntervals, Entry, Features, PairPosition, RecordPairs};

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

    pub fn filter(&self, ctx: &mut Context, record: &Record) -> io::Result<bool> {
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

    pub fn filter_pair(&self, ctx: &mut Context, r1: &Record, r2: &Record) -> io::Result<bool> {
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
    references: &[Reference],
    filter: &Filter,
    strand_irrelevant: bool,
) -> io::Result<Context>
where
    I: Iterator<Item = io::Result<Record>>,
{
    let mut ctx = Context::default();

    for result in records {
        let record = result?;

        count_single_end_record(
            &mut ctx,
            features,
            references,
            filter,
            strand_irrelevant,
            &record,
        )?;
    }

    Ok(ctx)
}

pub fn count_single_end_record(
    ctx: &mut Context,
    features: &Features,
    references: &[Reference],
    filter: &Filter,
    strand_irrelevant: bool,
    record: &Record,
) -> io::Result<()> {
    if filter.filter(ctx, record)? {
        return Ok(());
    }

    let reference = get_reference(&references, record.ref_id())?;

    let cigar = record.cigar();
    let start = record.pos() as u64;
    let flag = record.flag();
    let intervals = CigarToIntervals::new(&cigar, start, flag, false);

    let name = reference.name();
    let tree = match features.get(name) {
        Some(t) => t,
        None => {
            ctx.no_feature += 1;
            return Ok(());
        }
    };

    let set = find(tree, intervals, strand_irrelevant);

    update_intersections(ctx, set);

    Ok(())
}

pub fn count_paired_end_records<I>(
    records: I,
    features: Features,
    references: Vec<Reference>,
    filter: Filter,
    strand_irrelevant: bool,
) -> io::Result<Context>
where
    I: Iterator<Item = io::Result<Record>>,
{
    let mut ctx = Context::default();

    let primary_only = !filter.with_secondary_records && !filter.with_supplementary_records;
    let mut pairs = RecordPairs::new(records, primary_only);

    for pair in &mut pairs {
        let (r1, r2) = pair?;

        if filter.filter_pair(&mut ctx, &r1, &r2)? {
            continue;
        }

        let reference = get_reference(&references, r1.ref_id())?;

        let cigar = r1.cigar();
        let start = r1.pos() as u64;
        let f1 = r1.flag();
        let intervals = CigarToIntervals::new(&cigar, start, f1, false);

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            }
        };

        let mut set = find(tree, intervals, strand_irrelevant);

        let reference = get_reference(&references, r2.ref_id())?;

        let cigar = r2.cigar();
        let start = r2.pos() as u64;
        let f2 = r2.flag();
        let intervals = CigarToIntervals::new(&cigar, start, f2, true);

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            }
        };

        let set2 = find(tree, intervals, strand_irrelevant);

        set.extend(set2.into_iter());

        update_intersections(&mut ctx, set);
    }

    for record in pairs.singletons() {
        if filter.filter(&mut ctx, &record)? {
            continue;
        }

        let cigar = record.cigar();
        let start = record.pos() as u64;

        let reverse = match PairPosition::from(&record) {
            PairPosition::First => false,
            PairPosition::Second => true,
        };

        let flag = record.flag();
        let intervals = CigarToIntervals::new(&cigar, start, flag, reverse);

        let reference = get_reference(&references, record.ref_id())?;

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            }
        };

        let set = find(tree, intervals, strand_irrelevant);

        update_intersections(&mut ctx, set);
    }

    Ok(ctx)
}

fn find(
    tree: &IntervalTree<u64, Entry>,
    intervals: CigarToIntervals,
    strand_irrelevant: bool,
) -> HashSet<String> {
    let mut set = HashSet::new();

    for (interval, is_reverse) in intervals {
        for entry in tree.find(interval.clone()) {
            let gene_name = &entry.value.0;
            let strand = &entry.value.1;

            if strand_irrelevant
                || (strand == &gff::Strand::Reverse && is_reverse)
                || (strand == &gff::Strand::Forward && !is_reverse)
            {
                set.insert(gene_name.to_string());
            }
        }
    }

    set
}

fn is_nonunique_record(record: &Record) -> io::Result<bool> {
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

fn get_reference<'a>(references: &'a [Reference], ref_id: i32) -> io::Result<&'a Reference> {
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
