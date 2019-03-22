use std::collections::{HashMap, HashSet};
use std::io::{self, Read};

use interval_tree::IntervalTree;
use noodles::formats::{
    bam::{self, Record, Reference},
    gff,
};

use crate::{CigarToIntervals, Entry, Features, PairPosition, RecordPairs};

#[derive(Default)]
pub struct Context {
    pub counts: HashMap<String, u64>,
    pub no_feature: u64,
    pub ambiguous: u64,
    pub low_quality: u64,
    pub unmapped: u64,
    pub nonunique: u64,
}

pub fn count_single_end_records<R>(
    mut reader: bam::Reader<R>,
    features: Features,
    references: Vec<Reference>,
    min_mapq: u8,
    with_secondary_records: bool,
    with_supplementary_records: bool,
    with_nonunique_records: bool,
    strand_irrelevant: bool,
) -> io::Result<Context>
where
    R: Read,
{
    let mut ctx = Context::default();
    let mut record = Record::new();

    loop {
        let bytes_read = reader.read_record(&mut record)?;

        if bytes_read == 0 {
            break;
        }

        let flag = record.flag();

        if flag.is_unmapped() {
            ctx.unmapped += 1;
            continue;
        }

        if (!with_secondary_records && flag.is_secondary())
            || (!with_supplementary_records && flag.is_supplementary())
        {
            continue;
        }

        if !with_nonunique_records && is_nonunique_record(&record)? {
            ctx.nonunique += 1;
            continue;
        }

        if record.mapq() < min_mapq {
            ctx.low_quality += 1;
            continue;
        }

        let reference = get_reference(&references, record.ref_id())?;

        let cigar = record.cigar();
        let start = record.pos() as u64;
        let intervals = CigarToIntervals::new(&cigar, start, flag, false);

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            }
        };

        let set = find(tree, intervals, strand_irrelevant);

        if set.is_empty() {
            ctx.no_feature += 1;
        } else if set.len() == 1 {
            for gene_name in set {
                let count = ctx.counts.entry(gene_name).or_insert(0);
                *count += 1;
            }
        } else if set.len() > 1 {
            ctx.ambiguous += 1;
        }
    }

    Ok(ctx)
}

pub fn count_paired_end_records<R>(
    reader: bam::Reader<R>,
    features: Features,
    references: Vec<Reference>,
    min_mapq: u8,
    with_secondary_records: bool,
    with_supplementary_records: bool,
    with_nonunique_records: bool,
    strand_irrelevant: bool,
) -> io::Result<Context>
where
    R: Read,
{
    let mut ctx = Context::default();

    let primary_only = !with_secondary_records && !with_supplementary_records;
    let mut pairs = RecordPairs::new(reader, primary_only);

    for pair in &mut pairs {
        let (r1, r2) = pair?;

        let f1 = r1.flag();
        let f2 = r2.flag();

        if f1.is_unmapped() && f2.is_unmapped() {
            ctx.unmapped += 1;
            continue;
        }

        if (!with_secondary_records && (f1.is_secondary() || f2.is_secondary()))
            || (!with_supplementary_records && (f1.is_supplementary() || f2.is_supplementary()))
        {
            continue;
        }

        if !with_nonunique_records && (is_nonunique_record(&r1)? || is_nonunique_record(&r2)?) {
            ctx.nonunique += 1;
            continue;
        }

        if r1.mapq() < min_mapq || r2.mapq() < min_mapq {
            ctx.low_quality += 1;
            continue;
        }

        let reference = get_reference(&references, r1.ref_id())?;

        let cigar = r1.cigar();
        let start = r1.pos() as u64;
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

        if set.is_empty() {
            ctx.no_feature += 1;
        } else if set.len() == 1 {
            for gene_name in set {
                let count = ctx.counts.entry(gene_name).or_insert(0);
                *count += 1;
            }
        } else if set.len() > 1 {
            ctx.ambiguous += 1;
        }
    }

    for record in pairs.singletons() {
        let flag = record.flag();

        if flag.is_unmapped() {
            ctx.unmapped += 1;
            continue;
        }

        if (!with_secondary_records && flag.is_secondary())
            || (!with_supplementary_records && flag.is_supplementary())
        {
            continue;
        }

        if !with_nonunique_records && is_nonunique_record(&record)? {
            ctx.nonunique += 1;
            continue;
        }

        if record.mapq() < min_mapq {
            ctx.low_quality += 1;
            continue;
        }

        let cigar = record.cigar();
        let start = record.pos() as u64;

        let reverse = match PairPosition::from(&record) {
            PairPosition::First => false,
            PairPosition::Second => true,
        };

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

        if set.is_empty() {
            ctx.no_feature += 1;
        } else if set.len() == 1 {
            for gene_name in set {
                let count = ctx.counts.entry(gene_name).or_insert(0);
                *count += 1;
            }
        } else if set.len() > 1 {
            ctx.ambiguous += 1;
        }
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

#[cfg(test)]
mod tests {
    use noodles::formats::bam;

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
