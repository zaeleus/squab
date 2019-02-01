use std::collections::{HashMap, HashSet};
use std::io::{self, Read};

use interval_tree::IntervalTree;
use noodles::formats::bam::{self, ByteRecord, Cigar, Flag, Reference};

use crate::{CigarToIntervals, Entry, Features, PairPosition, RecordPairs, Strand};

#[derive(Default)]
pub struct Context {
    pub counts: HashMap<String, u64>,
    pub no_feature: u64,
    pub ambiguous: u64,
    pub low_quality: u64,
    pub unmapped: u64,
}

pub fn count_single_end_records<R>(
    mut reader: bam::Reader<R>,
    features: Features,
    references: Vec<Reference>,
    min_mapq: u8,
    with_secondary_records: bool,
    with_supplementary_records: bool,
    strand_irrelevant: bool,
) -> io::Result<Context>
where
    R: Read,
{
    let mut ctx = Context::default();
    let mut record = ByteRecord::new();

    loop {
        let bytes_read = reader.read_byte_record(&mut record)?;

        if bytes_read == 0 {
            break;
        }

        let flag = Flag::new(record.flag());

        if flag.is_unmapped() {
            ctx.unmapped += 1;
            continue;
        }

        if (!with_secondary_records && flag.is_secondary())
                || (!with_supplementary_records && flag.is_supplementary()) {
            continue;
        }

        if record.mapq() < min_mapq {
            ctx.low_quality += 1;
            continue;
        }

        let ref_id = record.ref_id();

        if ref_id < 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("expected ref id >= 0, got {}", ref_id),
            ));
        }

        let reference = &references[ref_id as usize];

        let cigar = Cigar::from_bytes(record.cigar());
        let start = record.pos() as u64;
        let intervals = CigarToIntervals::new(&cigar, start, flag, false);

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            },
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
    strand_irrelevant: bool,
) -> io::Result<Context>
where
    R: Read,
{
    let mut ctx = Context::default();
    let mut pairs = RecordPairs::new(reader);

    for pair in &mut pairs {
        let (r1, r2) = pair?;

        let f1 = Flag::new(r1.flag());
        let f2 = Flag::new(r2.flag());

        if f1.is_unmapped() && f2.is_unmapped() {
            ctx.unmapped += 1;
            continue;
        }

        if (!with_secondary_records && (f1.is_secondary() || f2.is_secondary()))
                || (!with_supplementary_records && (f1.is_supplementary() || f2.is_supplementary())) {
            continue;
        }

        if r1.mapq() < min_mapq || r2.mapq() < min_mapq {
            ctx.low_quality += 1;
            continue;
        }

        let ref_id = r1.ref_id();

        if ref_id < 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("expected ref id >= 0, got {}", ref_id),
            ));
        }

        let reference = &references[ref_id as usize];

        let cigar = Cigar::from_bytes(r1.cigar());
        let start = r1.pos() as u64;
        let intervals = CigarToIntervals::new(&cigar, start, f1, false);

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            },
        };

        let mut set = find(tree, intervals, strand_irrelevant);

        let ref_id = r2.ref_id();

        if ref_id < 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("expected ref id >= 0, got {}", ref_id),
            ));
        }

        let reference = &references[ref_id as usize];

        let cigar = Cigar::from_bytes(r2.cigar());
        let start = r2.pos() as u64;
        let intervals = CigarToIntervals::new(&cigar, start, f2, true);

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            },
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
        let flag = Flag::new(record.flag());

        if flag.is_unmapped() {
            ctx.unmapped += 1;
            continue;
        }

        if (!with_secondary_records && flag.is_secondary())
                || (!with_supplementary_records && flag.is_supplementary()) {
            continue;
        }

        if record.mapq() < min_mapq {
            ctx.low_quality += 1;
            continue;
        }

        let cigar = Cigar::from_bytes(record.cigar());
        let start = record.pos() as u64;

        let reverse = match PairPosition::from(&record) {
            PairPosition::First => false,
            PairPosition::Second => true,
        };

        let intervals = CigarToIntervals::new(&cigar, start, flag, reverse);

        let ref_id = record.ref_id();

        if ref_id < 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("expected ref id >= 0, got {}", ref_id),
            ));
        }

        let reference = &references[ref_id as usize];

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            },
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
                    || (strand == &Strand::Reverse && is_reverse)
                    || (strand == &Strand::Forward && !is_reverse) {
                set.insert(gene_name.to_string());
            }
        }
    }

    set
}
