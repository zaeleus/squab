use std::collections::{HashMap, HashSet};
use std::io::{self, Read};
use std::ops::Range;

use interval_tree::IntervalTree;
use noodles::formats::bam::{self, ByteRecord, Flag, Reference};
use {Entry, Features, PairPosition, RecordPairs, Strand, cigar_to_intervals};

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

        let intervals = cigar_to_intervals(&record, false);

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            },
        };

        let set = find(tree, &intervals);

        if set.len() == 0 {
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

        let intervals = cigar_to_intervals(&r1, false);

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            },
        };

        let mut set = find(tree, &intervals);

        let ref_id = r2.ref_id();

        if ref_id < 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("expected ref id >= 0, got {}", ref_id),
            ));
        }

        let reference = &references[ref_id as usize];

        let intervals = cigar_to_intervals(&r2, true);

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            },
        };

        let set2 = find(tree, &intervals);

        set.extend(set2.into_iter());

        if set.len() == 0 {
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

        let intervals = match PairPosition::from(&record) {
            PairPosition::First => cigar_to_intervals(&record, false),
            PairPosition::Second => cigar_to_intervals(&record, true),
        };

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

        let set = find(tree, &intervals);

        if set.len() == 0 {
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
    intervals: &[(Range<u64>, bool)],
) -> HashSet<String> {
    let mut set = HashSet::new();

    for (interval, is_reverse) in intervals {
        for entry in tree.find(interval.clone()) {
            let gene_name = &entry.value.0;
            let strand = &entry.value.1;

            if (strand == &Strand::Reverse && *is_reverse)
                    || (strand == &Strand::Forward && !*is_reverse) {
                set.insert(gene_name.to_string());
            }
        }
    }

    set
}
