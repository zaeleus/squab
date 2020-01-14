use std::{fs::File, io, path::Path};

use interval_tree::IntervalTree;
use noodles_bam as bam;
use noodles_gff as gff;

use crate::{count::get_tree, Context, Entry, Features, PairPosition, StrandSpecification};

const MAX_RECORDS: usize = 524288;
const STRANDEDNESS_THRESHOLD: f64 = 0.75;

#[derive(Debug, Default)]
struct Counts {
    paired: u64,
    matches: u64,
    forward: u64,
    reverse: u64,
}

fn count_paired_end_record(
    counts: &mut Counts,
    tree: &IntervalTree<u64, Entry>,
    record: &bam::Record,
) {
    let flag = record.flag();
    let start = record.pos() as u64;
    let end = start + record.cigar().mapped_len() as u64;
    let pair_position = PairPosition::from(record);

    for entry in tree.find(start..end) {
        let strand = &entry.value.1;

        match (pair_position, flag.is_reverse(), strand) {
            (PairPosition::First, false, gff::Strand::Forward)
            | (PairPosition::First, true, gff::Strand::Reverse)
            | (PairPosition::Second, false, gff::Strand::Reverse)
            | (PairPosition::Second, true, gff::Strand::Forward) => {
                counts.forward += 1;
                counts.matches += 1;
            }
            (PairPosition::First, false, gff::Strand::Reverse)
            | (PairPosition::First, true, gff::Strand::Forward)
            | (PairPosition::Second, false, gff::Strand::Forward)
            | (PairPosition::Second, true, gff::Strand::Reverse) => {
                counts.reverse += 1;
                counts.matches += 1;
            }
            (_, _, _) => unimplemented!(),
        }
    }
}

fn count_single_end_record(
    counts: &mut Counts,
    tree: &IntervalTree<u64, Entry>,
    record: &bam::Record,
) {
    let flag = record.flag();
    let start = record.pos() as u64;
    let end = start + record.cigar().mapped_len() as u64;

    for entry in tree.find(start..end) {
        let strand = &entry.value.1;

        match (flag.is_reverse(), strand) {
            (false, gff::Strand::Forward) | (true, gff::Strand::Reverse) => {
                counts.forward += 1;
                counts.matches += 1;
            }
            (false, gff::Strand::Reverse) | (true, gff::Strand::Forward) => {
                counts.reverse += 1;
                counts.matches += 1;
            }
            (_, _) => unimplemented!(),
        }
    }
}

pub fn detect_specification<P>(
    src: P,
    references: &[bam::Reference],
    features: &Features,
) -> io::Result<(bool, StrandSpecification, f64)>
where
    P: AsRef<Path>,
{
    let file = File::open(src)?;
    let mut reader = bam::Reader::new(file);
    reader.meta()?;

    let mut counts = Counts::default();
    let mut _ctx = Context::default();

    for result in reader.records().take(MAX_RECORDS) {
        let record = result?;
        let flag = record.flag();

        if flag.is_unmapped() || flag.is_secondary() || flag.is_supplementary() {
            continue;
        }

        let ref_id = record.ref_id();
        let tree = match get_tree(&mut _ctx, features, references, ref_id)? {
            Some(t) => t,
            None => continue,
        };

        if flag.is_paired() {
            counts.paired += 1;
            count_paired_end_record(&mut counts, tree, &record);
        } else {
            count_single_end_record(&mut counts, tree, &record);
        }
    }

    let is_paired_end = counts.paired > 0;

    let matches = counts.matches as f64;
    let forward_pct = counts.forward as f64 / matches;
    let reverse_pct = counts.reverse as f64 / matches;

    let (strand_specification, strandedness_confidence) = if forward_pct > STRANDEDNESS_THRESHOLD {
        (StrandSpecification::Forward, forward_pct)
    } else if reverse_pct > STRANDEDNESS_THRESHOLD {
        (StrandSpecification::Reverse, reverse_pct)
    } else {
        (StrandSpecification::None, forward_pct + reverse_pct)
    };

    Ok((is_paired_end, strand_specification, strandedness_confidence))
}
