use std::{convert::TryFrom, fs::File, io, path::Path};

use interval_tree::IntervalTree;
use noodles_bam as bam;
use noodles_gff as gff;

use crate::{count::get_tree, Context, Entry, Features, PairPosition, StrandSpecification};

const MAX_RECORDS: usize = 524_288;
const STRANDEDNESS_THRESHOLD: f64 = 0.75;

#[derive(Debug, Default)]
struct Counts {
    paired: u64,
    matches: u64,
    forward: u64,
    reverse: u64,
}

#[derive(Clone, Copy, Debug)]
pub enum LibraryLayout {
    SingleEnd,
    PairedEnd,
}

#[derive(Clone, Copy, Debug)]
enum Strand {
    Forward,
    Reverse,
}

impl From<bam::Flag> for Strand {
    fn from(flag: bam::Flag) -> Self {
        if flag.is_reverse() {
            Self::Reverse
        } else {
            Self::Forward
        }
    }
}

impl TryFrom<gff::Strand> for Strand {
    type Error = ();

    fn try_from(strand: gff::Strand) -> Result<Self, Self::Error> {
        match strand {
            gff::Strand::Forward => Ok(Self::Forward),
            gff::Strand::Reverse => Ok(Self::Reverse),
            _ => Err(()),
        }
    }
}

fn invalid_record_pair(_: ()) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, "invalid record pair")
}

fn count_paired_end_record(
    counts: &mut Counts,
    tree: &IntervalTree<u64, Entry>,
    record: &bam::Record,
) -> io::Result<()> {
    let start = record.pos() as u64;
    let end = start + record.cigar().mapped_len() as u64;

    let pair_position = PairPosition::try_from(record).map_err(invalid_record_pair)?;
    let record_strand = Strand::from(record.flag());

    for entry in tree.find(start..end) {
        let strand = entry.value.1;

        let feature_strand = match Strand::try_from(strand) {
            Ok(s) => s,
            Err(_) => continue,
        };

        match (pair_position, record_strand, feature_strand) {
            (PairPosition::First, Strand::Forward, Strand::Forward)
            | (PairPosition::First, Strand::Reverse, Strand::Reverse)
            | (PairPosition::Second, Strand::Forward, Strand::Reverse)
            | (PairPosition::Second, Strand::Reverse, Strand::Forward) => {
                counts.forward += 1;
            }
            (PairPosition::First, Strand::Forward, Strand::Reverse)
            | (PairPosition::First, Strand::Reverse, Strand::Forward)
            | (PairPosition::Second, Strand::Forward, Strand::Forward)
            | (PairPosition::Second, Strand::Reverse, Strand::Reverse) => {
                counts.reverse += 1;
            }
        }

        counts.matches += 1;
    }

    Ok(())
}

fn count_single_end_record(
    counts: &mut Counts,
    tree: &IntervalTree<u64, Entry>,
    record: &bam::Record,
) -> io::Result<()> {
    let start = record.pos() as u64;
    let end = start + record.cigar().mapped_len() as u64;

    let record_strand = Strand::from(record.flag());

    for entry in tree.find(start..end) {
        let strand = entry.value.1;

        let feature_strand = match Strand::try_from(strand) {
            Ok(s) => s,
            Err(_) => continue,
        };

        match (record_strand, feature_strand) {
            (Strand::Forward, Strand::Forward) | (Strand::Reverse, Strand::Reverse) => {
                counts.forward += 1;
            }
            (Strand::Forward, Strand::Reverse) | (Strand::Reverse, Strand::Forward) => {
                counts.reverse += 1;
            }
        }

        counts.matches += 1;
    }

    Ok(())
}

pub fn detect_specification<P>(
    src: P,
    references: &[bam::Reference],
    features: &Features,
) -> io::Result<(LibraryLayout, StrandSpecification, f64)>
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
            count_paired_end_record(&mut counts, tree, &record)?;
        } else {
            count_single_end_record(&mut counts, tree, &record)?;
        }
    }

    let library_layout = if counts.paired > 0 {
        LibraryLayout::PairedEnd
    } else {
        LibraryLayout::SingleEnd
    };

    let matches = counts.matches as f64;
    let forward_pct = counts.forward as f64 / matches;
    let reverse_pct = counts.reverse as f64 / matches;

    let (strand_specification, strandedness_confidence) = if forward_pct > STRANDEDNESS_THRESHOLD {
        (StrandSpecification::Forward, forward_pct)
    } else if reverse_pct > STRANDEDNESS_THRESHOLD {
        (StrandSpecification::Reverse, reverse_pct)
    } else {
        let confidence = (0.5 - (forward_pct - reverse_pct).abs()) / 0.5;
        (StrandSpecification::None, confidence)
    };

    Ok((
        library_layout,
        strand_specification,
        strandedness_confidence,
    ))
}
