use std::{convert::TryFrom, path::Path};

use futures::TryStreamExt;
use interval_tree::IntervalTree;
use noodles::{
    bam, gff,
    sam::{self, header::ReferenceSequences},
};
use tokio::{fs::File, io};

use crate::{count::get_tree, Entry, Features, PairPosition, StrandSpecification};

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

impl From<sam::record::Flags> for Strand {
    fn from(flags: sam::record::Flags) -> Self {
        if flags.is_reverse_complemented() {
            Self::Reverse
        } else {
            Self::Forward
        }
    }
}

impl TryFrom<gff::record::Strand> for Strand {
    type Error = ();

    fn try_from(strand: gff::record::Strand) -> Result<Self, Self::Error> {
        match strand {
            gff::record::Strand::Forward => Ok(Self::Forward),
            gff::record::Strand::Reverse => Ok(Self::Reverse),
            _ => Err(()),
        }
    }
}

fn count_paired_end_record(
    counts: &mut Counts,
    tree: &IntervalTree<u64, Entry>,
    record: &bam::Record,
) -> io::Result<()> {
    let start = record
        .position()
        .map(i32::from)
        .map(|n| n as u64)
        .expect("record cannot be unmapped");
    let reference_len = record.cigar().reference_len()? as u64;
    let end = start + reference_len - 1;

    let pair_position = PairPosition::try_from(record)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let record_strand = Strand::from(record.flags());

    for entry in tree.find(start..=end) {
        let strand = entry.get().1;

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
    let start = record
        .position()
        .map(i32::from)
        .map(|n| n as u64)
        .expect("record cannot be unmapped");
    let reference_len = record.cigar().reference_len()? as u64;
    let end = start + reference_len - 1;

    let record_strand = Strand::from(record.flags());

    for entry in tree.find(start..=end) {
        let strand = entry.get().1;

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

pub async fn detect_specification<P>(
    src: P,
    reference_sequences: &ReferenceSequences,
    features: &Features,
) -> io::Result<(LibraryLayout, StrandSpecification, f64)>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).await.map(bam::AsyncReader::new)?;
    reader.read_header().await?;
    reader.read_reference_sequences().await?;

    let mut counts = Counts::default();

    let mut records = reader.records();
    let mut n = 0;

    while let Some(record) = records.try_next().await? {
        if n >= MAX_RECORDS {
            break;
        }

        let flags = record.flags();

        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        let tree = match get_tree(
            features,
            reference_sequences,
            record.reference_sequence_id(),
        )? {
            Some(t) => t,
            None => continue,
        };

        if flags.is_paired() {
            counts.paired += 1;
            count_paired_end_record(&mut counts, tree, &record)?;
        } else {
            count_single_end_record(&mut counts, tree, &record)?;
        }

        n += 1;
    }

    let library_layout = if counts.paired > 0 {
        LibraryLayout::PairedEnd
    } else {
        LibraryLayout::SingleEnd
    };

    if counts.matches == 0 {
        return Ok((library_layout, StrandSpecification::None, 0.0));
    }

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
