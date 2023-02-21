use std::{fs::File, io, path::Path};

use interval_tree::IntervalTree;
use noodles::{
    bam,
    core::Position,
    gff,
    sam::{self, header::ReferenceSequences},
};

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
    tree: &IntervalTree<Position, Entry>,
    flags: sam::record::Flags,
    start: Position,
    end: Position,
) -> io::Result<()> {
    let pair_position =
        PairPosition::try_from(flags).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let record_strand = Strand::from(flags);

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
    tree: &IntervalTree<Position, Entry>,
    flags: sam::record::Flags,
    start: Position,
    end: Position,
) -> io::Result<()> {
    let record_strand = Strand::from(flags);

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

pub fn detect_specification<P>(
    src: P,
    reference_sequences: &ReferenceSequences,
    features: &Features,
) -> io::Result<(LibraryLayout, StrandSpecification, f64)>
where
    P: AsRef<Path>,
{
    use crate::record::alignment_end;

    let mut reader = File::open(src).map(bam::Reader::new)?;
    reader.read_header()?;
    reader.read_reference_sequences()?;

    let mut counts = Counts::default();

    for result in reader.lazy_records().take(MAX_RECORDS) {
        let record = result?;

        let flags = record.flags()?;

        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        let reference_sequence_id = record.reference_sequence_id()?;

        let tree = match get_tree(features, reference_sequences, reference_sequence_id)? {
            Some(t) => t,
            None => continue,
        };

        let alignment_start = record.alignment_start()?;
        let cigar = sam::record::Cigar::try_from(record.cigar())?;

        let start = alignment_start.expect("missing alignment start");
        let end = alignment_end(alignment_start, &cigar).expect("missing alignment end");

        if flags.is_segmented() {
            counts.paired += 1;
            count_paired_end_record(&mut counts, tree, flags, start, end)?;
        } else {
            count_single_end_record(&mut counts, tree, flags, start, end)?;
        }
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
