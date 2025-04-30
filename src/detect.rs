use std::{fmt, fs::File, io, path::Path};

use noodles::{
    bam,
    core::Position,
    gff,
    sam::{self, alignment::Record},
};

use crate::{
    Entry, IntervalTrees, SegmentPosition, StrandSpecification, collections::IntervalTree,
};

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

impl fmt::Display for LibraryLayout {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SingleEnd => "single end".fmt(f),
            Self::PairedEnd => "paired end".fmt(f),
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum Strand {
    Forward,
    Reverse,
}

impl From<sam::alignment::record::Flags> for Strand {
    fn from(flags: sam::alignment::record::Flags) -> Self {
        if flags.is_reverse_complemented() {
            Self::Reverse
        } else {
            Self::Forward
        }
    }
}

impl TryFrom<gff::feature::record::Strand> for Strand {
    type Error = ();

    fn try_from(strand: gff::feature::record::Strand) -> Result<Self, Self::Error> {
        match strand {
            gff::feature::record::Strand::Forward => Ok(Self::Forward),
            gff::feature::record::Strand::Reverse => Ok(Self::Reverse),
            _ => Err(()),
        }
    }
}

fn count_paired_end_record(
    counts: &mut Counts,
    interval_tree: &IntervalTree<Position, Entry>,
    flags: sam::alignment::record::Flags,
    start: Position,
    end: Position,
) -> io::Result<()> {
    let segment_position = SegmentPosition::try_from(flags)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let record_strand = Strand::from(flags);

    for (_, (_, strand)) in interval_tree.find(start..=end) {
        let feature_strand = match Strand::try_from(*strand) {
            Ok(s) => s,
            Err(_) => continue,
        };

        match (segment_position, record_strand, feature_strand) {
            (SegmentPosition::First, Strand::Forward, Strand::Forward)
            | (SegmentPosition::First, Strand::Reverse, Strand::Reverse)
            | (SegmentPosition::Last, Strand::Forward, Strand::Reverse)
            | (SegmentPosition::Last, Strand::Reverse, Strand::Forward) => {
                counts.forward += 1;
            }
            (SegmentPosition::First, Strand::Forward, Strand::Reverse)
            | (SegmentPosition::First, Strand::Reverse, Strand::Forward)
            | (SegmentPosition::Last, Strand::Forward, Strand::Forward)
            | (SegmentPosition::Last, Strand::Reverse, Strand::Reverse) => {
                counts.reverse += 1;
            }
        }

        counts.matches += 1;
    }

    Ok(())
}

fn count_single_end_record(
    counts: &mut Counts,
    interval_tree: &IntervalTree<Position, Entry>,
    flags: sam::alignment::record::Flags,
    start: Position,
    end: Position,
) {
    let record_strand = Strand::from(flags);

    for (_, (_, strand)) in interval_tree.find(start..=end) {
        let feature_strand = match Strand::try_from(*strand) {
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
}

pub fn detect_specification<P>(
    src: P,
    interval_trees: &IntervalTrees,
) -> io::Result<(LibraryLayout, StrandSpecification, f64)>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    reader.read_header()?;

    let mut record = bam::Record::default();
    let mut n = 0;
    let mut counts = Counts::default();

    while n < MAX_RECORDS && reader.read_record(&mut record)? != 0 {
        n += 1;

        let flags = record.flags();

        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        let reference_sequence_id = record
            .reference_sequence_id()
            .transpose()?
            .expect("missing reference sequence ID");

        let Some(interval_tree) = interval_trees.get(reference_sequence_id) else {
            continue;
        };

        let alignment_start = record.alignment_start().transpose()?;
        let alignment_end = record.alignment_end().transpose()?;

        let (Some(start), Some(end)) = (alignment_start, alignment_end) else {
            panic!("missing alignment context");
        };

        if flags.is_segmented() {
            counts.paired += 1;
            count_paired_end_record(&mut counts, interval_tree, flags, start, end)?;
        } else {
            count_single_end_record(&mut counts, interval_tree, flags, start, end);
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
