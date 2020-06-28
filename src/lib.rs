pub use self::{
    count::{count_paired_end_records, count_single_end_records, Context},
    feature::Feature,
    record_pairs::{PairPosition, RecordPairs},
};

pub mod commands;
pub mod count;
pub mod detect;
pub mod feature;
mod gff;
pub mod normalization;
pub mod record_pairs;

use std::{
    collections::{HashMap, HashSet},
    hash::BuildHasher,
    io,
    ops::RangeInclusive,
    path::Path,
    str::FromStr,
};

use interval_tree::IntervalTree;
use log::info;
use noodles_bam::{self as bam, record::cigar};
use noodles_sam as sam;

pub type Features = HashMap<String, IntervalTree<u64, Entry>>;

#[derive(Default)]
pub struct Entry(pub String, pub noodles_gff::record::Strand);

pub fn read_features<P>(
    src: P,
    feature_type: &str,
    feature_id: &str,
) -> io::Result<HashMap<String, Vec<Feature>>>
where
    P: AsRef<Path>,
{
    let mut reader = gff::open(src)?;
    let mut features: HashMap<String, Vec<Feature>> = HashMap::new();

    info!("reading features");

    for result in reader.records() {
        let record = result?;

        let ty = record.ty();

        if ty != feature_type {
            continue;
        }

        let reference_sequence_name = record.reference_sequence_name();
        let start = record.start();
        let end = record.end();

        let strand = record.strand();

        let attributes = record.attributes();
        let id = attributes
            .iter()
            .find(|e| e.key() == feature_id)
            .map(|e| e.value())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("missing attribute '{}'", feature_id),
                )
            })?;

        let list = features.entry(id.into()).or_default();

        let feature = Feature::new(
            reference_sequence_name.into(),
            start as u64,
            end as u64,
            strand,
        );

        list.push(feature);
    }

    info!("read {} unique features", features.len());

    Ok(features)
}

pub fn build_interval_trees<S: BuildHasher>(
    feature_map: &HashMap<String, Vec<Feature>, S>,
) -> (Features, HashSet<String>) {
    let mut interval_trees = Features::new();
    let mut names = HashSet::new();

    for (id, features) in feature_map {
        for feature in features {
            let reference_name = feature.reference_name();

            let start = feature.start() - 1;
            let end = feature.end();

            let strand = feature.strand();

            let tree = interval_trees.entry(reference_name.into()).or_default();
            tree.insert(start..=end, Entry(id.into(), strand));
        }

        names.insert(id.into());
    }

    (interval_trees, names)
}

pub struct CigarToIntervals<'a> {
    ops: cigar::Ops<'a>,
    start: u64,
    is_reverse: bool,
}

impl<'a> CigarToIntervals<'a> {
    fn new(
        cigar: &'a bam::record::Cigar,
        start: u64,
        flag: sam::record::Flags,
        reverse: bool,
    ) -> CigarToIntervals<'a> {
        let ops = cigar.ops();

        let is_reverse = if reverse {
            !flag.is_reverse()
        } else {
            flag.is_reverse()
        };

        CigarToIntervals {
            ops,
            start,
            is_reverse,
        }
    }
}

impl<'a> Iterator for CigarToIntervals<'a> {
    type Item = (RangeInclusive<u64>, bool);

    fn next(&mut self) -> Option<Self::Item> {
        use sam::record::cigar::op::Kind;

        loop {
            let op = self.ops.next()?;

            let len = u64::from(op.len());

            match op.kind() {
                Kind::Match | Kind::SeqMatch | Kind::SeqMismatch => {
                    let start = self.start;
                    let end = start + len;
                    self.start += len;
                    return Some((start..=end, self.is_reverse));
                }
                Kind::Deletion | Kind::Skip => {}
                _ => continue,
            }

            self.start += len;
        }
    }
}

#[cfg(test)]
mod tests {
    use noodles_bam::{self as bam, record::cigar};
    use noodles_sam::{self as sam, record::cigar::op};

    use super::CigarToIntervals;

    fn build_raw_cigar() -> Vec<u8> {
        let ops = [
            u32::from(cigar::Op::new(op::Kind::Match, 1)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::Insertion, 2)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::Deletion, 3)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::Skip, 5)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::SoftClip, 8)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::HardClip, 13)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::Pad, 21)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::SeqMatch, 34)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::SeqMismatch, 55)).to_le_bytes(),
        ];

        ops.iter().flat_map(|u| u).cloned().collect()
    }

    #[test]
    fn test_new() {
        let raw_cigar = build_raw_cigar();
        let cigar = bam::record::Cigar::new(&raw_cigar);

        let start = 0;

        let flags = sam::record::Flags::from(99);
        let it = CigarToIntervals::new(&cigar, start, flags, false);
        assert!(!it.is_reverse);

        let flags = sam::record::Flags::from(99);
        let it = CigarToIntervals::new(&cigar, start, flags, true);
        assert!(it.is_reverse);

        let flags = sam::record::Flags::from(147);
        let it = CigarToIntervals::new(&cigar, start, flags, false);
        assert!(it.is_reverse);

        let flags = sam::record::Flags::from(147);
        let it = CigarToIntervals::new(&cigar, start, flags, true);
        assert!(!it.is_reverse);
    }

    #[test]
    fn test_next() {
        let raw_cigar = build_raw_cigar();
        let cigar = bam::record::Cigar::new(&raw_cigar);

        let flags = sam::record::Flags::from(99);
        let mut it = CigarToIntervals::new(&cigar, 0, flags, false);

        let (interval, is_reverse) = it.next().unwrap();
        assert_eq!(interval, 0..=1);
        assert!(!is_reverse);

        let (interval, is_reverse) = it.next().unwrap();
        assert_eq!(interval, 9..=43);
        assert!(!is_reverse);

        let (interval, is_reverse) = it.next().unwrap();
        assert_eq!(interval, 43..=98);
        assert!(!is_reverse);

        assert!(it.next().is_none());
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum StrandSpecificationOption {
    None,
    Forward,
    Reverse,
    Auto,
}

impl FromStr for StrandSpecificationOption {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "none" => Ok(Self::None),
            "forward" => Ok(Self::Forward),
            "reverse" => Ok(Self::Reverse),
            "auto" => Ok(Self::Auto),
            _ => Err(()),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum StrandSpecification {
    None,
    Forward,
    Reverse,
}
