pub use self::{
    count::{count_paired_end_records, count_single_end_records, Context},
    feature::Feature,
    record_pairs::{PairPosition, RecordPairs},
};

pub mod count;
pub mod detect;
pub mod feature;
pub mod normalization;
pub mod record_pairs;
pub mod writer;

use std::{
    collections::{HashMap, HashSet},
    convert::TryFrom,
    io,
    ops::Range,
    path::Path,
    str::FromStr,
};

use interval_tree::IntervalTree;
use log::info;
use noodles_bam::{self as bam, cigar};
use noodles_gff as gff;

pub type Features = HashMap<String, IntervalTree<u64, Entry>>;

#[derive(Default)]
pub struct Entry(pub String, pub gff::Strand);

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
        let row = result?;
        let record = gff::Record::new(row);

        let ty = record.feature().map_err(invalid_data)?;

        if ty != feature_type {
            continue;
        }

        let seq_name = record.seq_name().map_err(invalid_data)?;
        let start = record.start().map_err(invalid_data)?;
        let end = record.end().map_err(invalid_data)?;

        let strand = record.strand().map_err(invalid_data)?;

        let attributes = record.attributes().map_err(invalid_data)?;
        let id = attributes.get(feature_id).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing attribute '{}'", feature_id),
            )
        })?;

        let list = features.entry(id.into()).or_default();
        let feature = Feature::new(seq_name.into(), start, end, strand);
        list.push(feature);
    }

    info!("read {} unique features", features.len());

    Ok(features)
}

pub fn build_interval_trees(
    feature_map: &HashMap<String, Vec<Feature>>,
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
            tree.insert(start..end, Entry(id.into(), strand));
        }

        names.insert(id.into());
    }

    (interval_trees, names)
}

fn invalid_data(e: gff::record::Error) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, e)
}

pub struct CigarToIntervals<'a> {
    ops: cigar::Ops<'a>,
    start: u64,
    is_reverse: bool,
}

impl<'a> CigarToIntervals<'a> {
    fn new(
        cigar: &'a bam::Cigar,
        start: u64,
        flag: bam::Flag,
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
    type Item = (Range<u64>, bool);

    fn next(&mut self) -> Option<Self::Item> {
        use bam::cigar::op::Kind;

        loop {
            let op = self.ops.next()?;

            let len = u64::from(op.len());

            match op.kind() {
                Kind::Match | Kind::SeqMatch | Kind::SeqMismatch => {
                    let start = self.start;
                    let end = start + len;
                    self.start += len;
                    return Some((start..end, self.is_reverse));
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
    use noodles_bam::{
        self as bam,
        cigar::{self, op},
    };

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
        let cigar = bam::Cigar::new(&raw_cigar);

        let start = 0;

        let flag = bam::Flag::from(99);
        let it = CigarToIntervals::new(&cigar, start, flag, false);
        assert!(!it.is_reverse);

        let flag = bam::Flag::from(99);
        let it = CigarToIntervals::new(&cigar, start, flag, true);
        assert!(it.is_reverse);

        let flag = bam::Flag::from(147);
        let it = CigarToIntervals::new(&cigar, start, flag, false);
        assert!(it.is_reverse);

        let flag = bam::Flag::from(147);
        let it = CigarToIntervals::new(&cigar, start, flag, true);
        assert!(!it.is_reverse);
    }

    #[test]
    fn test_next() {
        let raw_cigar = build_raw_cigar();
        let cigar = bam::Cigar::new(&raw_cigar);

        let flag = bam::Flag::from(99);
        let mut it = CigarToIntervals::new(&cigar, 0, flag, false);

        let (interval, is_reverse) = it.next().unwrap();
        assert_eq!(interval, 0..1);
        assert!(!is_reverse);

        let (interval, is_reverse) = it.next().unwrap();
        assert_eq!(interval, 9..43);
        assert!(!is_reverse);

        let (interval, is_reverse) = it.next().unwrap();
        assert_eq!(interval, 43..98);
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

impl TryFrom<StrandSpecificationOption> for StrandSpecification {
    type Error = ();

    fn try_from(option: StrandSpecificationOption) -> Result<Self, Self::Error> {
        match option {
            StrandSpecificationOption::None => Ok(Self::None),
            StrandSpecificationOption::Forward => Ok(Self::Forward),
            StrandSpecificationOption::Reverse => Ok(Self::Reverse),
            _ => Err(()),
        }
    }
}
