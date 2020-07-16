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
    io::{self, BufRead},
    ops::RangeInclusive,
    str::FromStr,
};

use interval_tree::IntervalTree;
use log::info;
use noodles_bam::{self as bam, record::cigar};
use noodles_sam as sam;

pub type Features = HashMap<String, IntervalTree<u64, Entry>>;

#[derive(Default)]
pub struct Entry(pub String, pub noodles_gff::record::Strand);

pub fn read_features<R>(
    reader: &mut noodles_gff::Reader<R>,
    feature_type: &str,
    feature_id: &str,
) -> io::Result<HashMap<String, Vec<Feature>>>
where
    R: BufRead,
{
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
            let reference_sequence_name = feature.reference_sequence_name();

            let start = feature.start();
            let end = feature.end();

            let strand = feature.strand();

            let tree = interval_trees
                .entry(reference_sequence_name.into())
                .or_default();
            tree.insert(start..=end, Entry(id.into(), strand));
        }

        names.insert(id.into());
    }

    (interval_trees, names)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_features() -> io::Result<()> {
        use noodles_gff::record::Strand;

        let data = b"##gff-version 3
sq0\t.\texon\t1\t10\t.\t+\t.\tID=exon0;gene_id=gene0;gene_name=NDLS_gene0
sq0\t.\texon\t21\t30\t.\t+\t.\tID=exon1;gene_id=gene0;gene_name=NDLS_gene0
sq1\t.\texon\t41\t50\t.\t-\t.\tID=exon3;gene_id=gene1;gene_name=NDLS_gene1
";
        let mut reader = noodles_gff::Reader::new(&data[..]);

        let features = read_features(&mut reader, "exon", "gene_id")?;

        assert_eq!(features.len(), 2);
        assert_eq!(
            features["gene0"],
            [
                Feature::new(String::from("sq0"), 1, 10, Strand::Forward),
                Feature::new(String::from("sq0"), 21, 30, Strand::Forward),
            ]
        );
        assert_eq!(
            features["gene1"],
            [Feature::new(String::from("sq1"), 41, 50, Strand::Reverse),]
        );

        Ok(())
    }
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
                    let end = start + len - 1;
                    self.start += len;
                    return Some((start..=end, self.is_reverse));
                }
                Kind::Deletion | Kind::Skip => {
                    self.start += len;
                }
                _ => continue,
            }
        }
    }
}

#[cfg(test)]
mod test_cigar_to_intervals {
    use noodles_bam::{self as bam, record::cigar};
    use noodles_sam::{self as sam, record::cigar::op};

    use super::CigarToIntervals;

    fn build_raw_cigar() -> Vec<u8> {
        let ops = [
            u32::from(cigar::Op::new(op::Kind::Match, 1)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::Insertion, 2)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::Deletion, 3)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::Skip, 4)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::SoftClip, 5)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::HardClip, 6)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::Pad, 7)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::SeqMatch, 8)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::SeqMismatch, 9)).to_le_bytes(),
        ];

        ops.iter().flatten().copied().collect()
    }

    #[test]
    fn test_new() {
        use sam::record::Flags;

        let raw_cigar = build_raw_cigar();
        let cigar = bam::record::Cigar::new(&raw_cigar);

        let start = 1;

        let flags = Flags::PAIRED | Flags::PROPER_PAIR | Flags::MATE_REVERSE | Flags::READ_1;
        let it = CigarToIntervals::new(&cigar, start, flags, false);
        assert!(!it.is_reverse);

        let flags = Flags::PAIRED | Flags::PROPER_PAIR | Flags::MATE_REVERSE | Flags::READ_1;
        let it = CigarToIntervals::new(&cigar, start, flags, true);
        assert!(it.is_reverse);

        let flags = Flags::PAIRED | Flags::PROPER_PAIR | Flags::REVERSE | Flags::READ_2;
        let it = CigarToIntervals::new(&cigar, start, flags, false);
        assert!(it.is_reverse);

        let flags = Flags::PAIRED | Flags::PROPER_PAIR | Flags::REVERSE | Flags::READ_2;
        let it = CigarToIntervals::new(&cigar, start, flags, true);
        assert!(!it.is_reverse);
    }

    #[test]
    fn test_next() {
        use sam::record::Flags;

        let raw_cigar = build_raw_cigar();
        let cigar = bam::record::Cigar::new(&raw_cigar);

        let start = 1;
        let flags = Flags::PAIRED | Flags::PROPER_PAIR | Flags::MATE_REVERSE | Flags::READ_1;
        let mut it = CigarToIntervals::new(&cigar, start, flags, false);

        let (interval, is_reverse) = it.next().unwrap();
        assert_eq!(interval, 1..=1);
        assert!(!is_reverse);

        let (interval, is_reverse) = it.next().unwrap();
        assert_eq!(interval, 9..=16);
        assert!(!is_reverse);

        let (interval, is_reverse) = it.next().unwrap();
        assert_eq!(interval, 17..=25);
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
