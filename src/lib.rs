pub use self::{
    count::{count_paired_end_records, count_single_end_records, Context},
    record_pairs::{PairPosition, RecordPairs},
};

pub mod count;
pub mod record_pairs;

use std::{
    collections::{HashMap, HashSet},
    io,
    ops::Range,
    path::Path,
    str::FromStr,
};

use interval_tree::IntervalTree;
use log::info;
use noodles_bam::{cigar, Cigar, Flag};
use noodles_gff as gff;

pub type Features = HashMap<String, IntervalTree<u64, Entry>>;

#[derive(Default)]
pub struct Entry(pub String, pub gff::Strand);

pub fn read_features<P>(
    src: P,
    feature_type: &str,
    feature_id: &str,
) -> io::Result<(Features, HashSet<String>)>
where
    P: AsRef<Path>,
{
    let mut reader = gff::open(src)?;
    let mut features = Features::new();
    let mut names = HashSet::new();

    info!("reading features");

    for result in reader.records() {
        let row = result?;
        let record = gff::Record::new(row);

        let ty = record.feature().map_err(invalid_data)?;

        if ty != feature_type {
            continue;
        }

        let seq_name = record.seq_name().map_err(invalid_data)?;
        let start = record.start().map_err(invalid_data)? - 1;
        let end = record.end().map_err(invalid_data)?;

        let strand = record.strand().map_err(invalid_data)?;

        let attributes = record.attributes().map_err(invalid_data)?;
        let id = attributes.get(feature_id).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing attribute '{}'", feature_id),
            )
        })?;

        let tree = features.entry(seq_name.to_string()).or_default();
        tree.insert(start..end, Entry(id.to_string(), strand));

        names.insert(id.to_string());
    }

    info!("read {} features", names.len());

    Ok((features, names))
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
    fn new(cigar: &'a Cigar, start: u64, flag: Flag, reverse: bool) -> CigarToIntervals<'a> {
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
        loop {
            let op = self.ops.next()?;

            let len = u64::from(op.len());

            match op {
                cigar::Op::Match(_) | cigar::Op::SeqMatch(_) | cigar::Op::SeqMismatch(_) => {
                    let start = self.start;
                    let end = start + len;
                    self.start += len;
                    return Some((start..end, self.is_reverse));
                }
                cigar::Op::Deletion(_) | cigar::Op::Skip(_) => {}
                _ => continue,
            }

            self.start += len;
        }
    }
}

#[cfg(test)]
mod tests {
    use noodles_bam::{cigar, Cigar, Flag};

    use super::CigarToIntervals;

    fn build_raw_cigar() -> Vec<u8> {
        let ops = [
            u32::from(cigar::Op::Match(1)).to_le_bytes(),
            u32::from(cigar::Op::Insertion(2)).to_le_bytes(),
            u32::from(cigar::Op::Deletion(3)).to_le_bytes(),
            u32::from(cigar::Op::Skip(5)).to_le_bytes(),
            u32::from(cigar::Op::SoftClip(8)).to_le_bytes(),
            u32::from(cigar::Op::HardClip(13)).to_le_bytes(),
            u32::from(cigar::Op::Pad(21)).to_le_bytes(),
            u32::from(cigar::Op::SeqMatch(34)).to_le_bytes(),
            u32::from(cigar::Op::SeqMismatch(55)).to_le_bytes(),
        ];

        ops.iter().flat_map(|u| u).cloned().collect()
    }

    #[test]
    fn test_new() {
        let raw_cigar = build_raw_cigar();
        let cigar = Cigar::new(&raw_cigar);

        let start = 0;

        let flag = Flag::from(99);
        let it = CigarToIntervals::new(&cigar, start, flag, false);
        assert!(!it.is_reverse);

        let flag = Flag::from(99);
        let it = CigarToIntervals::new(&cigar, start, flag, true);
        assert!(it.is_reverse);

        let flag = Flag::from(147);
        let it = CigarToIntervals::new(&cigar, start, flag, false);
        assert!(it.is_reverse);

        let flag = Flag::from(147);
        let it = CigarToIntervals::new(&cigar, start, flag, true);
        assert!(!it.is_reverse);
    }

    #[test]
    fn test_next() {
        let raw_cigar = build_raw_cigar();
        let cigar = Cigar::new(&raw_cigar);

        let flag = Flag::from(99);
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

#[derive(Clone, Copy, Debug)]
pub enum StrandSpecification {
    None,
    Forward,
}

impl FromStr for StrandSpecification {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "none" => Ok(Self::None),
            "forward" => Ok(Self::Forward),
            _ => Err(()),
        }
    }
}
