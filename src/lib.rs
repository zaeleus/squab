pub use self::count::{Context, count_paired_end_records, count_single_end_records};
pub use self::record_pairs::{PairPosition, RecordPairs};

pub mod count;
pub mod record_pairs;

use std::collections::{HashMap, HashSet};
use std::io;
use std::ops::Range;
use std::path::Path;

use interval_tree::IntervalTree;
use log::info;
use noodles::formats::bam::{Cigar, Flag, cigar};
use noodles::formats::gff;

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
        let is_reverse = if reverse { !flag.is_reverse() } else { flag.is_reverse() };
        CigarToIntervals { ops, start, is_reverse, }
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
                    let end = self.start + len;
                    return Some((self.start..end, self.is_reverse));
                },
                cigar::Op::Deletion(_) | cigar::Op::Skip(_) => {},
                _ => continue,
            }

            self.start += len;
        }
    }
}
