extern crate interval_tree;
#[macro_use] extern crate log;
extern crate noodles;

pub use self::record_pairs::RecordPairs;
pub use self::strand::Strand;

pub mod record_pairs;
pub mod strand;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io;
use std::ops::Range;
use std::path::Path;

use interval_tree::IntervalTree;
use noodles::formats::bam::{ByteRecord, Cigar, Flag, cigar};
use noodles::formats::gff;

pub type Features = HashMap<String, IntervalTree<u64, Entry>>;

#[derive(Default)]
pub struct Entry(pub String, pub Strand);

pub fn read_features<P>(
    src: P,
    feature_type: &str,
    id: &str,
) -> io::Result<(Features, HashSet<String>)>
where
    P: AsRef<Path>,
{
    let mut reader = gff::Reader::<File>::open(src)?;
    let mut features: Features = HashMap::new();
    let mut names = HashSet::new();

    info!("reading features");

    for record in reader.records().filter_map(Result::ok) {
        let ty = &record[2];

        if ty != feature_type {
            continue;
        }

        let seq_name = &record[0];

        let start = {
            let s = record[3].parse::<u64>().map_err(|e| {
                io::Error::new(io::ErrorKind::InvalidInput, format!("{}", e))
            })?;

            s - 1
        };

        let end = record[4].parse::<u64>().map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidInput, format!("{}", e))
        })?;

        let strand = record[6].parse::<Strand>().map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidInput, e)
        })?;

        let attrs = &record[8];

        for attr in attrs.split(';').map(str::trim_left) {
            if attr.starts_with(id) {
                let gene_name = attr.split('"').nth(1).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("could not parse attribute '{}'", attr),
                    )
                })?;

                let tree = features.entry(seq_name.to_string()).or_default();
                tree.insert(start..end, Entry(gene_name.to_string(), strand));

                names.insert(gene_name.to_string());
            }
        }
    }

    info!("read {} features", names.len());

    Ok((features, names))
}

pub fn cigar_to_intervals(record: &ByteRecord, reverse: bool) -> Vec<(Range<u64>, bool)> {
    let mut start = record.pos() as u64;
    let cigar = Cigar::from_bytes(record.cigar());

    let flag = Flag::new(record.flag());
    let is_reverse = if reverse { !flag.is_reverse() } else { flag.is_reverse() };

    let mut intervals = Vec::with_capacity(cigar.len());

    for op in cigar.ops() {
        let len = op.len() as u64;

        match op {
            cigar::Op::Match(_) | cigar::Op::SeqMatch(_) | cigar::Op::SeqMismatch(_) => {
                let end = start + len;
                intervals.push((start..end, is_reverse));
            },
            cigar::Op::Deletion(_) | cigar::Op::Skip(_) => {},
            _ => continue,
        }

        start += len;
    }

    intervals
}
