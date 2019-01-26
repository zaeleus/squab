pub use self::count::{Context, count_paired_end_records, count_single_end_records};
pub use self::record_pairs::{PairPosition, RecordPairs};
pub use self::strand::Strand;

pub mod count;
pub mod record_pairs;
pub mod strand;

use std::collections::{HashMap, HashSet};
use std::io;
use std::ops::Range;
use std::path::Path;

use csv::StringRecord;
use interval_tree::IntervalTree;
use log::info;
use noodles::formats::bam::{ByteRecord, Cigar, Flag, cigar};
use noodles::formats::gff;

const GFF_SEQ_NAME_INDEX: usize = 0;
const GFF_TYPE_INDEX: usize = 2;
const GFF_START_INDEX: usize = 3;
const GFF_END_INDEX: usize = 4;
const GFF_STRAND_INDEX: usize = 6;
const GFF_ATTRS_INDEX: usize = 8;

pub type Features = HashMap<String, IntervalTree<u64, Entry>>;

#[derive(Default)]
pub struct Entry(pub String, pub Strand);

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
        let record = result?;

        let ty = parse_type(&record)?;

        if ty != feature_type {
            continue;
        }

        let seq_name = parse_seq_name(&record)?;
        let start = parse_start(&record)?;
        let end = parse_end(&record)?;
        let strand = parse_strand(&record)?;
        let id = parse_attrs_and_get(&record, feature_id)?;

        let tree = features.entry(seq_name.to_string()).or_default();
        tree.insert(start..end, Entry(id.to_string(), strand));

        names.insert(id.to_string());
    }

    info!("read {} features", names.len());

    Ok((features, names))
}

fn parse_type(record: &StringRecord) -> io::Result<&str> {
    record.get(GFF_TYPE_INDEX).ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidData, "invalid type")
    })
}

fn parse_seq_name(record: &StringRecord) -> io::Result<&str> {
    record.get(GFF_SEQ_NAME_INDEX).ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidData, "invalid seq name")
    })
}

fn parse_start(record: &StringRecord) -> io::Result<u64> {
    record
        .get(GFF_START_INDEX)
        .and_then(|s| s.parse::<u64>().ok())
        .map(|s| s - 1)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid start"))
}

fn parse_end(record: &StringRecord) -> io::Result<u64> {
    record
        .get(GFF_END_INDEX)
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid end"))
}

fn parse_strand(record: &StringRecord) -> io::Result<Strand> {
    record
        .get(GFF_STRAND_INDEX)
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid strand"))
}

fn parse_attrs_and_get<'a>(record: &'a StringRecord, key: &str) -> io::Result<&'a str> {
    let attrs = record
        .get(GFF_ATTRS_INDEX)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid attrs"))?;

    for attr in attrs.split(';').map(str::trim_left) {
        if attr.starts_with(key) {
            return attr.split('"').nth(1).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("could not parse attribute '{}'", attr),
                )
            });
        }
    }

    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        format!("missing attribute '{}'", key),
    ))
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

#[cfg(test)]
mod tests {
    use csv::StringRecord;

    use crate::Strand;

    use super::*;

    fn build_record() -> StringRecord {
        StringRecord::from(vec![
           "chr1",
           "HAVANA",
           "gene",
           "11869",
           "14409",
           ".",
           "+",
           ".",
           r#"gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";"#
        ])
    }

    #[test]
    pub fn test_parse_type() {
        let record = build_record();
        assert_eq!(parse_type(&record).unwrap(), "gene");
    }

    #[test]
    pub fn test_parse_seq_name() {
        let record = build_record();
        assert_eq!(parse_seq_name(&record).unwrap(), "chr1");
    }

    #[test]
    pub fn test_parse_start() {
        let record = build_record();
        assert_eq!(parse_start(&record).unwrap(), 11868);
    }

    #[test]
    pub fn test_parse_end() {
        let record = build_record();
        assert_eq!(parse_end(&record).unwrap(), 14409);
    }

    #[test]
    pub fn test_parse_strand() {
        let record = build_record();
        assert_eq!(parse_strand(&record).unwrap(), Strand::Forward);
    }

    #[test]
    fn test_parse_attrs_and_get() {
        let record = build_record();

        assert_eq!(parse_attrs_and_get(&record, "gene_id").unwrap(), "ENSG00000223972.5");
        assert_eq!(parse_attrs_and_get(&record, "gene_type").unwrap(), "transcribed_unprocessed_pseudogene");
        assert_eq!(parse_attrs_and_get(&record, "gene_name").unwrap(), "DDX11L1");
        assert_eq!(parse_attrs_and_get(&record, "havana_gene").unwrap(), "OTTHUMG00000000961.2");

        assert!(parse_attrs_and_get(&record, "level").is_err());
        assert!(parse_attrs_and_get(&record, "dne").is_err());
    }
}
