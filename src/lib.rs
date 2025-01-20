pub use self::{
    commands::StrandSpecificationOption,
    count::{count_paired_end_records, count_single_end_records, Context},
    feature::Feature,
    match_intervals::MatchIntervals,
    record_pairs::{RecordPairs, SegmentPosition},
};

pub mod cli;
pub mod commands;
pub mod count;
pub mod detect;
pub mod feature;
mod gff;
mod match_intervals;
pub mod normalization;
pub mod record_pairs;

pub use self::cli::Cli;

use std::{
    collections::{HashMap, HashSet},
    fmt,
    io::{self, BufRead},
};

use bstr::BString;
use interval_tree::IntervalTree;
use noodles::core::{self as core, Position};
use thiserror::Error;
use tracing::info;

pub type Entry = (String, noodles::gff::record::Strand);
pub type Features = HashMap<BString, IntervalTree<Position, Entry>>;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum StrandSpecification {
    None,
    Forward,
    Reverse,
}

impl fmt::Display for StrandSpecification {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::None => "none".fmt(f),
            Self::Forward => "forward".fmt(f),
            Self::Reverse => "reverse".fmt(f),
        }
    }
}

#[derive(Debug, Error)]
pub enum ReadFeaturesError {
    #[error("invalid position")]
    InvalidPosition(#[from] core::position::ParseError),
    #[error("missing attribute: {0}")]
    MissingAttribute(String),
    #[error("invalid attribute: {0}")]
    InvalidAttribute(String),
    #[error("I/O error")]
    Io(#[from] io::Error),
}

pub fn read_features<R>(
    reader: &mut noodles::gff::Reader<R>,
    feature_type: &str,
    feature_id: &str,
) -> Result<HashMap<String, Vec<Feature>>, ReadFeaturesError>
where
    R: BufRead,
{
    use noodles::gff::{record::attributes::field::Value, Line};

    let mut features: HashMap<String, Vec<Feature>> = HashMap::new();

    info!("reading features");

    let mut line = Line::default();

    while reader.read_line(&mut line)? != 0 {
        let Some(record) = line.as_record().transpose()? else {
            continue;
        };

        if record.ty() != feature_type {
            continue;
        }

        let reference_sequence_name = record.reference_sequence_name();
        let start = record.start()?;
        let end = record.end()?;
        let strand = record.strand()?;

        let feature = Feature::new(reference_sequence_name.into(), start, end, strand);

        let attributes = record.attributes();
        let id = attributes
            .get(feature_id)
            .ok_or_else(|| ReadFeaturesError::MissingAttribute(feature_id.into()))?
            .map_err(|_| ReadFeaturesError::InvalidAttribute(feature_id.into()))
            .and_then(|value| match value {
                Value::String(s) => Ok(s),
                Value::Array(_) => Err(ReadFeaturesError::InvalidAttribute(feature_id.into())),
            })?;

        let list = features.entry(id.into()).or_default();

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
            let reference_sequence_name = feature.reference_sequence_name();

            let start = feature.start();
            let end = feature.end();

            let strand = feature.strand();

            let tree = interval_trees
                .entry(reference_sequence_name.into())
                .or_insert_with(IntervalTree::new);

            tree.insert(start..=end, (id.into(), strand));
        }

        names.insert(id.into());
    }

    (interval_trees, names)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_features() -> anyhow::Result<()> {
        use noodles::gff::record::Strand;

        let data = b"##gff-version 3
sq0\t.\texon\t1\t10\t.\t+\t.\tID=exon0;gene_id=gene0;gene_name=NDLS_gene0
sq0\t.\texon\t21\t30\t.\t+\t.\tID=exon1;gene_id=gene0;gene_name=NDLS_gene0
sq1\t.\texon\t41\t50\t.\t-\t.\tID=exon3;gene_id=gene1;gene_name=NDLS_gene1
";
        let mut reader = noodles::gff::Reader::new(&data[..]);

        let features = read_features(&mut reader, "exon", "gene_id")?;

        assert_eq!(features.len(), 2);
        assert_eq!(
            features["gene0"],
            [
                Feature::new(
                    String::from("sq0"),
                    Position::try_from(1)?,
                    Position::try_from(10)?,
                    Strand::Forward
                ),
                Feature::new(
                    String::from("sq0"),
                    Position::try_from(21)?,
                    Position::try_from(30)?,
                    Strand::Forward
                ),
            ]
        );
        assert_eq!(
            features["gene1"],
            [Feature::new(
                String::from("sq1"),
                Position::try_from(41)?,
                Position::try_from(50)?,
                Strand::Reverse
            ),]
        );

        Ok(())
    }
}
