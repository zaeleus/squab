pub mod cli;
pub mod collections;
pub mod commands;
pub mod count;
pub mod counts;
pub mod detect;
pub mod feature;
mod gff;
mod match_intervals;
pub mod normalization;
pub mod record_pairs;

use std::{
    collections::HashMap,
    fmt,
    io::{self, BufRead},
};

use bstr::{BStr, BString};
use indexmap::IndexSet;
use noodles::{
    core::{self as core, Position},
    sam,
};
use thiserror::Error;

use self::collections::IntervalTree;
pub use self::{
    cli::Cli,
    commands::StrandSpecificationOption,
    count::{Context, count_paired_end_records, count_single_end_records},
    feature::Feature,
    match_intervals::MatchIntervals,
    record_pairs::{RecordPairs, SegmentPosition},
};

pub type ReferenceSequenceNames = IndexSet<BString>;
pub type Features = HashMap<BString, Vec<Feature>>;

pub type Entry<'f> = (&'f BStr, noodles::gff::feature::record::Strand);
pub type IntervalTrees<'f> = Vec<IntervalTree<Position, Entry<'f>>>;

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
    MissingAttribute(BString),
    #[error("invalid attribute: {0}")]
    InvalidAttribute(BString),
    #[error("I/O error")]
    Io(#[from] io::Error),
}

pub fn read_features<R>(
    reader: &mut noodles::gff::Reader<R>,
    feature_type: &str,
    feature_id: &str,
) -> Result<(ReferenceSequenceNames, Features), ReadFeaturesError>
where
    R: BufRead,
{
    use noodles::gff::{Line, record::attributes::field::Value};

    let mut reference_sequence_names = ReferenceSequenceNames::new();
    let mut features = Features::new();

    let mut line = Line::default();

    while reader.read_line(&mut line)? != 0 {
        let Some(record) = line.as_record().transpose()? else {
            continue;
        };

        if record.ty() != feature_type {
            continue;
        }

        let reference_sequence_name = record.reference_sequence_name();

        let reference_sequence_id = match reference_sequence_names
            .get_index_of(reference_sequence_name)
        {
            Some(id) => id,
            None => {
                let (id, _) = reference_sequence_names.insert_full(reference_sequence_name.into());
                id
            }
        };

        let start = record.start()?;
        let end = record.end()?;
        let strand = record.strand()?;

        let feature = Feature::new(reference_sequence_id, start, end, strand);

        let attributes = record.attributes();
        let id = attributes
            .get(feature_id.as_ref())
            .ok_or_else(|| ReadFeaturesError::MissingAttribute(feature_id.into()))?
            .map_err(|_| ReadFeaturesError::InvalidAttribute(feature_id.into()))
            .and_then(|value| match value {
                Value::String(s) => Ok(s),
                Value::Array(_) => Err(ReadFeaturesError::InvalidAttribute(feature_id.into())),
            })?;

        let segments = features.entry(id.into_owned()).or_default();
        segments.push(feature);
    }

    Ok((reference_sequence_names, features))
}

pub fn build_interval_trees<'f>(
    header: &sam::Header,
    reference_sequence_names: &ReferenceSequenceNames,
    features: &'f Features,
) -> IntervalTrees<'f> {
    let reference_sequences = header.reference_sequences();

    let mut raw_entries = vec![Vec::new(); reference_sequences.len()];

    for (name, segments) in features {
        for feature in segments {
            let reference_sequence_id = feature.reference_sequence_id();
            let start = feature.start();
            let end = feature.end();
            let strand = feature.strand();

            let reference_sequence_name = reference_sequence_names
                .get_index(reference_sequence_id)
                .unwrap();

            let Some(i) = reference_sequences.get_index_of(reference_sequence_name) else {
                continue;
            };

            // SAFETY: `intervals_trees.len() == reference_sequences.len()`
            let entries = &mut raw_entries[i];
            entries.push((start..=end, (name.as_ref(), strand)));
        }
    }

    raw_entries
        .into_iter()
        .map(|entries| entries.into_iter().collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_features() -> anyhow::Result<()> {
        use noodles::gff::feature::record::Strand;

        let data = b"##gff-version 3
sq0\t.\texon\t1\t10\t.\t+\t.\tID=exon0;gene_id=gene0;gene_name=NDLS_gene0
sq0\t.\texon\t21\t30\t.\t+\t.\tID=exon1;gene_id=gene0;gene_name=NDLS_gene0
sq1\t.\texon\t41\t50\t.\t-\t.\tID=exon3;gene_id=gene1;gene_name=NDLS_gene1
";
        let mut reader = noodles::gff::Reader::new(&data[..]);

        let (reference_sequence_names, features) = read_features(&mut reader, "exon", "gene_id")?;

        assert_eq!(
            reference_sequence_names,
            [BString::from("sq0"), BString::from("sq1")]
                .into_iter()
                .collect::<IndexSet<_>>()
        );

        assert_eq!(features.len(), 2);
        assert_eq!(
            features[BStr::new("gene0")],
            [
                Feature::new(
                    0,
                    Position::try_from(1)?,
                    Position::try_from(10)?,
                    Strand::Forward
                ),
                Feature::new(
                    0,
                    Position::try_from(21)?,
                    Position::try_from(30)?,
                    Strand::Forward
                ),
            ]
        );
        assert_eq!(
            features[BStr::new("gene1")],
            [Feature::new(
                1,
                Position::try_from(41)?,
                Position::try_from(50)?,
                Strand::Reverse
            ),]
        );

        Ok(())
    }
}
