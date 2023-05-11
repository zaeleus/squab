pub use self::{
    commands::StrandSpecificationOption,
    count::{count_paired_end_records, count_single_end_records, Context},
    feature::Feature,
    match_intervals::MatchIntervals,
    record::Record,
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
pub mod record;
pub mod record_pairs;

pub use self::cli::Cli;

use std::{
    collections::{HashMap, HashSet},
    io::{self, BufRead},
};

use interval_tree::IntervalTree;
use noodles::core::Position;
use tracing::info;

pub type Entry = (String, noodles::gff::record::Strand);
pub type Features = HashMap<String, IntervalTree<Position, Entry>>;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum StrandSpecification {
    None,
    Forward,
    Reverse,
}

pub fn read_features<R>(
    reader: &mut noodles::gff::Reader<R>,
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

        let feature = Feature::new(reference_sequence_name.into(), start, end, strand);

        let id = record
            .attributes()
            .iter()
            .find(|e| e.key() == feature_id)
            .map(|e| e.value())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("missing attribute '{feature_id}'"),
                )
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
