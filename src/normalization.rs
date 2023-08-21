mod fpkm;
mod method;
mod tpm;
mod writer;

pub use self::{fpkm::calculate_fpkms, method::Method, tpm::calculate_tpms, writer::Writer};

use std::collections::HashMap;

use thiserror::Error;

use crate::Feature;

type FeatureMap = HashMap<String, Vec<Feature>>;

#[derive(Debug, Error)]
pub enum Error {
    #[error("missing feature: {0}")]
    MissingFeature(String),
}

fn sum_nonoverlapping_feature_lengths(features: &[Feature]) -> usize {
    merge_features(features).iter().map(|f| f.len()).sum()
}

/// Merges a list of overlapping features into a list of non-overlapping intervals.
///
/// The intervals are assumed to be inclusive.
///
/// # Panics
///
/// Panics when the input list of features is empty.
///
/// # Example
///
/// Given the list of intervals [2, 5], [3, 4], [5, 7], [9, 12], [10, 15], and
/// [16, 21], the following overlap and are combined into single features.
///
///   * [2, 5], [3, 4], [5, 7] => [2, 7]
///   * [9, 12], [10, 15] => [9, 15]
///   * [16, 21] => [16, 21]
fn merge_features(features: &[Feature]) -> Vec<Feature> {
    assert!(!features.is_empty());

    let mut features = features.to_vec();
    features.sort_unstable_by_key(|i| i.start());

    let mut merged_features = Vec::with_capacity(features.len());
    merged_features.push(features[0].clone());

    for b in features {
        let a = merged_features.last_mut().expect("list cannot be empty");

        if b.start() > a.end() {
            merged_features.push(b);
            continue;
        }

        if a.end() < b.end() {
            *a.end_mut() = b.end();
        }
    }

    merged_features
}

#[cfg(test)]
mod tests {
    use crate::Feature;

    use noodles::{core::Position, gff};

    use super::*;

    #[test]
    fn test_sum_nonoverlapping_feature_lengths(
    ) -> Result<(), noodles::core::position::TryFromIntError> {
        let reference_name = String::from("chr1");
        let strand = gff::record::Strand::Forward;

        let features = [
            Feature::new(
                reference_name.clone(),
                Position::try_from(2)?,
                Position::try_from(5)?,
                strand,
            ),
            Feature::new(
                reference_name.clone(),
                Position::try_from(3)?,
                Position::try_from(4)?,
                strand,
            ),
            Feature::new(
                reference_name.clone(),
                Position::try_from(5)?,
                Position::try_from(7)?,
                strand,
            ),
            Feature::new(
                reference_name.clone(),
                Position::try_from(9)?,
                Position::try_from(12)?,
                strand,
            ),
            Feature::new(
                reference_name.clone(),
                Position::try_from(10)?,
                Position::try_from(15)?,
                strand,
            ),
            Feature::new(
                reference_name,
                Position::try_from(16)?,
                Position::try_from(21)?,
                strand,
            ),
        ];

        let len = sum_nonoverlapping_feature_lengths(&features);
        assert_eq!(len, 19);

        Ok(())
    }

    #[test]
    fn test_merge_features() -> Result<(), noodles::core::position::TryFromIntError> {
        let reference_name = String::from("chr1");
        let strand = gff::record::Strand::Forward;

        let features = [
            Feature::new(
                reference_name.clone(),
                Position::try_from(2)?,
                Position::try_from(5)?,
                strand,
            ),
            Feature::new(
                reference_name.clone(),
                Position::try_from(3)?,
                Position::try_from(4)?,
                strand,
            ),
            Feature::new(
                reference_name.clone(),
                Position::try_from(5)?,
                Position::try_from(7)?,
                strand,
            ),
            Feature::new(
                reference_name.clone(),
                Position::try_from(9)?,
                Position::try_from(12)?,
                strand,
            ),
            Feature::new(
                reference_name.clone(),
                Position::try_from(10)?,
                Position::try_from(15)?,
                strand,
            ),
            Feature::new(
                reference_name.clone(),
                Position::try_from(16)?,
                Position::try_from(21)?,
                strand,
            ),
        ];

        let actual = merge_features(&features);
        let expected = [
            Feature::new(
                reference_name.clone(),
                Position::try_from(2)?,
                Position::try_from(7)?,
                strand,
            ),
            Feature::new(
                reference_name.clone(),
                Position::try_from(9)?,
                Position::try_from(15)?,
                strand,
            ),
            Feature::new(
                reference_name,
                Position::try_from(16)?,
                Position::try_from(21)?,
                strand,
            ),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
