pub mod fpkm;
mod method;
pub mod tpm;

pub use self::method::Method;

use bstr::{BString, ByteSlice};
use thiserror::Error;

use crate::{Feature, Features};

#[derive(Debug, Error)]
pub enum Error {
    #[error("missing feature: {0}")]
    MissingFeature(BString),
}

pub fn calculate_feature_lengths(
    features: &Features,
    names: &[BString],
) -> Result<Vec<usize>, Error> {
    names
        .iter()
        .map(|name| {
            features
                .get(name)
                .map(|segments| sum_nonoverlapping_feature_lengths(segments))
                .ok_or_else(|| Error::MissingFeature(name.as_bstr().into()))
        })
        .collect()
}

fn sum_nonoverlapping_feature_lengths(features: &[Feature]) -> usize {
    merge_features(features).iter().map(|f| f.length()).sum()
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
    use noodles::{core::Position, gff::feature::record::Strand};

    use super::*;
    use crate::Feature;

    const STRAND: Strand = Strand::Forward;

    macro_rules! position {
        ($n:literal) => {
            const { Position::new($n).unwrap() }
        };
    }

    fn build_features() -> [Feature; 6] {
        //          111111111122222222222
        // 123456789012345678901234567890
        //  [--]   [--]   [----]
        //   [][-]  [----]

        [
            Feature::new(0, position!(2), position!(5), STRAND),
            Feature::new(0, position!(3), position!(4), STRAND),
            Feature::new(0, position!(5), position!(7), STRAND),
            Feature::new(0, position!(9), position!(12), STRAND),
            Feature::new(0, position!(10), position!(15), STRAND),
            Feature::new(0, position!(16), position!(21), STRAND),
        ]
    }

    #[test]
    fn test_sum_nonoverlapping_feature_lengths() {
        let features = build_features();
        let len = sum_nonoverlapping_feature_lengths(&features);
        assert_eq!(len, 19);
    }

    #[test]
    fn test_merge_features() {
        let features = build_features();
        let actual = merge_features(&features);

        let expected = [
            Feature::new(0, position!(2), position!(7), STRAND),
            Feature::new(0, position!(9), position!(15), STRAND),
            Feature::new(0, position!(16), position!(21), STRAND),
        ];

        assert_eq!(actual, expected);
    }
}
