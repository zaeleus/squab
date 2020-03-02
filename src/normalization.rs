use std::{collections::HashMap, error, fmt, str::FromStr};

use crate::Feature;

type Counts = HashMap<String, u64>;
type FeatureMap = HashMap<String, Vec<Feature>>;

/// Normalization method
#[derive(Debug)]
pub enum Method {
    /// fragments per kilobase per million mapped reads
    Fpkm,
    /// transcripts per million
    Tpm,
}

#[derive(Debug, Eq, PartialEq)]
pub struct ParseMethodError(String);

impl fmt::Display for ParseMethodError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid method: {}", self.0)
    }
}

impl error::Error for ParseMethodError {}

impl FromStr for Method {
    type Err = ParseMethodError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "fpkm" => Ok(Self::Fpkm),
            "tpm" => Ok(Self::Tpm),
            _ => Err(ParseMethodError(s.into())),
        }
    }
}

#[derive(Debug)]
pub enum Error {
    MissingFeature(String),
}

pub fn calculate_tpms(
    counts: &Counts,
    feature_map: &FeatureMap,
) -> Result<HashMap<String, f64>, Error> {
    let cpbs: HashMap<String, f64> = counts
        .iter()
        .map(|(name, &count)| {
            feature_map
                .get(name)
                .map(|features| {
                    let len = sum_nonoverlapping_feature_lengths(features);
                    let cpb = count as f64 / len as f64;
                    (name.clone(), cpb)
                })
                .ok_or_else(|| Error::MissingFeature(name.clone()))
        })
        .collect::<Result<_, _>>()?;

    let cpbs_sum = cpbs.values().sum();

    let tpms = cpbs
        .iter()
        .map(|(name, &cpb)| (name.clone(), calculate_tpm(cpb, cpbs_sum)))
        .collect();

    Ok(tpms)
}

pub fn calculate_fpkms(
    counts: &Counts,
    feature_map: &FeatureMap,
) -> Result<HashMap<String, f64>, Error> {
    let counts_sum = sum_counts(counts);

    counts
        .iter()
        .map(|(name, &count)| {
            feature_map
                .get(name)
                .map(|features| {
                    let len = sum_nonoverlapping_feature_lengths(features);
                    let fpkm = calculate_fpkm(count, len, counts_sum);
                    (name.clone(), fpkm)
                })
                .ok_or_else(|| Error::MissingFeature(name.clone()))
        })
        .collect()
}

fn calculate_tpm(cpb: f64, cpbs_sum: f64) -> f64 {
    cpb * 1e6 / cpbs_sum
}

fn sum_nonoverlapping_feature_lengths(features: &[Feature]) -> u64 {
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

fn sum_counts(counts: &Counts) -> u64 {
    counts.values().sum()
}

fn calculate_fpkm(count: u64, len: u64, counts_sum: u64) -> f64 {
    (count as f64 * 1e9) / (len as f64 * counts_sum as f64)
}

#[cfg(test)]
mod tests {
    use std::f64::EPSILON;

    use crate::feature::Feature;

    use noodles_gff as gff;

    use super::*;

    fn build_counts() -> HashMap<String, u64> {
        let counts = [
            (String::from("AAAS"), 645),
            (String::from("AC009952.3"), 1),
            (String::from("RPL37AP1"), 5714),
        ];

        counts.iter().cloned().collect()
    }

    fn build_feature_map() -> FeatureMap {
        let reference_name = String::from("chr1");
        let strand = gff::Strand::Forward;

        let features = [
            (
                String::from("AAAS"),
                vec![Feature::new(
                    reference_name.clone(),
                    53307456,
                    53324864,
                    strand,
                )],
            ),
            (
                String::from("AC009952.3"),
                vec![Feature::new(
                    reference_name.clone(),
                    9189629,
                    9204611,
                    strand,
                )],
            ),
            (
                String::from("RPL37AP1"),
                vec![Feature::new(
                    reference_name.clone(),
                    44466564,
                    44466842,
                    strand,
                )],
            ),
        ];

        features.iter().cloned().collect()
    }

    #[test]
    fn test_calculate_fpkms() {
        let counts = build_counts();
        let feature_map = build_feature_map();

        let fpkms = calculate_fpkms(&counts, &feature_map).unwrap();

        assert_eq!(fpkms.len(), 3);

        let a = fpkms["AAAS"];
        let b = 5825.440538780093;
        assert!((a - b).abs() < EPSILON);

        let a = fpkms["AC009952.3"];
        let b = 10.494073576888189;
        assert!((a - b).abs() < EPSILON);

        let a = fpkms["RPL37AP1"];
        let b = 3220170.8708099457;
        assert!((a - b).abs() < EPSILON);
    }

    #[test]
    fn test_calculate_fpkms_with_missing_feature() {
        let counts = build_counts();

        let mut feature_map = build_feature_map();
        feature_map.remove("AC009952.3");

        assert!(calculate_fpkms(&counts, &feature_map).is_err());
    }

    #[test]
    fn test_calculate_tpms() {
        let counts = build_counts();
        let feature_map = build_feature_map();

        let tpms = calculate_tpms(&counts, &feature_map).unwrap();

        assert_eq!(tpms.len(), 3);

        let a = tpms["AAAS"];
        let b = 1805.7744109493626;
        assert!((a - b).abs() < EPSILON);

        let a = tpms["AC009952.3"];
        let b = 3.252960768479983;
        assert!((a - b).abs() < EPSILON);

        let a = tpms["RPL37AP1"];
        let b = 998190.972628282;
        assert!((a - b).abs() < EPSILON);
    }

    #[test]
    fn test_calculate_tpms_with_missing_feature() {
        let counts = build_counts();

        let mut feature_map = build_feature_map();
        feature_map.remove("AC009952.3");

        assert!(calculate_tpms(&counts, &feature_map).is_err());
    }

    #[test]
    fn test_sum_nonoverlapping_feature_lengths() {
        let reference_name = String::from("chr1");
        let strand = gff::Strand::Forward;

        let features = [
            Feature::new(reference_name.clone(), 2, 5, strand),
            Feature::new(reference_name.clone(), 3, 4, strand),
            Feature::new(reference_name.clone(), 5, 7, strand),
            Feature::new(reference_name.clone(), 9, 12, strand),
            Feature::new(reference_name.clone(), 10, 15, strand),
            Feature::new(reference_name.clone(), 16, 21, strand),
        ];

        let len = sum_nonoverlapping_feature_lengths(&features);
        assert_eq!(len, 19);
    }

    #[test]
    fn test_merge_features() {
        let reference_name = String::from("chr1");
        let strand = gff::Strand::Forward;

        let features = [
            Feature::new(reference_name.clone(), 2, 5, strand),
            Feature::new(reference_name.clone(), 3, 4, strand),
            Feature::new(reference_name.clone(), 5, 7, strand),
            Feature::new(reference_name.clone(), 9, 12, strand),
            Feature::new(reference_name.clone(), 10, 15, strand),
            Feature::new(reference_name.clone(), 16, 21, strand),
        ];

        let actual = merge_features(&features);
        let expected = [
            Feature::new(reference_name.clone(), 2, 7, strand),
            Feature::new(reference_name.clone(), 9, 15, strand),
            Feature::new(reference_name.clone(), 16, 21, strand),
        ];

        assert_eq!(actual, expected);
    }
}
