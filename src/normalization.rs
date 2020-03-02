use std::{collections::HashMap, error, fmt, str::FromStr};

use crate::Feature;

type Counts = HashMap<String, u64>;
type Features = HashMap<String, Vec<Feature>>;

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

pub fn calculate_tpms(counts: &Counts, features: &Features) -> Result<HashMap<String, f64>, Error> {
    let cpbs: HashMap<String, f64> = counts
        .iter()
        .map(|(name, &count)| {
            features
                .get(name)
                .map(|intervals| {
                    let len = sum_nonoverlapping_interval_lengths(intervals);
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
    features: &Features,
) -> Result<HashMap<String, f64>, Error> {
    let counts_sum = sum_counts(counts);

    counts
        .iter()
        .map(|(name, &count)| {
            features
                .get(name)
                .map(|intervals| {
                    let len = sum_nonoverlapping_interval_lengths(intervals);
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

fn sum_nonoverlapping_interval_lengths(features: &[Feature]) -> u64 {
    merge_features(features).iter().map(|f| f.len()).sum()
}

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
