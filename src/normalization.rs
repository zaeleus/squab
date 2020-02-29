use std::collections::HashMap;

use crate::Feature;

type Counts = HashMap<String, u64>;
type Features = HashMap<String, Vec<Feature>>;

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
