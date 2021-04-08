use std::collections::HashMap;

use super::{sum_nonoverlapping_feature_lengths, Counts, Error, FeatureMap};

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

fn calculate_tpm(cpb: f64, cpbs_sum: f64) -> f64 {
    cpb * 1e6 / cpbs_sum
}

#[cfg(test)]
mod tests {
    use crate::Feature;

    use noodles::gff;

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
        let strand = gff::record::Strand::Forward;

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
    fn test_calculate_tpms() -> Result<(), Error> {
        let counts = build_counts();
        let feature_map = build_feature_map();

        let tpms = calculate_tpms(&counts, &feature_map)?;

        assert_eq!(tpms.len(), 3);

        let a = tpms["AAAS"];
        let b = 1805.7744109493626;
        assert!((a - b).abs() < f64::EPSILON);

        let a = tpms["AC009952.3"];
        let b = 3.252960768479983;
        assert!((a - b).abs() < f64::EPSILON);

        let a = tpms["RPL37AP1"];
        let b = 998190.972628282;
        assert!((a - b).abs() < f64::EPSILON);

        Ok(())
    }

    #[test]
    fn test_calculate_tpms_with_missing_feature() {
        let counts = build_counts();

        let mut feature_map = build_feature_map();
        feature_map.remove("AC009952.3");

        assert!(calculate_tpms(&counts, &feature_map).is_err());
    }
}
