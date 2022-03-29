use std::collections::HashMap;

use super::{sum_nonoverlapping_feature_lengths, Error, FeatureMap};
use crate::count::context::Counts;

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

fn sum_counts(counts: &Counts) -> u64 {
    counts.values().sum()
}

fn calculate_fpkm(count: u64, len: usize, counts_sum: u64) -> f64 {
    (count as f64 * 1e9) / (len as f64 * counts_sum as f64)
}

#[cfg(test)]
mod tests {
    use crate::Feature;

    use noodles::{core::Position, gff};

    use super::*;

    fn build_counts() -> HashMap<String, u64> {
        let counts = [
            (String::from("AAAS"), 645),
            (String::from("AC009952.3"), 1),
            (String::from("RPL37AP1"), 5714),
        ];

        counts.iter().cloned().collect()
    }

    fn build_feature_map() -> Result<FeatureMap, noodles::core::position::TryFromIntError> {
        let reference_name = String::from("chr1");
        let strand = gff::record::Strand::Forward;

        let features = [
            (
                String::from("AAAS"),
                vec![Feature::new(
                    reference_name.clone(),
                    Position::try_from(53307456)?,
                    Position::try_from(53324864)?,
                    strand,
                )],
            ),
            (
                String::from("AC009952.3"),
                vec![Feature::new(
                    reference_name.clone(),
                    Position::try_from(9189629)?,
                    Position::try_from(9204611)?,
                    strand,
                )],
            ),
            (
                String::from("RPL37AP1"),
                vec![Feature::new(
                    reference_name,
                    Position::try_from(44466564)?,
                    Position::try_from(44466842)?,
                    strand,
                )],
            ),
        ];

        Ok(features.iter().cloned().collect())
    }

    #[test]
    fn test_calculate_fpkms() -> Result<(), Box<dyn std::error::Error>> {
        let counts = build_counts();

        let feature_map = build_feature_map()?;
        let fpkms = calculate_fpkms(&counts, &feature_map)?;

        assert_eq!(fpkms.len(), 3);

        let a = fpkms["AAAS"];
        let b = 5825.440538780093;
        assert!((a - b).abs() < f64::EPSILON);

        let a = fpkms["AC009952.3"];
        let b = 10.494073576888189;
        assert!((a - b).abs() < f64::EPSILON);

        let a = fpkms["RPL37AP1"];
        let b = 3220170.8708099457;
        assert!((a - b).abs() < f64::EPSILON);

        Ok(())
    }

    #[test]
    fn test_calculate_fpkms_with_missing_feature(
    ) -> Result<(), noodles::core::position::TryFromIntError> {
        let counts = build_counts();

        let mut feature_map = build_feature_map()?;
        feature_map.remove("AC009952.3");

        assert!(calculate_fpkms(&counts, &feature_map).is_err());

        Ok(())
    }
}
