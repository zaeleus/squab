use std::collections::HashMap;

use super::{sum_nonoverlapping_feature_lengths, Counts, Error};
use crate::Features;

pub fn calculate_tpms<'f>(
    counts: &Counts,
    features: &'f Features,
) -> Result<HashMap<&'f str, f64>, Error> {
    let cpbs: HashMap<_, _> = counts
        .iter()
        .map(|(key, &count)| {
            features
                .get_key_value(key)
                .map(|(name, features)| {
                    let len = sum_nonoverlapping_feature_lengths(features);
                    let cpb = count as f64 / len as f64;
                    (name.as_str(), cpb)
                })
                .ok_or_else(|| Error::MissingFeature(key.into()))
        })
        .collect::<Result<_, _>>()?;

    let cpbs_sum = cpbs.values().sum();

    let tpms = cpbs
        .into_iter()
        .map(|(name, cpb)| (name, calculate_tpm(cpb, cpbs_sum)))
        .collect();

    Ok(tpms)
}

fn calculate_tpm(cpb: f64, cpbs_sum: f64) -> f64 {
    cpb * 1e6 / cpbs_sum
}

#[cfg(test)]
mod tests {
    use crate::Feature;

    use noodles::{core::Position, gff};

    use super::*;

    fn build_counts() -> Counts {
        [
            (String::from("AAAS"), 645),
            (String::from("AC009952.3"), 1),
            (String::from("RPL37AP1"), 5714),
        ]
        .into_iter()
        .collect()
    }

    fn build_features() -> Result<Features, noodles::core::position::TryFromIntError> {
        let strand = gff::record::Strand::Forward;

        let features = [
            (
                String::from("AAAS"),
                vec![Feature::new(
                    0,
                    Position::try_from(53307456)?,
                    Position::try_from(53324864)?,
                    strand,
                )],
            ),
            (
                String::from("AC009952.3"),
                vec![Feature::new(
                    0,
                    Position::try_from(9189629)?,
                    Position::try_from(9204611)?,
                    strand,
                )],
            ),
            (
                String::from("RPL37AP1"),
                vec![Feature::new(
                    0,
                    Position::try_from(44466564)?,
                    Position::try_from(44466842)?,
                    strand,
                )],
            ),
        ]
        .into_iter()
        .collect();

        Ok(features)
    }
    #[test]
    fn test_calculate_tpms() -> Result<(), Box<dyn std::error::Error>> {
        let counts = build_counts();

        let features = build_features()?;
        let tpms = calculate_tpms(&counts, &features)?;

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
    fn test_calculate_tpms_with_missing_feature(
    ) -> Result<(), noodles::core::position::TryFromIntError> {
        let counts = build_counts();

        let mut features = build_features()?;
        features.remove("AC009952.3");

        assert!(calculate_tpms(&counts, &features).is_err());

        Ok(())
    }
}
