use std::collections::HashMap;

use super::{Counts, Error};

pub fn calculate_tpms<'f>(
    lengths: &'f HashMap<String, usize>,
    counts: &Counts,
) -> Result<HashMap<&'f str, f64>, Error> {
    let cpbs: HashMap<_, _> = counts
        .iter()
        .map(|(key, &count)| {
            lengths
                .get_key_value(key)
                .map(|(name, len)| {
                    let cpb = (count as f64) / (*len as f64);
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

    fn build_lengths() -> HashMap<String, usize> {
        [
            (String::from("AAAS"), 17409),
            (String::from("AC009952.3"), 14983),
            (String::from("RPL37AP1"), 279),
        ]
        .into_iter()
        .collect()
    }

    #[test]
    fn test_calculate_tpms() -> Result<(), Box<dyn std::error::Error>> {
        let lengths = build_lengths();
        let counts = build_counts();
        let tpms = calculate_tpms(&lengths, &counts)?;

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
        let lengths = HashMap::new();
        let counts = build_counts();
        assert!(calculate_tpms(&lengths, &counts).is_err());
        Ok(())
    }
}
