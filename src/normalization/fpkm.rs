use std::collections::HashMap;

use super::{Counts, Error};

pub fn calculate_fpkms<'f>(
    lengths: &'f HashMap<String, usize>,
    counts: &Counts,
) -> Result<HashMap<&'f str, f64>, Error> {
    let counts_sum = sum_counts(counts);

    counts
        .iter()
        .map(|(key, &count)| {
            lengths
                .get_key_value(key)
                .map(|(name, len)| {
                    let fpkm = calculate_fpkm(count, *len, counts_sum);
                    (name.as_str(), fpkm)
                })
                .ok_or_else(|| Error::MissingFeature(key.into()))
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
    fn test_calculate_fpkms() -> Result<(), Error> {
        let lengths = build_lengths();
        let counts = build_counts();
        let fpkms = calculate_fpkms(&lengths, &counts)?;

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
    fn test_calculate_fpkms_with_missing_feature() {
        let lengths = HashMap::new();
        let counts = build_counts();
        assert!(calculate_fpkms(&lengths, &counts).is_err());
    }
}
