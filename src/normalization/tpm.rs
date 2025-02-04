pub fn normalize(lengths: &[usize], counts: &[u32]) -> Vec<f64> {
    let length_normalized_counts: Vec<_> = lengths
        .iter()
        .zip(counts)
        .map(|(length, count)| {
            assert!(*length > 0);
            f64::from(*count) / (*length as f64)
        })
        .collect();

    let sum = length_normalized_counts.iter().sum();

    length_normalized_counts
        .into_iter()
        .map(|n| calculate_tpm(n, sum))
        .collect()
}

fn calculate_tpm(n: f64, sum: f64) -> f64 {
    n * 1e6 / sum
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize() {
        let counts = [645, 1, 5714];
        let lengths = [17409, 14983, 279];

        let tpms = normalize(&lengths, &counts);

        assert_eq!(tpms.len(), 3);

        let a = tpms[0];
        let b = 1805.7744109493626;
        assert!((a - b).abs() < f64::EPSILON);

        let a = tpms[1];
        let b = 3.252960768479983;
        assert!((a - b).abs() < f64::EPSILON);

        let a = tpms[2];
        let b = 998190.972628282;
        assert!((a - b).abs() < f64::EPSILON);
    }
}
