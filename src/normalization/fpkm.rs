pub fn normalize(lengths: &[u32], counts: &[u32]) -> Vec<f64> {
    let sum = counts.iter().copied().map(f64::from).sum();

    lengths
        .iter()
        .zip(counts)
        .map(|(length, count)| calculate_fpkm(*count, *length, sum))
        .collect()
}

fn calculate_fpkm(count: u32, length: u32, sum: f64) -> f64 {
    assert!(length > 0);
    (f64::from(count) * 1e9) / (f64::from(length) * sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize() {
        let counts = [645, 1, 5714];
        let lengths = [17409, 14983, 279];

        let fpkms = normalize(&lengths, &counts);

        assert_eq!(fpkms.len(), 3);

        let a = fpkms[0];
        let b = 5825.440538780093;
        assert!((a - b).abs() < f64::EPSILON);

        let a = fpkms[1];
        let b = 10.494073576888189;
        assert!((a - b).abs() < f64::EPSILON);

        let a = fpkms[2];
        let b = 3220170.8708099457;
        assert!((a - b).abs() < f64::EPSILON);
    }
}
