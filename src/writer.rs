use std::{
    collections::HashMap,
    hash::BuildHasher,
    io::{self, Write},
};

pub fn write_normalized_count_values<W, S: BuildHasher>(
    writer: &mut W,
    values: &HashMap<String, f64, S>,
    ids: &[String],
) -> io::Result<()>
where
    W: Write,
{
    for id in ids {
        let value = values.get(id).unwrap_or(&0.0);
        writeln!(writer, "{}\t{}", id, value)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_normalized_count_values() {
        let values: HashMap<String, f64> = [
            (String::from("AADAT"), 30.2),
            (String::from("CLN3"), 3.7),
            (String::from("PAK4"), 14.5),
        ]
        .iter()
        .cloned()
        .collect();

        let ids = vec![
            String::from("AADAT"),
            String::from("CLN3"),
            String::from("NEO1"),
            String::from("PAK4"),
        ];

        let mut buf = Vec::new();

        write_normalized_count_values(&mut buf, &values, &ids).unwrap();

        let actual = String::from_utf8(buf).unwrap();

        let expected = "\
AADAT\t30.2
CLN3\t3.7
NEO1\t0
PAK4\t14.5
";

        assert_eq!(actual, expected);
    }
}
