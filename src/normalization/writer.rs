use std::{
    collections::HashMap,
    io::{self, Write},
};

pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W>
where
    W: Write,
{
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    pub fn write_values(&mut self, ids: &[String], values: &HashMap<String, f64>) -> io::Result<()>
    where
        W: Write,
    {
        for id in ids {
            let value = values.get(id).unwrap_or(&0.0);
            writeln!(self.inner, "{}\t{}", id, value)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_values() -> io::Result<()> {
        let values: HashMap<String, f64> = vec![
            (String::from("AADAT"), 30.2),
            (String::from("CLN3"), 3.7),
            (String::from("PAK4"), 14.5),
        ]
        .into_iter()
        .collect();

        let ids = [
            String::from("AADAT"),
            String::from("CLN3"),
            String::from("NEO1"),
            String::from("PAK4"),
        ];

        let mut writer = Writer::new(Vec::new());
        writer.write_values(&ids, &values)?;

        let actual = writer.get_ref();
        let expected = b"\
AADAT\t30.2
CLN3\t3.7
NEO1\t0
PAK4\t14.5
";

        assert_eq!(&actual[..], &expected[..]);

        Ok(())
    }
}
