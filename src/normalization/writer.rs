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

    pub fn write_values(&mut self, names: &[String], values: &HashMap<&str, f64>) -> io::Result<()>
    where
        W: Write,
    {
        for name in names {
            let value = values.get(name.as_str()).unwrap_or(&0.0);
            writeln!(self.inner, "{name}\t{value}")?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_values() -> io::Result<()> {
        let names = [
            String::from("AADAT"),
            String::from("CLN3"),
            String::from("NEO1"),
            String::from("PAK4"),
        ];

        let values = [("AADAT", 30.2), ("CLN3", 3.7), ("PAK4", 14.5)]
            .into_iter()
            .collect();

        let mut writer = Writer::new(Vec::new());
        writer.write_values(&names, &values)?;

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
