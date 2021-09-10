use std::{
    collections::HashMap,
    io::{self, BufRead},
};

use super::context::Counts;

const DELIMITER: char = '\t';
static HTSEQ_COUNT_META_PREFIX: &str = "__";

pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    pub fn read_counts(&mut self) -> io::Result<Counts> {
        let mut counts = HashMap::default();
        let mut buf = String::new();

        loop {
            buf.clear();

            match self.inner.read_line(&mut buf) {
                Ok(0) => break,
                Ok(_) => {}
                Err(e) => return Err(e),
            }

            buf.pop();

            let mut fields = buf.split(DELIMITER);

            let id = parse_string(&mut fields)?;

            if id.starts_with(HTSEQ_COUNT_META_PREFIX) {
                break;
            }

            let count = parse_u64(&mut fields)?;

            counts.insert(id.into(), count);
        }

        Ok(counts)
    }
}

fn parse_string<'a, I>(fields: &mut I) -> io::Result<&'a str>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData))
}

fn parse_u64<'a, I>(fields: &mut I) -> io::Result<u64>
where
    I: Iterator<Item = &'a str>,
{
    parse_string(fields).and_then(|s| {
        s.parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_counts() -> io::Result<()> {
        let data = b"\
AADAT\t302
CLN3\t37
PAK4\t145
__no_feature\t0
__ambiguous\t0
__too_low_aQual\t0
__not_aligned\t0
__alignment_not_unique\t0
";

        let mut reader = Reader::new(&data[..]);
        let counts = reader.read_counts()?;

        assert_eq!(counts.len(), 3);
        assert_eq!(counts["AADAT"], 302);
        assert_eq!(counts["CLN3"], 37);
        assert_eq!(counts["PAK4"], 145);

        Ok(())
    }
}
