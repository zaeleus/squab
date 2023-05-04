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

            if read_line(&mut self.inner, &mut buf)? == 0 {
                break;
            }

            let (id, count) = parse_line(&buf)?;

            if id.starts_with(HTSEQ_COUNT_META_PREFIX) {
                break;
            }

            counts.insert(id, count);
        }

        Ok(counts)
    }
}

fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: char = '\n';
    const CARRIAGE_RETURN: char = '\r';

    reader.read_line(buf).map(|n| {
        if n > 0 && buf.ends_with(LINE_FEED) {
            buf.pop();

            if buf.ends_with(CARRIAGE_RETURN) {
                buf.pop();
            }
        }

        n
    })
}

fn parse_line(s: &str) -> io::Result<(String, u64)> {
    let (raw_id, raw_count) = s
        .split_once(DELIMITER)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "malformed record"))?;

    let count = raw_count
        .parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((raw_id.into(), count))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_line() -> io::Result<()> {
        fn t(buf: &mut String, mut reader: &[u8], expected: &str) -> io::Result<()> {
            buf.clear();
            read_line(&mut reader, buf)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = String::new();

        t(&mut buf, b"squab\n", "squab")?;
        t(&mut buf, b"squab\r\n", "squab")?;
        t(&mut buf, b"squab", "squab")?;

        Ok(())
    }

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
