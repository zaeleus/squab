use std::{
    collections::HashMap,
    io::{self, BufRead},
    num,
};

use thiserror::Error;

static HTSEQ_COUNT_META_PREFIX: &str = "__";

#[derive(Debug, Error)]
pub enum ReadCountsError {
    #[error("invalid record")]
    InvalidRecord,
    #[error("invalid count")]
    InvalidCount(num::ParseIntError),
    #[error("I/O error")]
    Io(#[from] io::Error),
}

pub fn read<R>(reader: &mut R) -> Result<HashMap<String, u64>, ReadCountsError>
where
    R: BufRead,
{
    let mut line = String::new();
    let mut counts = HashMap::new();

    loop {
        line.clear();

        if read_line(reader, &mut line)? == 0 {
            break;
        }

        let (name, count) = parse_line(&line)?;

        if name.starts_with(HTSEQ_COUNT_META_PREFIX) {
            break;
        }

        counts.insert(name.into(), count);
    }

    Ok(counts)
}

fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: char = '\n';
    const CARRIAGE_RETURN: char = '\r';

    match reader.read_line(buf)? {
        0 => Ok(0),
        n => {
            if buf.ends_with(LINE_FEED) {
                buf.pop();

                if buf.ends_with(CARRIAGE_RETURN) {
                    buf.pop();
                }
            }

            Ok(n)
        }
    }
}

fn parse_line(s: &str) -> Result<(&str, u64), ReadCountsError> {
    const DELIMITER: char = '\t';

    let (name, raw_count) = s
        .split_once(DELIMITER)
        .ok_or(ReadCountsError::InvalidRecord)?;

    let count = raw_count.parse().map_err(ReadCountsError::InvalidCount)?;

    Ok((name, count))
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
    fn test_read_counts() -> Result<(), ReadCountsError> {
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

        let mut reader = &data[..];
        let counts = read(&mut reader)?;

        assert_eq!(counts.len(), 3);
        assert_eq!(counts["AADAT"], 302);
        assert_eq!(counts["CLN3"], 37);
        assert_eq!(counts["PAK4"], 145);

        Ok(())
    }
}
