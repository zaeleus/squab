use std::{
    io::{self, BufRead},
    num,
};

use bstr::BString;
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

pub fn read<R>(reader: &mut R) -> Result<Vec<(BString, u32)>, ReadCountsError>
where
    R: BufRead,
{
    let mut line = String::new();
    let mut counts = Vec::new();

    loop {
        line.clear();

        if read_line(reader, &mut line)? == 0 {
            break;
        }

        if line.starts_with(HTSEQ_COUNT_META_PREFIX) {
            break;
        }

        let (name, count) = parse_line(&line)?;
        counts.push((name.into(), count));
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

fn parse_line(s: &str) -> Result<(&str, u32), ReadCountsError> {
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
        let actual = read(&mut reader)?;

        let expected = [
            (BString::from("AADAT"), 302),
            (BString::from("CLN3"), 37),
            (BString::from("PAK4"), 145),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
