use std::{
    collections::HashMap,
    io::{self, Write},
};

use crate::count::Context;

pub fn write_counts<W>(
    writer: &mut W,
    counts: &HashMap<String, u64>,
    ids: &[String],
) -> io::Result<()>
where
    W: Write,
{
    for id in ids {
        let count = counts.get(id).unwrap_or(&0);
        writeln!(writer, "{}\t{}", id, count)?;
    }

    Ok(())
}

pub fn write_stats<W>(writer: &mut W, ctx: &Context) -> io::Result<()>
where
    W: Write,
{
    writeln!(writer, "__no_feature\t{}", ctx.no_feature)?;
    writeln!(writer, "__ambiguous\t{}", ctx.ambiguous)?;
    writeln!(writer, "__too_low_aQual\t{}", ctx.low_quality)?;
    writeln!(writer, "__not_aligned\t{}", ctx.unmapped)?;
    writeln!(writer, "__alignment_not_unique\t{}", ctx.nonunique)?;
    Ok(())
}

pub fn write_normalized_count_values<W>(
    writer: &mut W,
    values: &HashMap<String, f64>,
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
    fn test_write_counts() {
        let counts: HashMap<String, u64> = [
            (String::from("AADAT"), 302),
            (String::from("CLN3"), 37),
            (String::from("PAK4"), 145),
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

        write_counts(&mut buf, &counts, &ids).unwrap();

        let actual = String::from_utf8(buf).unwrap();

        let expected = "\
AADAT\t302
CLN3\t37
NEO1\t0
PAK4\t145
";

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_write_stats() {
        let mut buf = Vec::new();

        let mut ctx = Context::default();
        ctx.no_feature = 735;
        ctx.ambiguous = 5;
        ctx.low_quality = 60;
        ctx.unmapped = 8;
        ctx.nonunique = 13;

        write_stats(&mut buf, &ctx).unwrap();

        let actual = String::from_utf8(buf).unwrap();

        let expected = "\
__no_feature\t735
__ambiguous\t5
__too_low_aQual\t60
__not_aligned\t8
__alignment_not_unique\t13
";

        assert_eq!(actual, expected);
    }

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
