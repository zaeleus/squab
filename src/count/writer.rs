use std::io::{self, Write};

use super::{context::Counts, Context};

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

    pub fn write_counts(&mut self, ids: &[&String], counts: &Counts) -> io::Result<()> {
        for id in ids {
            let count = counts.get(id.as_str()).unwrap_or(&0);
            writeln!(self.inner, "{id}\t{count}")?;
        }

        Ok(())
    }

    pub fn write_stats(&mut self, ctx: &Context) -> io::Result<()> {
        writeln!(self.inner, "__no_feature\t{}", ctx.miss)?;
        writeln!(self.inner, "__ambiguous\t{}", ctx.ambiguous)?;
        writeln!(self.inner, "__too_low_aQual\t{}", ctx.low_quality)?;
        writeln!(self.inner, "__not_aligned\t{}", ctx.unmapped)?;
        writeln!(self.inner, "__alignment_not_unique\t{}", ctx.nonunique)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_counts() -> io::Result<()> {
        let counts = [("AADAT", 302), ("CLN3", 37), ("PAK4", 145)]
            .into_iter()
            .collect();

        let ids = [
            String::from("AADAT"),
            String::from("CLN3"),
            String::from("NEO1"),
            String::from("PAK4"),
        ];

        let ids_refs: Vec<_> = ids.iter().collect();

        let mut writer = Writer::new(Vec::new());
        writer.write_counts(&ids_refs, &counts)?;

        let actual = writer.get_ref();
        let expected = b"\
AADAT\t302
CLN3\t37
NEO1\t0
PAK4\t145
";

        assert_eq!(&actual[..], &expected[..]);

        Ok(())
    }

    #[test]
    fn test_write_stats() -> io::Result<()> {
        let ctx = Context {
            miss: 735,
            ambiguous: 5,
            low_quality: 60,
            unmapped: 8,
            nonunique: 13,
            ..Default::default()
        };

        let mut writer = Writer::new(Vec::new());
        writer.write_stats(&ctx)?;

        let actual = writer.get_ref();
        let expected = b"\
__no_feature\t735
__ambiguous\t5
__too_low_aQual\t60
__not_aligned\t8
__alignment_not_unique\t13
";

        assert_eq!(&actual[..], &expected[..]);

        Ok(())
    }
}
