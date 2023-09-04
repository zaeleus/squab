use std::{io, ops::RangeInclusive};

use noodles::{bam::lazy::record::Cigar, core::Position, sam::record::cigar::Op};

pub struct MatchIntervals<'a> {
    ops: Box<dyn Iterator<Item = io::Result<Op>> + 'a>,
    prev_alignment_start: Position,
}

impl<'a> MatchIntervals<'a> {
    pub fn new(cigar: &'a Cigar<'a>, alignment_start: Position) -> Self {
        Self {
            ops: Box::new(cigar.iter()),
            prev_alignment_start: alignment_start,
        }
    }
}

impl<'a> Iterator for MatchIntervals<'a> {
    type Item = io::Result<RangeInclusive<Position>>;

    fn next(&mut self) -> Option<Self::Item> {
        use noodles::sam::record::cigar::op::Kind;

        loop {
            let op = match self.ops.next()? {
                Ok(op) => op,
                Err(e) => return Some(Err(e)),
            };

            let len = op.len();

            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    let start = self.prev_alignment_start;

                    let end = start
                        .checked_add(len - 1)
                        .expect("attempt to add with overflow");

                    self.prev_alignment_start = self
                        .prev_alignment_start
                        .checked_add(len)
                        .expect("attempt to add with overflow");

                    return Some(Ok(start..=end));
                }
                Kind::Deletion | Kind::Skip => {
                    self.prev_alignment_start = self
                        .prev_alignment_start
                        .checked_add(len)
                        .expect("attempt to add with overflow");
                }
                _ => continue,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use noodles::{bam, sam};

    use super::*;

    #[test]
    fn test_next() -> Result<(), Box<dyn std::error::Error>> {
        let record = sam::alignment::Record::builder()
            .set_alignment_start(Position::MIN)
            .set_cigar("1M2I3D4N5S6H7P8=9X".parse()?)
            .build();

        let mut writer = bam::Writer::from(Vec::new());
        writer.write_record(&sam::Header::default(), &record)?;
        let buf = writer.into_inner()[4..].to_vec();

        let raw_record = bam::lazy::Record::try_from(buf)?;

        let cigar = raw_record.cigar();
        let start = raw_record
            .alignment_start()?
            .expect("missing alignment start");
        let mut it = MatchIntervals::new(&cigar, start);

        assert_eq!(
            it.next().transpose()?,
            Some(Position::try_from(1)?..=Position::try_from(1)?)
        );
        assert_eq!(
            it.next().transpose()?,
            Some(Position::try_from(9)?..=Position::try_from(16)?)
        );
        assert_eq!(
            it.next().transpose()?,
            Some(Position::try_from(17)?..=Position::try_from(25)?)
        );
        assert!(it.next().is_none());

        Ok(())
    }
}
