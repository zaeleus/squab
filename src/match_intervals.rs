use std::{io, ops::RangeInclusive};

use noodles::{
    bam::record::Cigar,
    core::Position,
    sam::alignment::record::cigar::{op::Kind, Op},
};

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
    use noodles::{
        bam,
        sam::{self, alignment::io::Write},
    };

    use super::*;

    #[test]
    fn test_next() -> Result<(), Box<dyn std::error::Error>> {
        let record = sam::alignment::RecordBuf::builder()
            .set_alignment_start(Position::MIN)
            .set_cigar(
                [
                    Op::new(Kind::Match, 1),
                    Op::new(Kind::Insertion, 2),
                    Op::new(Kind::Deletion, 3),
                    Op::new(Kind::Skip, 4),
                    Op::new(Kind::SoftClip, 5),
                    Op::new(Kind::HardClip, 6),
                    Op::new(Kind::Pad, 7),
                    Op::new(Kind::SequenceMatch, 8),
                    Op::new(Kind::SequenceMismatch, 9),
                ]
                .into_iter()
                .collect(),
            )
            .build();

        let mut writer = bam::io::Writer::from(Vec::new());
        writer.write_alignment_record(&sam::Header::default(), &record)?;

        let mut reader = bam::io::Reader::from(writer.get_ref().as_slice());

        let mut record = bam::Record::default();
        reader.read_record(&mut record)?;

        let cigar = record.cigar();
        let start = record.alignment_start().expect("missing alignment start")?;
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
