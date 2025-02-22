use std::{io, ops::RangeInclusive};

use noodles::{
    core::Position,
    sam::alignment::record::cigar::{Op, op::Kind},
};

pub struct MatchIntervals<'r> {
    ops: &'r mut dyn Iterator<Item = io::Result<Op>>,
    prev_alignment_start: Position,
}

impl<'r> MatchIntervals<'r> {
    pub fn new(
        ops: &'r mut dyn Iterator<Item = io::Result<Op>>,
        initial_alignment_start: Position,
    ) -> Self {
        Self {
            ops,
            prev_alignment_start: initial_alignment_start,
        }
    }
}

impl Iterator for MatchIntervals<'_> {
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
    use super::*;

    #[test]
    fn test_next() -> Result<(), Box<dyn std::error::Error>> {
        let mut ops = [
            Ok(Op::new(Kind::Match, 1)),
            Ok(Op::new(Kind::Insertion, 2)),
            Ok(Op::new(Kind::Deletion, 3)),
            Ok(Op::new(Kind::Skip, 4)),
            Ok(Op::new(Kind::SoftClip, 5)),
            Ok(Op::new(Kind::HardClip, 6)),
            Ok(Op::new(Kind::Pad, 7)),
            Ok(Op::new(Kind::SequenceMatch, 8)),
            Ok(Op::new(Kind::SequenceMismatch, 9)),
        ]
        .into_iter();

        let intervals = MatchIntervals::new(&mut ops, Position::MIN);
        let actual: Vec<_> = intervals.collect::<io::Result<_>>()?;

        let expected = [
            Position::try_from(1)?..=Position::try_from(1)?,
            Position::try_from(9)?..=Position::try_from(16)?,
            Position::try_from(17)?..=Position::try_from(25)?,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
