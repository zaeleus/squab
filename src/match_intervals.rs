use std::{io, num::NonZero, ops::RangeInclusive};

use noodles::{
    core::Position,
    sam::alignment::record::cigar::{Op, op::Kind},
};

pub struct MatchIntervals<'r> {
    ops: &'r mut dyn Iterator<Item = io::Result<Op>>,
    reference_position: Position,
}

impl<'r> MatchIntervals<'r> {
    pub fn new(
        ops: &'r mut dyn Iterator<Item = io::Result<Op>>,
        alignment_start: Position,
    ) -> Self {
        Self {
            ops,
            reference_position: alignment_start,
        }
    }

    fn consume_reference(&mut self, len: NonZero<usize>) {
        self.reference_position = self
            .reference_position
            .checked_add(len.get())
            .expect("attempt to add with overflow");
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

            let Some(len) = NonZero::new(op.len()) else {
                continue;
            };

            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    let interval = build_interval(self.reference_position, len);
                    self.consume_reference(len);
                    return Some(Ok(interval));
                }
                Kind::Deletion | Kind::Skip => self.consume_reference(len),
                _ => continue,
            }
        }
    }
}

fn build_interval(start: Position, len: NonZero<usize>) -> RangeInclusive<Position> {
    let end = start
        .checked_add(len.get() - 1)
        .expect("attempt to add with overflow");

    start..=end
}

#[cfg(test)]
mod tests {
    use std::num;

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

    #[test]
    fn test_build_interval() -> Result<(), num::TryFromIntError> {
        assert_eq!(
            build_interval(Position::MIN, NonZero::<usize>::MIN),
            Position::MIN..=Position::MIN
        );

        assert_eq!(
            build_interval(Position::MIN, NonZero::try_from(4)?),
            Position::MIN..=Position::try_from(4)?
        );

        assert_eq!(
            build_interval(Position::try_from(3)?, NonZero::try_from(5)?),
            Position::try_from(3)?..=Position::try_from(7)?
        );

        Ok(())
    }
}
