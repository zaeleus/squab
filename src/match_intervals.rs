use std::{ops::RangeInclusive, slice};

use noodles::{
    core::Position,
    sam::record::{cigar::Op, Cigar},
};

pub struct MatchIntervals<'a> {
    ops: slice::Iter<'a, Op>,
    prev_alignment_start: Position,
}

impl<'a> MatchIntervals<'a> {
    pub fn new(cigar: &'a Cigar, alignment_start: Position) -> Self {
        Self {
            ops: cigar.iter(),
            prev_alignment_start: alignment_start,
        }
    }
}

impl<'a> Iterator for MatchIntervals<'a> {
    type Item = RangeInclusive<Position>;

    fn next(&mut self) -> Option<Self::Item> {
        use noodles::sam::record::cigar::op::Kind;

        loop {
            let op = self.ops.next()?;
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

                    return Some(start..=end);
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
        let cigar = "1M2I3D4N5S6H7P8=9X".parse()?;

        let start = Position::MIN;
        let mut it = MatchIntervals::new(&cigar, start);

        assert_eq!(
            it.next(),
            Some(Position::try_from(1)?..=Position::try_from(1)?)
        );
        assert_eq!(
            it.next(),
            Some(Position::try_from(9)?..=Position::try_from(16)?)
        );
        assert_eq!(
            it.next(),
            Some(Position::try_from(17)?..=Position::try_from(25)?)
        );
        assert!(it.next().is_none());

        Ok(())
    }
}
