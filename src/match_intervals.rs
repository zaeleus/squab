use std::{io, ops::RangeInclusive};

use noodles_bam::{self as bam, record::cigar};
use noodles_sam as sam;

pub struct MatchIntervals<'a> {
    ops: cigar::Ops<'a>,
    prev_start: u64,
}

impl<'a> MatchIntervals<'a> {
    pub fn new(cigar: &'a bam::record::Cigar, initial_start: u64) -> Self {
        Self {
            ops: cigar.ops(),
            prev_start: initial_start,
        }
    }
}

impl<'a> Iterator for MatchIntervals<'a> {
    type Item = io::Result<RangeInclusive<u64>>;

    fn next(&mut self) -> Option<Self::Item> {
        use sam::record::cigar::op::Kind;

        loop {
            let op = match self.ops.next() {
                Some(Ok(o)) => o,
                Some(Err(e)) => return Some(Err(e)),
                None => return None,
            };

            let len = u64::from(op.len());

            match op.kind() {
                Kind::Match | Kind::SeqMatch | Kind::SeqMismatch => {
                    let start = self.prev_start;
                    let end = start + len - 1;
                    self.prev_start += len;
                    return Some(Ok(start..=end));
                }
                Kind::Deletion | Kind::Skip => {
                    self.prev_start += len;
                }
                _ => continue,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use noodles_bam::{self as bam, record::cigar};
    use noodles_sam::record::cigar::op;

    use super::*;

    fn build_raw_cigar() -> Vec<u8> {
        let ops = [
            u32::from(cigar::Op::new(op::Kind::Match, 1)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::Insertion, 2)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::Deletion, 3)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::Skip, 4)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::SoftClip, 5)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::HardClip, 6)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::Pad, 7)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::SeqMatch, 8)).to_le_bytes(),
            u32::from(cigar::Op::new(op::Kind::SeqMismatch, 9)).to_le_bytes(),
        ];

        ops.iter().flatten().copied().collect()
    }

    #[test]
    fn test_next() -> io::Result<()> {
        let raw_cigar = build_raw_cigar();
        let cigar = bam::record::Cigar::new(&raw_cigar);

        let start = 1;
        let mut it = MatchIntervals::new(&cigar, start);

        assert_eq!(it.next().transpose()?, Some(1..=1));
        assert_eq!(it.next().transpose()?, Some(9..=16));
        assert_eq!(it.next().transpose()?, Some(17..=25));
        assert_eq!(it.next().transpose()?, None);

        Ok(())
    }
}
