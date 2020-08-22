use std::ops::RangeInclusive;

use noodles_bam::{self as bam, record::cigar};
use noodles_sam as sam;

pub struct CigarToIntervals<'a> {
    ops: cigar::Ops<'a>,
    prev_start: u64,
}

impl<'a> CigarToIntervals<'a> {
    pub fn new(cigar: &'a bam::record::Cigar, initial_start: u64) -> CigarToIntervals<'a> {
        let ops = cigar.ops();

        CigarToIntervals {
            ops,
            prev_start: initial_start,
        }
    }
}

impl<'a> Iterator for CigarToIntervals<'a> {
    type Item = RangeInclusive<u64>;

    fn next(&mut self) -> Option<Self::Item> {
        use sam::record::cigar::op::Kind;

        loop {
            let op = self.ops.next()?;
            let len = u64::from(op.len());

            match op.kind() {
                Kind::Match | Kind::SeqMatch | Kind::SeqMismatch => {
                    let start = self.prev_start;
                    let end = start + len - 1;
                    self.prev_start += len;
                    return Some(start..=end);
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

    use super::CigarToIntervals;

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
    fn test_next() {
        let raw_cigar = build_raw_cigar();
        let cigar = bam::record::Cigar::new(&raw_cigar);

        let start = 1;
        let mut it = CigarToIntervals::new(&cigar, start);

        assert_eq!(it.next(), Some(1..=1));
        assert_eq!(it.next(), Some(9..=16));
        assert_eq!(it.next(), Some(17..=25));
        assert!(it.next().is_none());
    }
}
