use std::{io, ops::RangeInclusive};

use noodles::{
    bam::{self, record::cigar},
    sam,
};

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
    use noodles::{bam::record::cigar, sam::record::cigar::op};

    use super::*;

    fn build_raw_cigar() -> io::Result<Vec<u8>> {
        fn build_raw_op(kind: op::Kind, len: u32) -> io::Result<[u8; 4]> {
            cigar::Op::new(kind, len)
                .map(u32::from)
                .map(|n| n.to_le_bytes())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        }

        let raw_ops = [
            build_raw_op(op::Kind::Match, 1)?,
            build_raw_op(op::Kind::Insertion, 2)?,
            build_raw_op(op::Kind::Deletion, 3)?,
            build_raw_op(op::Kind::Skip, 4)?,
            build_raw_op(op::Kind::SoftClip, 5)?,
            build_raw_op(op::Kind::HardClip, 6)?,
            build_raw_op(op::Kind::Pad, 7)?,
            build_raw_op(op::Kind::SeqMatch, 8)?,
            build_raw_op(op::Kind::SeqMismatch, 9)?,
        ];

        Ok(raw_ops.iter().flatten().copied().collect())
    }

    #[test]
    fn test_next() -> io::Result<()> {
        let raw_cigar = build_raw_cigar()?;
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
