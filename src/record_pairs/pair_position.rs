use std::convert::TryFrom;

use noodles_bam as bam;
use noodles_sam as sam;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum PairPosition {
    First,
    Second,
}

impl PairPosition {
    pub fn mate(self) -> PairPosition {
        match self {
            PairPosition::First => PairPosition::Second,
            PairPosition::Second => PairPosition::First,
        }
    }
}

impl<'a> TryFrom<&'a bam::Record> for PairPosition {
    type Error = ();

    fn try_from(record: &bam::Record) -> Result<Self, Self::Error> {
        Self::try_from(record.flags())
    }
}

impl TryFrom<sam::record::Flags> for PairPosition {
    type Error = ();

    fn try_from(flags: sam::record::Flags) -> Result<Self, Self::Error> {
        if flags.is_read_1() {
            Ok(PairPosition::First)
        } else if flags.is_read_2() {
            Ok(PairPosition::Second)
        } else {
            Err(())
        }
    }
}

#[cfg(test)]
mod tests {
    use std::convert::TryFrom;

    use noodles_sam as sam;

    use super::PairPosition;

    #[test]
    fn test_mate() {
        assert_eq!(PairPosition::First.mate(), PairPosition::Second);
        assert_eq!(PairPosition::Second.mate(), PairPosition::First);
    }

    #[test]
    fn test_try_from_flag() {
        let flags = sam::record::Flags::from(0x41);
        assert_eq!(PairPosition::try_from(flags), Ok(PairPosition::First));

        let flags = sam::record::Flags::from(0x81);
        assert_eq!(PairPosition::try_from(flags), Ok(PairPosition::Second));

        let flags = sam::record::Flags::from(0x01);
        assert!(PairPosition::try_from(flags).is_err());
    }
}
