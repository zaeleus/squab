use std::{error, fmt};

use noodles::{
    bam,
    sam::{self, AlignmentRecord},
};

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
    type Error = TryFromFlagsError;

    fn try_from(record: &bam::Record) -> Result<Self, Self::Error> {
        Self::try_from(record.flags())
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromFlagsError;

impl error::Error for TryFromFlagsError {}

impl fmt::Display for TryFromFlagsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("neither read 1 nor read 2 flag is set")
    }
}

impl TryFrom<sam::record::Flags> for PairPosition {
    type Error = TryFromFlagsError;

    fn try_from(flags: sam::record::Flags) -> Result<Self, Self::Error> {
        if flags.is_first_segment() {
            Ok(PairPosition::First)
        } else if flags.is_last_segment() {
            Ok(PairPosition::Second)
        } else {
            Err(TryFromFlagsError)
        }
    }
}

#[cfg(test)]
mod tests {
    use std::convert::TryFrom;

    use super::*;

    #[test]
    fn test_mate() {
        assert_eq!(PairPosition::First.mate(), PairPosition::Second);
        assert_eq!(PairPosition::Second.mate(), PairPosition::First);
    }

    #[test]
    fn test_try_from_flag() {
        use sam::record::Flags;

        let flags = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        assert_eq!(PairPosition::try_from(flags), Ok(PairPosition::First));

        let flags = Flags::SEGMENTED | Flags::LAST_SEGMENT;
        assert_eq!(PairPosition::try_from(flags), Ok(PairPosition::Second));

        let flags = Flags::SEGMENTED;
        assert_eq!(PairPosition::try_from(flags), Err(TryFromFlagsError));
    }
}
