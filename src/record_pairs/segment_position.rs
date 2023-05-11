use noodles::sam;
use thiserror::Error;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum SegmentPosition {
    First,
    Second,
}

impl SegmentPosition {
    pub fn mate(self) -> SegmentPosition {
        match self {
            Self::First => SegmentPosition::Second,
            Self::Second => SegmentPosition::First,
        }
    }
}

#[derive(Clone, Debug, Error, Eq, PartialEq)]
#[error("neither read 1 nor read 2 flag is set")]
pub struct TryFromFlagsError;

impl TryFrom<sam::record::Flags> for SegmentPosition {
    type Error = TryFromFlagsError;

    fn try_from(flags: sam::record::Flags) -> Result<Self, Self::Error> {
        if flags.is_first_segment() {
            Ok(SegmentPosition::First)
        } else if flags.is_last_segment() {
            Ok(SegmentPosition::Second)
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
        assert_eq!(SegmentPosition::First.mate(), SegmentPosition::Second);
        assert_eq!(SegmentPosition::Second.mate(), SegmentPosition::First);
    }

    #[test]
    fn test_try_from_flag() {
        use sam::record::Flags;

        let flags = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        assert_eq!(SegmentPosition::try_from(flags), Ok(SegmentPosition::First));

        let flags = Flags::SEGMENTED | Flags::LAST_SEGMENT;
        assert_eq!(
            SegmentPosition::try_from(flags),
            Ok(SegmentPosition::Second)
        );

        let flags = Flags::SEGMENTED;
        assert_eq!(SegmentPosition::try_from(flags), Err(TryFromFlagsError));
    }
}
