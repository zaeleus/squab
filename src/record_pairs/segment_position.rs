use noodles::sam;
use thiserror::Error;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum SegmentPosition {
    First,
    Last,
}

impl SegmentPosition {
    pub fn mate(self) -> SegmentPosition {
        match self {
            Self::First => SegmentPosition::Last,
            Self::Last => SegmentPosition::First,
        }
    }
}

#[derive(Clone, Debug, Error, Eq, PartialEq)]
#[error("neither read 1 nor read 2 flag is set")]
pub struct TryFromFlagsError;

impl TryFrom<sam::alignment::record::Flags> for SegmentPosition {
    type Error = TryFromFlagsError;

    fn try_from(flags: sam::alignment::record::Flags) -> Result<Self, Self::Error> {
        if flags.is_first_segment() {
            Ok(SegmentPosition::First)
        } else if flags.is_last_segment() {
            Ok(SegmentPosition::Last)
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
        assert_eq!(SegmentPosition::First.mate(), SegmentPosition::Last);
        assert_eq!(SegmentPosition::Last.mate(), SegmentPosition::First);
    }

    #[test]
    fn test_try_from_flag() {
        use sam::alignment::record::Flags;

        let flags = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        assert_eq!(SegmentPosition::try_from(flags), Ok(SegmentPosition::First));

        let flags = Flags::SEGMENTED | Flags::LAST_SEGMENT;
        assert_eq!(SegmentPosition::try_from(flags), Ok(SegmentPosition::Last));

        let flags = Flags::SEGMENTED;
        assert_eq!(SegmentPosition::try_from(flags), Err(TryFromFlagsError));
    }
}
