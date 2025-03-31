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
pub enum TryFromFlagsError {
    #[error("ambiguous segment position")]
    Ambiguous,
    #[error("missing segment position")]
    Missing,
}

impl TryFrom<sam::alignment::record::Flags> for SegmentPosition {
    type Error = TryFromFlagsError;

    fn try_from(flags: sam::alignment::record::Flags) -> Result<Self, Self::Error> {
        match (flags.is_first_segment(), flags.is_last_segment()) {
            (true, true) => Err(TryFromFlagsError::Ambiguous),
            (true, false) => Ok(SegmentPosition::First),
            (false, true) => Ok(SegmentPosition::Last),
            (false, false) => Err(TryFromFlagsError::Missing),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mate() {
        assert_eq!(SegmentPosition::First.mate(), SegmentPosition::Last);
        assert_eq!(SegmentPosition::Last.mate(), SegmentPosition::First);
    }

    #[test]
    fn test_try_from_flag() {
        use sam::alignment::record::Flags;

        assert_eq!(
            SegmentPosition::try_from(Flags::SEGMENTED | Flags::FIRST_SEGMENT),
            Ok(SegmentPosition::First)
        );

        assert_eq!(
            SegmentPosition::try_from(Flags::SEGMENTED | Flags::LAST_SEGMENT),
            Ok(SegmentPosition::Last)
        );

        assert_eq!(
            SegmentPosition::try_from(
                Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::LAST_SEGMENT
            ),
            Err(TryFromFlagsError::Ambiguous)
        );

        assert_eq!(
            SegmentPosition::try_from(Flags::SEGMENTED),
            Err(TryFromFlagsError::Missing)
        );
    }
}
