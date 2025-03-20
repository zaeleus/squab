use noodles::{core::Position, gff};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Feature {
    reference_sequence_id: usize,
    start: Position,
    end: Position,
    strand: gff::feature::record::Strand,
}

impl Feature {
    pub fn new(
        reference_sequence_id: usize,
        start: Position,
        end: Position,
        strand: gff::feature::record::Strand,
    ) -> Self {
        Self {
            reference_sequence_id,
            start,
            end,
            strand,
        }
    }

    pub fn reference_sequence_id(&self) -> usize {
        self.reference_sequence_id
    }

    pub fn start(&self) -> Position {
        self.start
    }

    pub fn end(&self) -> Position {
        self.end
    }

    pub fn end_mut(&mut self) -> &mut Position {
        &mut self.end
    }

    pub fn strand(&self) -> gff::feature::record::Strand {
        self.strand
    }

    pub fn length(&self) -> usize {
        usize::from(self.end) - usize::from(self.start) + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_feature() -> Result<Feature, noodles::core::position::TryFromIntError> {
        Ok(Feature::new(
            0,
            Position::try_from(8)?,
            Position::try_from(13)?,
            gff::feature::record::Strand::Forward,
        ))
    }

    #[test]
    fn test_reference_sequence_name() -> Result<(), noodles::core::position::TryFromIntError> {
        let feature = build_feature()?;
        assert_eq!(feature.reference_sequence_id(), 0);
        Ok(())
    }

    #[test]
    fn test_start() -> Result<(), noodles::core::position::TryFromIntError> {
        let feature = build_feature()?;
        assert_eq!(feature.start(), Position::try_from(8)?);
        Ok(())
    }

    #[test]
    fn test_end() -> Result<(), noodles::core::position::TryFromIntError> {
        let feature = build_feature()?;
        assert_eq!(feature.end(), Position::try_from(13)?);
        Ok(())
    }

    #[test]
    fn test_strand() -> Result<(), noodles::core::position::TryFromIntError> {
        let feature = build_feature()?;
        assert_eq!(feature.strand(), gff::feature::record::Strand::Forward);
        Ok(())
    }

    #[test]
    fn test_length() -> Result<(), noodles::core::position::TryFromIntError> {
        let feature = build_feature()?;
        assert_eq!(feature.length(), 6);
        Ok(())
    }
}
