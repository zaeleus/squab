use noodles::{core::Position, gff};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Feature {
    reference_sequence_name: String,
    start: Position,
    end: Position,
    strand: gff::record::Strand,
}

impl Feature {
    pub fn new(
        reference_sequence_name: String,
        start: Position,
        end: Position,
        strand: gff::record::Strand,
    ) -> Self {
        Self {
            reference_sequence_name,
            start,
            end,
            strand,
        }
    }

    pub fn reference_sequence_name(&self) -> &str {
        &self.reference_sequence_name
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

    pub fn strand(&self) -> gff::record::Strand {
        self.strand
    }

    pub fn len(&self) -> usize {
        usize::from(self.end) - usize::from(self.start) + 1
    }

    pub fn is_empty(&self) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_feature() -> Result<Feature, noodles::core::position::TryFromIntError> {
        Ok(Feature::new(
            String::from("sq0"),
            Position::try_from(8)?,
            Position::try_from(13)?,
            gff::record::Strand::Forward,
        ))
    }

    #[test]
    fn test_reference_sequence_name() -> Result<(), noodles::core::position::TryFromIntError> {
        let feature = build_feature()?;
        assert_eq!(feature.reference_sequence_name(), "sq0");
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
        assert_eq!(feature.strand(), gff::record::Strand::Forward);
        Ok(())
    }

    #[test]
    fn test_len() -> Result<(), noodles::core::position::TryFromIntError> {
        let feature = build_feature()?;
        assert_eq!(feature.len(), 6);
        Ok(())
    }

    #[test]
    fn test_is_empty() -> Result<(), noodles::core::position::TryFromIntError> {
        let feature = Feature::new(
            String::from("sq0"),
            Position::MIN,
            Position::MIN,
            gff::record::Strand::Forward,
        );
        assert!(!feature.is_empty());

        let feature = build_feature()?;
        assert!(!feature.is_empty());

        Ok(())
    }
}
