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

    const START: Position = Position::new(8).unwrap();
    const END: Position = Position::new(13).unwrap();

    fn build_feature() -> Feature {
        Feature::new(0, START, END, gff::feature::record::Strand::Forward)
    }

    #[test]
    fn test_reference_sequence_id() {
        let feature = build_feature();
        assert_eq!(feature.reference_sequence_id(), 0);
    }

    #[test]
    fn test_start() {
        let feature = build_feature();
        assert_eq!(feature.start(), START);
    }

    #[test]
    fn test_end() {
        let feature = build_feature();
        assert_eq!(feature.end(), END);
    }

    #[test]
    fn test_strand() {
        let feature = build_feature();
        assert_eq!(feature.strand(), gff::feature::record::Strand::Forward);
    }

    #[test]
    fn test_length() {
        let feature = build_feature();
        assert_eq!(feature.length(), 6);
    }
}
