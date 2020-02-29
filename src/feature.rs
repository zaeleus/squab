use noodles_gff as gff;

#[derive(Clone, Debug)]
pub struct Feature {
    reference_name: String,
    start: u64,
    end: u64,
    strand: gff::Strand,
}

impl Feature {
    pub fn new(reference_name: String, start: u64, end: u64, strand: gff::Strand) -> Self {
        Self {
            reference_name,
            start,
            end,
            strand,
        }
    }

    pub fn reference_name(&self) -> &str {
        &self.reference_name
    }

    pub fn start(&self) -> u64 {
        self.start
    }

    pub fn end(&self) -> u64 {
        self.end
    }

    pub fn end_mut(&mut self) -> &mut u64 {
        &mut self.end
    }

    pub fn strand(&self) -> gff::Strand {
        self.strand
    }

    pub fn len(&self) -> u64 {
        self.end - self.start + 1
    }
}
