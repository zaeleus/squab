use std::io;

use noodles::{
    bam,
    core::Position,
    sam::{
        self,
        record::{Cigar, Data, Flags, MappingQuality, ReadName},
    },
};

pub struct Record {
    read_name: Option<ReadName>,
    flags: Flags,
    reference_sequence_id: Option<usize>,
    alignment_start: Option<Position>,
    mapping_quality: Option<MappingQuality>,
    cigar: Cigar,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<Position>,
    template_length: i32,
    data: Data,
}

impl Record {
    pub fn read_name(&self) -> Option<&ReadName> {
        self.read_name.as_ref()
    }

    pub fn flags(&self) -> Flags {
        self.flags
    }

    pub fn reference_sequence_id(&self) -> Option<usize> {
        self.reference_sequence_id
    }

    pub fn alignment_start(&self) -> Option<Position> {
        self.alignment_start
    }

    pub fn mapping_quality(&self) -> Option<MappingQuality> {
        self.mapping_quality
    }

    pub fn cigar(&self) -> &Cigar {
        &self.cigar
    }

    pub fn mate_reference_sequence_id(&self) -> Option<usize> {
        self.mate_reference_sequence_id
    }

    pub fn mate_alignment_start(&self) -> Option<Position> {
        self.mate_alignment_start
    }

    pub fn template_length(&self) -> i32 {
        self.template_length
    }

    pub fn data(&self) -> &Data {
        &self.data
    }
}

impl From<sam::alignment::Record> for Record {
    fn from(record: sam::alignment::Record) -> Self {
        Self {
            read_name: record.read_name().cloned(),
            flags: record.flags(),
            reference_sequence_id: record.reference_sequence_id(),
            alignment_start: record.alignment_start(),
            mapping_quality: record.mapping_quality(),
            cigar: record.cigar().clone(),
            mate_reference_sequence_id: record.mate_reference_sequence_id(),
            mate_alignment_start: record.mate_alignment_start(),
            template_length: record.template_length(),
            data: record.data().clone(),
        }
    }
}

impl TryFrom<&bam::lazy::Record> for Record {
    type Error = io::Error;

    fn try_from(record: &bam::lazy::Record) -> Result<Self, Self::Error> {
        Ok(Self {
            read_name: record
                .read_name()
                .map(|read_name| read_name.try_into())
                .transpose()?,
            flags: record.flags()?,
            reference_sequence_id: record.reference_sequence_id()?,
            alignment_start: record.alignment_start()?,
            mapping_quality: record.mapping_quality(),
            cigar: sam::record::Cigar::try_from(record.cigar())?,
            mate_reference_sequence_id: record.mate_reference_sequence_id()?,
            mate_alignment_start: record.mate_alignment_start()?,
            template_length: record.template_length(),
            data: sam::record::Data::try_from(record.data())?,
        })
    }
}

pub fn alignment_end(alignment_start: Option<Position>, cigar: &Cigar) -> Option<Position> {
    alignment_start.and_then(|start| {
        let span = alignment_span(cigar);
        let end = usize::from(start) + span - 1;
        Position::new(end)
    })
}

fn alignment_span(cigar: &Cigar) -> usize {
    cigar.alignment_span()
}
