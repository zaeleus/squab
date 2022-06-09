use std::io;

use noodles::sam::{self, record::MappingQuality};

use super::context::Event;
use crate::Record;

#[derive(Clone)]
pub struct Filter {
    min_mapping_quality: MappingQuality,
    with_secondary_records: bool,
    with_supplementary_records: bool,
    with_nonunique_records: bool,
}

impl Filter {
    pub fn min_mapping_quality(&self) -> MappingQuality {
        self.min_mapping_quality
    }

    pub fn with_secondary_records(&self) -> bool {
        self.with_secondary_records
    }

    pub fn with_supplementary_records(&self) -> bool {
        self.with_supplementary_records
    }

    pub fn with_nonunique_records(&self) -> bool {
        self.with_nonunique_records
    }
}

impl Filter {
    pub fn new(
        min_mapping_quality: MappingQuality,
        with_secondary_records: bool,
        with_supplementary_records: bool,
        with_nonunique_records: bool,
    ) -> Filter {
        Self {
            min_mapping_quality,
            with_secondary_records,
            with_supplementary_records,
            with_nonunique_records,
        }
    }

    pub fn filter(&self, record: &Record) -> io::Result<Option<Event>> {
        let flags = record.flags();

        if flags.is_unmapped() {
            return Ok(Some(Event::Unmapped));
        }

        if (!self.with_secondary_records && flags.is_secondary())
            || (!self.with_supplementary_records && flags.is_supplementary())
        {
            return Ok(Some(Event::Skip));
        }

        if !self.with_nonunique_records && is_nonunique_record(record)? {
            return Ok(Some(Event::Nonunique));
        }

        if let Some(mapping_quality) = record.mapping_quality() {
            if mapping_quality < self.min_mapping_quality {
                return Ok(Some(Event::LowQuality));
            }
        }

        Ok(None)
    }

    pub fn filter_pair(&self, r1: &Record, r2: &Record) -> io::Result<Option<Event>> {
        let f1 = r1.flags();
        let f2 = r2.flags();

        if f1.is_unmapped() && f2.is_unmapped() {
            return Ok(Some(Event::Unmapped));
        }

        if (!self.with_secondary_records && (f1.is_secondary() || f2.is_secondary()))
            || (!self.with_supplementary_records
                && (f1.is_supplementary() || f2.is_supplementary()))
        {
            return Ok(Some(Event::Skip));
        }

        if !self.with_nonunique_records && (is_nonunique_record(r1)? || is_nonunique_record(r2)?) {
            return Ok(Some(Event::Nonunique));
        }

        if let Some(mapping_quality) = r1.mapping_quality() {
            if mapping_quality < self.min_mapping_quality() {
                return Ok(Some(Event::LowQuality));
            }
        }

        if let Some(mapping_quality) = r2.mapping_quality() {
            if mapping_quality < self.min_mapping_quality() {
                return Ok(Some(Event::LowQuality));
            }
        }

        Ok(None)
    }
}

fn is_nonunique_record(record: &Record) -> io::Result<bool> {
    use sam::record::data::field::{value::Type, Tag};

    let field = match record.data().get(Tag::AlignmentHitCount) {
        Some(f) => f,
        None => return Ok(false),
    };

    let value = field.value();

    match value.as_int() {
        Some(hits) => Ok(hits > 1),
        None => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid {} value type: expected {:?}, got {:?}",
                Tag::AlignmentHitCount,
                Type::Int32,
                value.ty()
            ),
        )),
    }
}
