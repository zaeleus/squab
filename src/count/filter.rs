use std::io;

use noodles_bam as bam;
use noodles_sam as sam;

use super::Context;

#[derive(Clone)]
pub struct Filter {
    min_mapping_quality: u8,
    with_secondary_records: bool,
    with_supplementary_records: bool,
    with_nonunique_records: bool,
}

impl Filter {
    pub fn min_mapping_quality(&self) -> u8 {
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
        min_mapping_quality: u8,
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

    pub fn filter(&self, ctx: &mut Context, record: &bam::Record) -> io::Result<bool> {
        let flags = record.flags();

        if flags.is_unmapped() {
            ctx.unmapped += 1;
            return Ok(true);
        }

        if (!self.with_secondary_records && flags.is_secondary())
            || (!self.with_supplementary_records && flags.is_supplementary())
        {
            return Ok(true);
        }

        if !self.with_nonunique_records && is_nonunique_record(&record)? {
            ctx.nonunique += 1;
            return Ok(true);
        }

        if u8::from(record.mapping_quality()) < self.min_mapping_quality {
            ctx.low_quality += 1;
            return Ok(true);
        }

        Ok(false)
    }

    pub fn filter_pair(
        &self,
        ctx: &mut Context,
        r1: &bam::Record,
        r2: &bam::Record,
    ) -> io::Result<bool> {
        let f1 = r1.flags();
        let f2 = r2.flags();

        if f1.is_unmapped() && f2.is_unmapped() {
            ctx.unmapped += 1;
            return Ok(true);
        }

        if (!self.with_secondary_records && (f1.is_secondary() || f2.is_secondary()))
            || (!self.with_supplementary_records
                && (f1.is_supplementary() || f2.is_supplementary()))
        {
            return Ok(true);
        }

        if !self.with_nonunique_records && (is_nonunique_record(&r1)? || is_nonunique_record(&r2)?)
        {
            ctx.nonunique += 1;
            return Ok(true);
        }

        if u8::from(r1.mapping_quality()) < self.min_mapping_quality
            || u8::from(r2.mapping_quality()) < self.min_mapping_quality
        {
            ctx.low_quality += 1;
            return Ok(true);
        }

        Ok(false)
    }
}

fn is_nonunique_record(record: &bam::Record) -> io::Result<bool> {
    use bam::record::data::field::Value;
    use sam::record::data::field::Tag;

    let data = record.data();

    for result in data.fields() {
        let field = result?;

        if field.tag() == &Tag::AlignmentHitCount {
            if let Value::Int8(n) = field.value() {
                return Ok(*n > 1);
            }
        }
    }

    Ok(false)
}
