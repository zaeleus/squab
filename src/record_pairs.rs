pub mod segment_position;

pub use self::segment_position::SegmentPosition;

use std::{
    collections::{hash_map::Drain, HashMap},
    io::{self, Read},
};

use noodles::{bam, core::Position};
use tracing::warn;

type RecordKey = (
    Option<Vec<u8>>,
    SegmentPosition,
    Option<usize>,
    Option<Position>,
    Option<usize>,
    Option<Position>,
    i32,
);

pub struct RecordPairs<R>
where
    R: Read,
{
    reader: bam::Reader<R>,
    buf: HashMap<RecordKey, bam::lazy::Record>,
    primary_only: bool,
}

impl<R> RecordPairs<R>
where
    R: Read,
{
    pub fn new(reader: bam::Reader<R>, primary_only: bool) -> Self {
        Self {
            reader,
            buf: HashMap::new(),
            primary_only,
        }
    }

    pub fn next_pair(&mut self) -> io::Result<Option<(bam::lazy::Record, bam::lazy::Record)>> {
        loop {
            let mut record = bam::lazy::Record::default();

            if self.reader.read_lazy_record(&mut record)? == 0 {
                if !self.buf.is_empty() {
                    warn!("{} records are singletons", self.buf.len());
                }

                return Ok(None);
            }

            if self.primary_only && is_not_primary(&record)? {
                continue;
            }

            let mate_key = mate_key(&record)?;

            if let Some(mate) = self.buf.remove(&mate_key) {
                return match mate_key.1 {
                    SegmentPosition::First => Ok(Some((mate, record))),
                    SegmentPosition::Last => Ok(Some((record, mate))),
                };
            }

            let key = key(&record)?;

            self.buf.insert(key, record);
        }
    }

    pub fn singletons(&mut self) -> Singletons {
        Singletons {
            drain: self.buf.drain(),
        }
    }
}

fn is_not_primary(record: &bam::lazy::Record) -> io::Result<bool> {
    let flags = record.flags()?;
    Ok(flags.is_secondary() || flags.is_supplementary())
}

fn key(record: &bam::lazy::Record) -> io::Result<RecordKey> {
    Ok((
        record.read_name().map(|buf| buf.as_ref().to_vec()),
        SegmentPosition::try_from(record.flags()?)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        record.reference_sequence_id()?,
        record.alignment_start()?,
        record.mate_reference_sequence_id()?,
        record.mate_alignment_start()?,
        record.template_length(),
    ))
}

fn mate_key(record: &bam::lazy::Record) -> io::Result<RecordKey> {
    Ok((
        record.read_name().map(|buf| buf.as_ref().to_vec()),
        SegmentPosition::try_from(record.flags()?)
            .map(|p| p.mate())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        record.mate_reference_sequence_id()?,
        record.mate_alignment_start()?,
        record.reference_sequence_id()?,
        record.alignment_start()?,
        -record.template_length(),
    ))
}

pub struct Singletons<'a> {
    drain: Drain<'a, RecordKey, bam::lazy::Record>,
}

impl<'a> Iterator for Singletons<'a> {
    type Item = bam::lazy::Record;

    fn next(&mut self) -> Option<Self::Item> {
        self.drain.next().map(|(_, r)| r)
    }
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use noodles::sam::{
        self,
        header::record::value::{map::ReferenceSequence, Map},
        record::{Flags, ReadName},
    };

    use super::*;

    fn build_record_pair(
    ) -> Result<(bam::lazy::Record, bam::lazy::Record), Box<dyn std::error::Error>> {
        let read_name: ReadName = "r0".parse()?;
        let reference_sequence_id = 0;
        let alignment_start = Position::try_from(8)?;
        let mate_reference_sequence_id = 1;
        let mate_position = Position::try_from(13)?;

        let r1 = sam::alignment::Record::builder()
            .set_read_name(read_name.clone())
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .set_reference_sequence_id(reference_sequence_id)
            .set_alignment_start(alignment_start)
            .set_mate_reference_sequence_id(mate_reference_sequence_id)
            .set_mate_alignment_start(mate_position)
            .set_template_length(144)
            .build();

        let r2 = sam::alignment::Record::builder()
            .set_read_name(read_name)
            .set_flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .set_reference_sequence_id(mate_reference_sequence_id)
            .set_alignment_start(mate_position)
            .set_mate_reference_sequence_id(reference_sequence_id)
            .set_mate_alignment_start(alignment_start)
            .set_template_length(-144)
            .build();

        let header = sam::Header::builder()
            .add_reference_sequence(
                "sq0".parse()?,
                Map::<ReferenceSequence>::new(NonZeroUsize::MIN),
            )
            .add_reference_sequence(
                "sq1".parse()?,
                Map::<ReferenceSequence>::new(NonZeroUsize::MIN),
            )
            .build();

        let mut writer = bam::Writer::from(Vec::new());
        writer.write_record(&header, &r1)?;
        let buf1 = writer.into_inner()[4..].to_vec();

        let mut writer = bam::Writer::from(Vec::new());
        writer.write_record(&header, &r2)?;
        let buf2 = writer.into_inner()[4..].to_vec();

        let lr1 = bam::lazy::Record::try_from(buf1)?;
        let lr2 = bam::lazy::Record::try_from(buf2)?;

        Ok((lr1, lr2))
    }

    #[test]
    fn test_key() -> Result<(), Box<dyn std::error::Error>> {
        let (r1, _) = build_record_pair()?;

        let actual = key(&r1)?;
        let expected = (
            r1.read_name().map(|buf| buf.as_ref().to_vec()),
            SegmentPosition::First,
            r1.reference_sequence_id()?,
            r1.alignment_start()?,
            r1.mate_reference_sequence_id()?,
            r1.mate_alignment_start()?,
            r1.template_length(),
        );

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_mate_key() -> Result<(), Box<dyn std::error::Error>> {
        let (r1, _) = build_record_pair()?;

        let actual = mate_key(&r1)?;
        let expected = (
            r1.read_name().map(|buf| buf.as_ref().to_vec()),
            SegmentPosition::Last,
            r1.mate_reference_sequence_id()?,
            r1.mate_alignment_start()?,
            r1.reference_sequence_id()?,
            r1.alignment_start()?,
            -r1.template_length(),
        );

        assert_eq!(actual, expected);

        Ok(())
    }
}
