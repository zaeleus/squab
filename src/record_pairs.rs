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
    reader: bam::io::Reader<R>,
    buf: HashMap<RecordKey, bam::Record>,
    primary_only: bool,
}

impl<R> RecordPairs<R>
where
    R: Read,
{
    pub fn new(reader: bam::io::Reader<R>, primary_only: bool) -> Self {
        Self {
            reader,
            buf: HashMap::new(),
            primary_only,
        }
    }

    pub fn next_pair(&mut self) -> io::Result<Option<(bam::Record, bam::Record)>> {
        loop {
            let mut record = bam::Record::default();

            if self.reader.read_record(&mut record)? == 0 {
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

fn is_not_primary(record: &bam::Record) -> io::Result<bool> {
    let flags = record.flags();
    Ok(flags.is_secondary() || flags.is_supplementary())
}

fn key(record: &bam::Record) -> io::Result<RecordKey> {
    Ok((
        record.name().map(|name| name.as_bytes().to_vec()),
        SegmentPosition::try_from(record.flags())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        record.reference_sequence_id().transpose()?,
        record.alignment_start().transpose()?,
        record.mate_reference_sequence_id().transpose()?,
        record.mate_alignment_start().transpose()?,
        record.template_length(),
    ))
}

fn mate_key(record: &bam::Record) -> io::Result<RecordKey> {
    Ok((
        record.name().map(|name| name.as_bytes().to_vec()),
        SegmentPosition::try_from(record.flags())
            .map(|p| p.mate())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        record.mate_reference_sequence_id().transpose()?,
        record.mate_alignment_start().transpose()?,
        record.reference_sequence_id().transpose()?,
        record.alignment_start().transpose()?,
        -record.template_length(),
    ))
}

pub struct Singletons<'a> {
    drain: Drain<'a, RecordKey, bam::Record>,
}

impl<'a> Iterator for Singletons<'a> {
    type Item = bam::Record;

    fn next(&mut self) -> Option<Self::Item> {
        self.drain.next().map(|(_, r)| r)
    }
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use noodles::sam::{
        self,
        alignment::{io::Write, record::Flags, record_buf::Name},
        header::record::value::{map::ReferenceSequence, Map},
    };

    use super::*;

    fn build_record_pair() -> Result<(bam::Record, bam::Record), Box<dyn std::error::Error>> {
        let name = Name::from(b"r0");
        let reference_sequence_id = 0;
        let alignment_start = Position::try_from(8)?;
        let mate_reference_sequence_id = 1;
        let mate_position = Position::try_from(13)?;

        let rb1 = sam::alignment::RecordBuf::builder()
            .set_name(name.clone())
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .set_reference_sequence_id(reference_sequence_id)
            .set_alignment_start(alignment_start)
            .set_mate_reference_sequence_id(mate_reference_sequence_id)
            .set_mate_alignment_start(mate_position)
            .set_template_length(144)
            .build();

        let rb2 = sam::alignment::RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .set_reference_sequence_id(mate_reference_sequence_id)
            .set_alignment_start(mate_position)
            .set_mate_reference_sequence_id(reference_sequence_id)
            .set_mate_alignment_start(alignment_start)
            .set_template_length(-144)
            .build();

        let header = sam::Header::builder()
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(NonZeroUsize::MIN))
            .add_reference_sequence("sq1", Map::<ReferenceSequence>::new(NonZeroUsize::MIN))
            .build();

        let mut writer = bam::io::Writer::from(Vec::new());
        writer.write_alignment_record(&header, &rb1)?;
        writer.write_alignment_record(&header, &rb2)?;

        let mut reader = bam::io::Reader::from(writer.get_ref().as_slice());

        let mut r1 = bam::Record::default();
        reader.read_record(&mut r1)?;

        let mut r2 = bam::Record::default();
        reader.read_record(&mut r2)?;

        Ok((r1, r2))
    }

    #[test]
    fn test_key() -> Result<(), Box<dyn std::error::Error>> {
        let (r1, _) = build_record_pair()?;

        let actual = key(&r1)?;
        let expected = (
            r1.name().map(|name| name.as_bytes().to_vec()),
            SegmentPosition::First,
            r1.reference_sequence_id().transpose()?,
            r1.alignment_start().transpose()?,
            r1.mate_reference_sequence_id().transpose()?,
            r1.mate_alignment_start().transpose()?,
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
            r1.name().map(|name| name.as_bytes().to_vec()),
            SegmentPosition::Last,
            r1.mate_reference_sequence_id().transpose()?,
            r1.mate_alignment_start().transpose()?,
            r1.reference_sequence_id().transpose()?,
            r1.alignment_start().transpose()?,
            -r1.template_length(),
        );

        assert_eq!(actual, expected);

        Ok(())
    }
}
