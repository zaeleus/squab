mod pair_position;

pub use self::pair_position::PairPosition;

use std::{
    collections::{hash_map::Drain, HashMap},
    io::{self, Read},
};

use noodles::{bam, core::Position, sam};
use tracing::warn;

use crate::Record;

type RecordKey = (
    Option<sam::record::ReadName>,
    PairPosition,
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
    buf: HashMap<RecordKey, Record>,
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

    pub fn next_pair(&mut self) -> io::Result<Option<(Record, Record)>> {
        loop {
            let mut lazy_record = bam::lazy::Record::default();

            match self.reader.read_lazy_record(&mut lazy_record) {
                Ok(0) => {
                    if !self.buf.is_empty() {
                        warn!("{} records are singletons", self.buf.len());
                    }

                    return Ok(None);
                }
                Ok(_) => {}
                Err(e) => return Err(e),
            }

            let record = Record::try_from(&lazy_record)?;

            if self.primary_only && is_not_primary(&record) {
                continue;
            }

            let mate_key = mate_key(&record)?;

            if let Some(mate) = self.buf.remove(&mate_key) {
                return match mate_key.1 {
                    PairPosition::First => Ok(Some((mate, record))),
                    PairPosition::Second => Ok(Some((record, mate))),
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

fn is_not_primary(record: &Record) -> bool {
    let flags = record.flags();
    flags.is_secondary() || flags.is_supplementary()
}

fn key(record: &Record) -> io::Result<RecordKey> {
    Ok((
        record.read_name().cloned(),
        PairPosition::try_from(record.flags())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        record.reference_sequence_id(),
        record.alignment_start(),
        record.mate_reference_sequence_id(),
        record.mate_alignment_start(),
        record.template_length(),
    ))
}

fn mate_key(record: &Record) -> io::Result<RecordKey> {
    Ok((
        record.read_name().cloned(),
        PairPosition::try_from(record.flags())
            .map(|p| p.mate())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        record.mate_reference_sequence_id(),
        record.mate_alignment_start(),
        record.reference_sequence_id(),
        record.alignment_start(),
        -record.template_length(),
    ))
}

pub struct Singletons<'a> {
    drain: Drain<'a, RecordKey, Record>,
}

impl<'a> Iterator for Singletons<'a> {
    type Item = Record;

    fn next(&mut self) -> Option<Self::Item> {
        self.drain.next().map(|(_, r)| r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_record_pair() -> Result<(Record, Record), Box<dyn std::error::Error>> {
        use noodles::sam::record::{Flags, ReadName};

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

        Ok((Record::from(r1), Record::from(r2)))
    }

    #[test]
    fn test_key() -> Result<(), Box<dyn std::error::Error>> {
        let (r1, _) = build_record_pair()?;

        let actual = key(&r1)?;
        let expected = (
            r1.read_name().cloned(),
            PairPosition::First,
            r1.reference_sequence_id(),
            r1.alignment_start(),
            r1.mate_reference_sequence_id(),
            r1.mate_alignment_start(),
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
            r1.read_name().cloned(),
            PairPosition::Second,
            r1.mate_reference_sequence_id(),
            r1.mate_alignment_start(),
            r1.reference_sequence_id(),
            r1.alignment_start(),
            -r1.template_length(),
        );

        assert_eq!(actual, expected);

        Ok(())
    }
}
