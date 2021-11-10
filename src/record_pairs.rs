mod pair_position;

pub use self::pair_position::PairPosition;

use std::{
    collections::{hash_map::Drain, HashMap},
    io,
};

use noodles::bam;
use tokio::io::AsyncRead;
use tracing::warn;

type RecordKey = (
    Vec<u8>,
    PairPosition,
    Option<i32>,
    Option<i32>,
    Option<i32>,
    Option<i32>,
    i32,
);

pub struct RecordPairs<R>
where
    R: AsyncRead,
{
    reader: bam::AsyncReader<R>,
    buf: HashMap<RecordKey, bam::Record>,
    primary_only: bool,
}

impl<R> RecordPairs<R>
where
    R: AsyncRead + Unpin,
{
    pub fn new(reader: bam::AsyncReader<R>, primary_only: bool) -> Self {
        Self {
            reader,
            buf: HashMap::new(),
            primary_only,
        }
    }

    pub async fn next_pair(&mut self) -> io::Result<Option<(bam::Record, bam::Record)>> {
        loop {
            let mut record = bam::Record::default();

            match self.reader.read_record(&mut record).await {
                Ok(0) => {
                    if !self.buf.is_empty() {
                        warn!("{} records are singletons", self.buf.len());
                    }

                    return Ok(None);
                }
                Ok(_) => {}
                Err(e) => return Err(e),
            }

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

fn is_not_primary(record: &bam::Record) -> bool {
    let flags = record.flags();
    flags.is_secondary() || flags.is_supplementary()
}

fn key(record: &bam::Record) -> io::Result<RecordKey> {
    Ok((
        record
            .read_name()
            .map(|s| s.to_bytes().to_vec())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        PairPosition::try_from(record)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        record.reference_sequence_id().map(i32::from),
        record.position().map(i32::from),
        record.mate_reference_sequence_id().map(i32::from),
        record.mate_position().map(i32::from),
        record.template_length(),
    ))
}

fn mate_key(record: &bam::Record) -> io::Result<RecordKey> {
    Ok((
        record
            .read_name()
            .map(|s| s.to_bytes().to_vec())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        PairPosition::try_from(record)
            .map(|p| p.mate())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        record.mate_reference_sequence_id().map(i32::from),
        record.mate_position().map(i32::from),
        record.reference_sequence_id().map(i32::from),
        record.position().map(i32::from),
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
    use noodles::sam::{
        self,
        header::ReferenceSequence,
        record::{Flags, Position, ReadName, ReferenceSequenceName},
    };

    use super::*;

    fn build_record_pair() -> Result<(bam::Record, bam::Record), Box<dyn std::error::Error>> {
        let read_name: ReadName = "r0".parse()?;
        let reference_sequence_name: ReferenceSequenceName = "sq0".parse()?;
        let position = Position::try_from(8)?;
        let mate_reference_sequence_name: ReferenceSequenceName = "sq1".parse()?;
        let mate_position = Position::try_from(13)?;

        let reference_sequences = vec![("sq0", 8), ("sq1", 13)]
            .into_iter()
            .map(|(name, len)| ReferenceSequence::new(name, len).map(|rs| (name.into(), rs)))
            .into_iter()
            .collect::<Result<_, _>>()?;

        let s1 = sam::Record::builder()
            .set_read_name(read_name.clone())
            .set_flags(Flags::PAIRED | Flags::READ_1)
            .set_reference_sequence_name(reference_sequence_name.clone())
            .set_position(position)
            .set_mate_reference_sequence_name(mate_reference_sequence_name.clone())
            .set_mate_position(mate_position)
            .set_template_length(144)
            .build()?;

        let r1 = bam::Record::try_from_sam_record(&reference_sequences, &s1)?;

        let s2 = sam::Record::builder()
            .set_read_name(read_name)
            .set_flags(Flags::PAIRED | Flags::READ_2)
            .set_reference_sequence_name(mate_reference_sequence_name)
            .set_position(mate_position)
            .set_mate_reference_sequence_name(reference_sequence_name)
            .set_mate_position(position)
            .set_template_length(-144)
            .build()?;

        let r2 = bam::Record::try_from_sam_record(&reference_sequences, &s2)?;

        Ok((r1, r2))
    }

    #[test]
    fn test_key() -> Result<(), Box<dyn std::error::Error>> {
        let (r1, _) = build_record_pair()?;

        let actual = key(&r1)?;
        let expected = (
            b"r0".to_vec(),
            PairPosition::First,
            Some(0),
            Some(8),
            Some(1),
            Some(13),
            144,
        );

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_mate_key() -> Result<(), Box<dyn std::error::Error>> {
        let (r1, _) = build_record_pair()?;

        let actual = mate_key(&r1)?;
        let expected = (
            b"r0".to_vec(),
            PairPosition::Second,
            Some(1),
            Some(13),
            Some(0),
            Some(8),
            -144,
        );

        assert_eq!(actual, expected);

        Ok(())
    }
}
