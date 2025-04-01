pub mod segment_position;

pub use self::segment_position::SegmentPosition;

use std::io::{self, Read};

use noodles::bam;
use rustc_hash::FxHashMap;

pub struct RecordPairs<R>
where
    R: Read,
{
    reader: bam::io::Reader<R>,
    cache: FxHashMap<Vec<u8>, Vec<bam::Record>>,
    primary_only: bool,
}

impl<R> RecordPairs<R>
where
    R: Read,
{
    pub fn new(reader: bam::io::Reader<R>, primary_only: bool) -> Self {
        Self {
            reader,
            cache: FxHashMap::default(),
            primary_only,
        }
    }

    pub fn next_pair(&mut self) -> io::Result<Option<(bam::Record, bam::Record)>> {
        use std::collections::hash_map::Entry;

        loop {
            let mut record = bam::Record::default();

            if self.reader.read_record(&mut record)? == 0 {
                return Ok(None);
            }

            if self.primary_only && is_not_primary(&record)? {
                continue;
            }

            let name = record.name().unwrap();

            match self.cache.entry(name.to_vec()) {
                Entry::Occupied(mut entry) => {
                    let records = entry.get_mut();

                    let Some(i) = find_mate(records, &record)? else {
                        records.push(record);
                        continue;
                    };

                    let mate = records.swap_remove(i);

                    if records.is_empty() {
                        entry.remove();
                    }

                    let segment_position = SegmentPosition::try_from(record.flags())
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                    return match segment_position {
                        SegmentPosition::First => Ok(Some((record, mate))),
                        SegmentPosition::Last => Ok(Some((mate, record))),
                    };
                }
                Entry::Vacant(entry) => {
                    entry.insert(vec![record]);
                }
            }
        }
    }

    pub fn unmatched_records(self) -> impl Iterator<Item = bam::Record> {
        self.cache.into_values().flatten()
    }
}

impl<R> Iterator for RecordPairs<R>
where
    R: Read,
{
    type Item = io::Result<(bam::Record, bam::Record)>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_pair() {
            Ok(Some(segments)) => Some(Ok(segments)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

fn is_not_primary(record: &bam::Record) -> io::Result<bool> {
    let flags = record.flags();
    Ok(flags.is_secondary() || flags.is_supplementary())
}

fn find_mate(records: &[bam::Record], record: &bam::Record) -> io::Result<Option<usize>> {
    for (i, mate) in records.iter().enumerate() {
        if is_mate(record, mate)? {
            return Ok(Some(i));
        }
    }

    Ok(None)
}

fn is_mate(a: &bam::Record, b: &bam::Record) -> io::Result<bool> {
    let a_fields = (
        SegmentPosition::try_from(a.flags())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        a.reference_sequence_id().transpose()?,
        a.alignment_start().transpose()?,
        a.mate_reference_sequence_id().transpose()?,
        a.mate_alignment_start().transpose()?,
        a.template_length(),
    );

    let b_fields = (
        SegmentPosition::try_from(b.flags())
            .map(|p| p.mate())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        b.mate_reference_sequence_id().transpose()?,
        b.mate_alignment_start().transpose()?,
        b.reference_sequence_id().transpose()?,
        b.alignment_start().transpose()?,
        -b.template_length(),
    );

    Ok(a_fields == b_fields)
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use bstr::BString;
    use noodles::{
        core::Position,
        sam::{
            self,
            alignment::{io::Write, record::Flags},
            header::record::value::{Map, map::ReferenceSequence},
        },
    };

    use super::*;

    fn build_record_pair() -> Result<(bam::Record, bam::Record), Box<dyn std::error::Error>> {
        let name = BString::from(b"r0");
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
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(NonZero::<usize>::MIN))
            .add_reference_sequence("sq1", Map::<ReferenceSequence>::new(NonZero::<usize>::MIN))
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
    fn test_is_mate() -> Result<(), Box<dyn std::error::Error>> {
        let (r1, r2) = build_record_pair()?;
        assert!(is_mate(&r1, &r2)?);
        Ok(())
    }
}
