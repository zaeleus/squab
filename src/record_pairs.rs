use std::{
    collections::{hash_map::Drain, HashMap},
    convert::TryFrom,
    io,
};

use log::warn;
use noodles_bam as bam;
use noodles_sam as sam;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum PairPosition {
    First,
    Second,
}

impl PairPosition {
    pub fn mate(self) -> PairPosition {
        match self {
            PairPosition::First => PairPosition::Second,
            PairPosition::Second => PairPosition::First,
        }
    }
}

impl<'a> TryFrom<&'a bam::Record> for PairPosition {
    type Error = ();

    fn try_from(record: &bam::Record) -> Result<Self, Self::Error> {
        Self::try_from(record.flags())
    }
}

impl TryFrom<sam::Flags> for PairPosition {
    type Error = ();

    fn try_from(flags: sam::Flags) -> Result<Self, Self::Error> {
        if flags.is_read_1() {
            Ok(PairPosition::First)
        } else if flags.is_read_2() {
            Ok(PairPosition::Second)
        } else {
            Err(())
        }
    }
}

#[cfg(test)]
mod pair_position_tests {
    use std::convert::TryFrom;

    use noodles_sam as sam;

    use super::PairPosition;

    #[test]
    fn test_mate() {
        assert_eq!(PairPosition::First.mate(), PairPosition::Second);
        assert_eq!(PairPosition::Second.mate(), PairPosition::First);
    }

    #[test]
    fn test_try_from_flag() {
        let flags = sam::Flags::from(0x41);
        assert_eq!(PairPosition::try_from(flags), Ok(PairPosition::First));

        let flags = sam::Flags::from(0x81);
        assert_eq!(PairPosition::try_from(flags), Ok(PairPosition::Second));

        let flags = sam::Flags::from(0x01);
        assert!(PairPosition::try_from(flags).is_err());
    }
}

type RecordKey = (Vec<u8>, PairPosition, i32, i32, i32, i32, i32);

pub struct RecordPairs<R: Iterator<Item = io::Result<bam::Record>>> {
    records: R,
    buf: HashMap<RecordKey, bam::Record>,
    primary_only: bool,
}

impl<R> RecordPairs<R>
where
    R: Iterator<Item = io::Result<bam::Record>>,
{
    pub fn new(records: R, primary_only: bool) -> RecordPairs<R> {
        RecordPairs {
            records,
            buf: HashMap::new(),
            primary_only,
        }
    }

    fn next_pair(&mut self) -> Option<io::Result<(bam::Record, bam::Record)>> {
        loop {
            let record = match self.records.next() {
                Some(result) => match result {
                    Ok(r) => r,
                    Err(e) => return Some(Err(e)),
                },
                None => {
                    if !self.buf.is_empty() {
                        warn!("{} records are singletons", self.buf.len());
                    }

                    return None;
                }
            };

            if self.primary_only && is_primary(&record) {
                continue;
            }

            let mate_key = mate_key(&record);

            if let Some(mate) = self.buf.remove(&mate_key) {
                return match mate_key.1 {
                    PairPosition::First => Some(Ok((mate, record))),
                    PairPosition::Second => Some(Ok((record, mate))),
                };
            }

            let key = key(&record);

            self.buf.insert(key, record.clone());
        }
    }

    pub fn singletons(&mut self) -> Singletons {
        Singletons {
            drain: self.buf.drain(),
        }
    }
}

impl<R> Iterator for RecordPairs<R>
where
    R: Iterator<Item = io::Result<bam::Record>>,
{
    type Item = io::Result<(bam::Record, bam::Record)>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_pair()
    }
}

fn is_primary(record: &bam::Record) -> bool {
    let flags = record.flags();
    flags.is_secondary() || flags.is_supplementary()
}

fn key(record: &bam::Record) -> RecordKey {
    (
        record.name().to_vec(),
        PairPosition::try_from(record).unwrap(),
        record.reference_sequence_id(),
        record.position(),
        record.mate_reference_sequence_id(),
        record.mate_position(),
        record.template_len(),
    )
}

fn mate_key(record: &bam::Record) -> RecordKey {
    (
        record.name().to_vec(),
        PairPosition::try_from(record).map(|p| p.mate()).unwrap(),
        record.mate_reference_sequence_id(),
        record.mate_position(),
        record.reference_sequence_id(),
        record.position(),
        -record.template_len(),
    )
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
