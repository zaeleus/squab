use std::{
    collections::{hash_map::Drain, HashMap},
    convert::TryFrom,
    io,
};

use log::warn;
use noodles_bam::{Flag, Record};

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

impl<'a> TryFrom<&'a Record> for PairPosition {
    type Error = ();

    fn try_from(record: &Record) -> Result<Self, Self::Error> {
        Self::try_from(record.flag())
    }
}

impl TryFrom<Flag> for PairPosition {
    type Error = ();

    fn try_from(flag: Flag) -> Result<Self, Self::Error> {
        if flag.is_read_1() {
            Ok(PairPosition::First)
        } else if flag.is_read_2() {
            Ok(PairPosition::Second)
        } else {
            Err(())
        }
    }
}

#[cfg(test)]
mod pair_position_tests {
    use std::convert::TryFrom;

    use noodles_bam::Flag;

    use super::PairPosition;

    #[test]
    fn test_mate() {
        assert_eq!(PairPosition::First.mate(), PairPosition::Second);
        assert_eq!(PairPosition::Second.mate(), PairPosition::First);
    }

    #[test]
    fn test_try_from_flag() {
        let flag = Flag::from(0x41);
        assert_eq!(PairPosition::try_from(flag), Ok(PairPosition::First));

        let flag = Flag::from(0x81);
        assert_eq!(PairPosition::try_from(flag), Ok(PairPosition::Second));

        let flag = Flag::from(0x01);
        assert!(PairPosition::try_from(flag).is_err());
    }
}

type RecordKey = (Vec<u8>, PairPosition, i32, i32, i32, i32, i32);

pub struct RecordPairs<R: Iterator<Item = io::Result<Record>>> {
    records: R,
    buf: HashMap<RecordKey, Record>,
    primary_only: bool,
}

impl<R> RecordPairs<R>
where
    R: Iterator<Item = io::Result<Record>>,
{
    pub fn new(records: R, primary_only: bool) -> RecordPairs<R> {
        RecordPairs {
            records,
            buf: HashMap::new(),
            primary_only,
        }
    }

    fn next_pair(&mut self) -> Option<io::Result<(Record, Record)>> {
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
    R: Iterator<Item = io::Result<Record>>,
{
    type Item = io::Result<(Record, Record)>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_pair()
    }
}

fn is_primary(record: &Record) -> bool {
    let flag = record.flag();
    flag.is_secondary() || flag.is_supplementary()
}

fn key(record: &Record) -> RecordKey {
    (
        record.read_name().to_vec(),
        PairPosition::try_from(record).unwrap(),
        record.ref_id(),
        record.pos(),
        record.next_ref_id(),
        record.next_pos(),
        record.tlen(),
    )
}

fn mate_key(record: &Record) -> RecordKey {
    (
        record.read_name().to_vec(),
        PairPosition::try_from(record).map(|p| p.mate()).unwrap(),
        record.next_ref_id(),
        record.next_pos(),
        record.ref_id(),
        record.pos(),
        -record.tlen(),
    )
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
