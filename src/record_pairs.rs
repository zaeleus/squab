use std::collections::hash_map::Drain;
use std::collections::HashMap;
use std::io::{self, Read};

use log::warn;
use noodles::formats::bam::{self, Flag, Record};

#[derive(Debug, Eq, Hash, PartialEq)]
pub enum PairPosition {
    First,
    Second,
}

impl PairPosition {
    pub fn mate(&self) -> PairPosition {
        match *self {
            PairPosition::First => PairPosition::Second,
            PairPosition::Second => PairPosition::First,
        }
    }
}

impl<'a> From<&'a Record> for PairPosition {
    fn from(record: &Record) -> PairPosition {
        record.flag().into()
    }
}

impl From<Flag> for PairPosition {
    fn from(flag: Flag) -> PairPosition {
        if flag.is_read_1() {
            PairPosition::First
        } else if flag.is_read_2() {
            PairPosition::Second
        } else {
            panic!("unknown pair position");
        }
    }
}

#[cfg(test)]
mod pair_position_tests {
    use noodles::formats::bam::Flag;

    use super::PairPosition;

    #[test]
    fn test_mate() {
        assert_eq!(PairPosition::First.mate(), PairPosition::Second);
        assert_eq!(PairPosition::Second.mate(), PairPosition::First);
    }

    #[test]
    fn test_from_flag() {
        let flag = Flag::from(0x41);
        assert_eq!(PairPosition::from(flag), PairPosition::First);

        let flag = Flag::from(0x81);
        assert_eq!(PairPosition::from(flag), PairPosition::Second);
    }

    #[test]
    #[should_panic]
    fn test_from_flag_with_invalid_flag() {
        let flag = Flag::from(0x01);
        PairPosition::from(flag);
    }
}

type RecordKey = (Vec<u8>, PairPosition, i32, i32, i32, i32, i32);

pub struct RecordPairs<R: Read> {
    reader: bam::Reader<R>,
    record: Record,
    buf: HashMap<RecordKey, Record>,
    primary_only: bool,
}

impl<R: Read> RecordPairs<R> {
    pub fn new(reader: bam::Reader<R>, primary_only: bool) -> RecordPairs<R> {
        RecordPairs {
            reader,
            record: Record::new(),
            buf: HashMap::new(),
            primary_only,
        }
    }

    fn next_pair(&mut self) -> Option<io::Result<(Record, Record)>> {
        loop {
            match self.reader.read_record(&mut self.record) {
                Ok(0) => {
                    if !self.buf.is_empty() {
                        warn!("{} records are singletons", self.buf.len());
                    }

                    return None;
                }
                Ok(_) => {}
                Err(e) => return Some(Err(e)),
            }

            if self.primary_only && is_primary(&self.record) {
                continue;
            }

            let mate_key = mate_key(&self.record);

            if let Some(mate) = self.buf.remove(&mate_key) {
                return match mate_key.1 {
                    PairPosition::First => Some(Ok((mate, self.record.clone()))),
                    PairPosition::Second => Some(Ok((self.record.clone(), mate))),
                };
            }

            let key = key(&self.record);

            self.buf.insert(key, self.record.clone());
        }
    }

    pub fn singletons(&mut self) -> Singletons {
        Singletons {
            drain: self.buf.drain(),
        }
    }
}

impl<R: Read> Iterator for RecordPairs<R> {
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
        PairPosition::from(record),
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
        PairPosition::from(record).mate(),
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
