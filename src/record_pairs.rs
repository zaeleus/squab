use std::collections::HashMap;
use std::collections::hash_map::Drain;
use std::io::Read;

use noodles::formats::bam::{self, ByteRecord, Flag};

type RecordKey = (Vec<u8>, u8, i32, i32, i32, i32, i32);

pub struct RecordPairs<R: Read> {
    reader: bam::Reader<R>,
    // record buffer
    record: ByteRecord,
    // name-record pairs
    buf: HashMap<RecordKey, ByteRecord>,
}

impl<R: Read> RecordPairs<R> {
    pub fn new(reader: bam::Reader<R>) -> RecordPairs<R> {
        RecordPairs {
            reader,
            record: ByteRecord::new(),
            buf: HashMap::new(),
        }
    }

    fn next_pair(&mut self) -> Option<(ByteRecord, ByteRecord)> {
        loop {
            match self.reader.read_byte_record(&mut self.record) {
                Ok(0) => {
                    if !self.buf.is_empty() {
                        warn!("{} mate records missing", self.buf.len());
                    }

                    return None;
                },
                Ok(_) => {},
                Err(e) => panic!("{}", e),
            }

            let mate_key = mate_key(&self.record);

            if let Some(mate) = self.buf.remove(&mate_key) {
                let flag = Flag::new(self.record.flag());

                if flag.is_read_1() {
                    return Some((self.record.clone(), mate));
                } else if flag.is_read_2() {
                    return Some((mate, self.record.clone()));
                } else {
                    panic!("pair: read flag not set");
                }
            }

            let key = key(&self.record);
            self.buf.insert(key, self.record.clone());
        }
    }

    pub fn orphan_pairs(&mut self) -> OrphanRecordPairs {
        OrphanRecordPairs { drain: self.buf.drain() }
    }
}

impl<R: Read> Iterator for RecordPairs<R> {
    type Item = (ByteRecord, ByteRecord);

    fn next(&mut self) -> Option<Self::Item> {
        self.next_pair()
    }
}

fn key(record: &ByteRecord) -> RecordKey {
    let flag = Flag::new(record.flag());

    (
        record.read_name().to_vec(),
        if flag.is_read_1() { 1 } else { 2 },
        record.ref_id(),
        record.pos(),
        record.next_ref_id(),
        record.next_pos(),
        record.tlen(),
    )
}

fn mate_key(record: &ByteRecord) -> RecordKey {
    let flag = Flag::new(record.flag());

    (
        record.read_name().to_vec(),
        if flag.is_read_1() { 2 } else { 1 },
        record.next_ref_id(),
        record.next_pos(),
        record.ref_id(),
        record.pos(),
        -record.tlen(),
    )
}

pub struct OrphanRecordPairs<'a> {
    drain: Drain<'a, RecordKey, ByteRecord>,
}

impl<'a> Iterator for OrphanRecordPairs<'a> {
    type Item = (Option<ByteRecord>, Option<ByteRecord>);

    fn next(&mut self) -> Option<Self::Item> {
        self.drain.next().map(|(_, r)| {
            let f = Flag::new(r.flag());

            if f.is_read_1() {
                (Some(r.clone()), None)
            } else if f.is_read_2() {
                (None, Some(r.clone()))
            } else {
                panic!("orphan: read flag not set");
            }
        })
    }
}
