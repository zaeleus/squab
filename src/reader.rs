use std::{
    collections::HashMap,
    io::{self, Read},
};

use csv::StringRecord;

const NAME_INDEX: usize = 0;
const COUNT_INDEX: usize = 1;

static HTSEQ_COUNT_META_PREFIX: &str = "__";

pub fn read_counts<R>(reader: &mut R) -> io::Result<HashMap<String, u64>>
where
    R: Read,
{
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(reader);

    let mut count_map = HashMap::new();

    for result in rdr.records() {
        let record = result?;

        let id = parse_id(&record)?;

        if id.starts_with(HTSEQ_COUNT_META_PREFIX) {
            break;
        }

        let count = parse_count(&record)?;

        count_map.insert(id.into(), count);
    }

    Ok(count_map)
}

fn parse_id(record: &StringRecord) -> io::Result<&str> {
    let cell = record.get(NAME_INDEX);

    cell.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid name: {:?}", cell),
        )
    })
}

fn parse_count(record: &StringRecord) -> io::Result<u64> {
    let cell = record.get(COUNT_INDEX);

    cell.and_then(|s| s.parse().ok()).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid count: {:?}", cell),
        )
    })
}
