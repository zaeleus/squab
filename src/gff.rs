use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
};

use flate2::read::MultiGzDecoder;

pub fn open<P>(src: P) -> io::Result<noodles_gff::Reader<Box<dyn BufRead>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let extension = path.extension();
    let file = File::open(path)?;

    match extension.and_then(|ext| ext.to_str()) {
        Some("gz") => {
            let decoder = MultiGzDecoder::new(file);
            let reader = BufReader::new(decoder);
            Ok(noodles_gff::Reader::new(Box::new(reader)))
        }
        _ => {
            let reader = BufReader::new(file);
            Ok(noodles_gff::Reader::new(Box::new(reader)))
        }
    }
}
