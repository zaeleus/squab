use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};

use flate2::read::MultiGzDecoder;

pub fn open<P>(src: P) -> io::Result<noodles_gff::Reader<Box<dyn Read>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let extension = path.extension();
    let file = File::open(path)?;

    match extension.and_then(|ext| ext.to_str()) {
        Some("gz") => {
            let decoder = MultiGzDecoder::new(file);
            Ok(noodles_gff::Reader::new(Box::new(decoder)))
        }
        _ => Ok(noodles_gff::Reader::new(Box::new(file))),
    }
}
