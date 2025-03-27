use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
};

use flate2::read::MultiGzDecoder;

pub fn open<P>(src: P) -> io::Result<noodles::gff::Reader<Box<dyn BufRead>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = File::open(path)?;

    let inner: Box<dyn BufRead> = if is_gzip(path) {
        let decoder = MultiGzDecoder::new(file);
        Box::new(BufReader::new(decoder))
    } else {
        Box::new(BufReader::new(file))
    };

    Ok(noodles::gff::Reader::new(Box::new(inner)))
}

fn is_gzip<P>(src: P) -> bool
where
    P: AsRef<Path>,
{
    src.as_ref()
        .extension()
        .map(|ext| ext == "gz")
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_gzip() {
        assert!(is_gzip("in.txt.gz"));
        assert!(!is_gzip("in.txt"));
    }
}
