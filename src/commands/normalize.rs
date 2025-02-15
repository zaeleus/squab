use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use thiserror::Error;
use tracing::info;

use crate::{
    counts,
    normalization::{self, fpkm, tpm},
    read_features, Feature,
};

#[derive(Debug, Error)]
pub enum NormalizeError {
    #[error("I/O error")]
    Io(#[from] io::Error),
    #[error("could not open file: {1}")]
    OpenFile(#[source] io::Error, PathBuf),
    #[error("invalid counts")]
    ReadCounts(#[source] counts::ReadCountsError),
    #[error("invalid annotations")]
    ReadAnnotations(#[source] crate::ReadFeaturesError),
    #[error("normalization error")]
    Normalization(#[source] normalization::Error),
}

pub fn normalize<P, Q, R>(
    src: P,
    annotations_src: Q,
    feature_type: &str,
    id: &str,
    method: normalization::Method,
    dst: Option<R>,
) -> Result<(), NormalizeError>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: AsRef<Path>,
{
    let annotations_src = annotations_src.as_ref();

    info!(src = ?annotations_src, feature_type, feature_id = id, "reading features");

    let mut gff_reader = crate::gff::open(annotations_src)
        .map_err(|e| NormalizeError::OpenFile(e, annotations_src.into()))?;

    let (_, features) = read_features(&mut gff_reader, feature_type, id)
        .map_err(NormalizeError::ReadAnnotations)?;

    info!(feature_count = features.len(), "read features");

    let counts = read_counts(src)?;

    let names: Vec<_> = counts.iter().map(|(name, _)| name.clone()).collect();
    let counts: Vec<_> = counts.into_iter().map(|(_, count)| count).collect();

    info!(normalization_method = ?method, "normalizing counts");

    let lengths = calculate_feature_lengths(&features, &names)?;

    let normalized_counts = match method {
        normalization::Method::Fpkm => fpkm::normalize(&lengths, &counts),
        normalization::Method::Tpm => tpm::normalize(&lengths, &counts),
    };

    let mut writer: Box<dyn Write> = if let Some(dst) = dst {
        File::create(dst).map(BufWriter::new).map(Box::new)?
    } else {
        let stdout = io::stdout().lock();
        Box::new(BufWriter::new(stdout))
    };

    write_normalized_counts(&mut writer, &names, &normalized_counts)?;

    Ok(())
}

fn read_counts<P>(src: P) -> Result<Vec<(String, u32)>, NormalizeError>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map_err(NormalizeError::Io)?;

    counts::read(&mut reader).map_err(NormalizeError::ReadCounts)
}

fn calculate_feature_lengths(
    features: &HashMap<String, Vec<Feature>>,
    names: &[String],
) -> io::Result<Vec<u32>> {
    normalization::calculate_feature_lengths(features, names)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
        .into_iter()
        .map(|n| u32::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
        .collect()
}

fn write_normalized_counts<W>(writer: &mut W, names: &[String], values: &[f64]) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: char = '\t';

    for (name, value) in names.iter().zip(values) {
        writeln!(writer, "{name}{SEPARATOR}{value}")?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_normalized_counts() -> io::Result<()> {
        let names = [
            String::from("AADAT"),
            String::from("CLN3"),
            String::from("NEO1"),
            String::from("PAK4"),
        ];

        let values = [30.2, 3.7, 0.0, 14.5];

        let mut buf = Vec::new();
        write_normalized_counts(&mut buf, &names, &values)?;

        let expected = b"AADAT\t30.2\nCLN3\t3.7\nNEO1\t0\nPAK4\t14.5\n";
        assert_eq!(buf, expected);

        Ok(())
    }
}
