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
    normalization::{self, calculate_feature_lengths, calculate_fpkms, calculate_tpms},
    read_features,
};

#[derive(Debug, Error)]
pub enum NormalizeError {
    #[error("I/O error")]
    Io(#[source] io::Error),
    #[error("could not open file: {1}")]
    OpenFile(#[source] io::Error, PathBuf),
    #[error("invalid counts")]
    ReadCounts(#[source] counts::ReadCountsError),
    #[error("invalid annotations")]
    ReadAnnotations(#[source] crate::ReadFeaturesError),
    #[error("normalization error")]
    Normalization(#[source] normalization::Error),
}

pub fn normalize<P, Q>(
    src: P,
    annotations_src: Q,
    feature_type: &str,
    id: &str,
    method: normalization::Method,
) -> Result<(), NormalizeError>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let annotations_src = annotations_src.as_ref();

    info!(src = ?annotations_src, feature_type, feature_id = id, "reading features");

    let mut gff_reader = crate::gff::open(annotations_src)
        .map_err(|e| NormalizeError::OpenFile(e, annotations_src.into()))?;

    let (_, features) = read_features(&mut gff_reader, feature_type, id)
        .map_err(NormalizeError::ReadAnnotations)?;

    info!(feature_count = features.len(), "read features");

    let counts = read_counts(src)?.into_iter().collect();

    info!(normalization_method = ?method, "normalizing counts");

    let lengths = calculate_feature_lengths(&features);

    let values = match method {
        normalization::Method::Fpkm => {
            calculate_fpkms(&lengths, &counts).map_err(NormalizeError::Normalization)?
        }
        normalization::Method::Tpm => {
            calculate_tpms(&lengths, &counts).map_err(NormalizeError::Normalization)?
        }
    };

    let mut feature_names: Vec<_> = features.keys().collect();
    feature_names.sort();

    let stdout = io::stdout().lock();
    let mut writer = BufWriter::new(stdout);

    write_normalized_counts(&mut writer, &feature_names, &values).map_err(NormalizeError::Io)?;

    Ok(())
}

fn read_counts<P>(src: P) -> Result<Vec<(String, u64)>, NormalizeError>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map_err(NormalizeError::Io)?;

    counts::read(&mut reader).map_err(NormalizeError::ReadCounts)
}

fn write_normalized_counts<W>(
    writer: &mut W,
    names: &[&String],
    values: &HashMap<&str, f64>,
) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: char = '\t';

    for name in names {
        let value = values.get(name.as_str()).unwrap_or(&0.0);
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
            &String::from("AADAT"),
            &String::from("CLN3"),
            &String::from("NEO1"),
            &String::from("PAK4"),
        ];

        let values = [("AADAT", 30.2), ("CLN3", 3.7), ("PAK4", 14.5)]
            .into_iter()
            .collect();

        let mut buf = Vec::new();
        write_normalized_counts(&mut buf, &names, &values)?;

        let expected = b"AADAT\t30.2
CLN3\t3.7
NEO1\t0
PAK4\t14.5
";

        assert_eq!(buf, expected);

        Ok(())
    }
}
