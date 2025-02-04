use std::{
    fs::File,
    io::{self, BufReader},
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
    let src = src.as_ref();
    let annotations_src = annotations_src.as_ref();

    let counts = read_counts(src)?.into_iter().collect();

    let mut gff_reader = crate::gff::open(annotations_src)
        .map_err(|e| NormalizeError::OpenFile(e, annotations_src.into()))?;

    info!(src = ?annotations_src, feature_type, feature_id = id, "reading features");

    let (_, features) = read_features(&mut gff_reader, feature_type, id)
        .map_err(NormalizeError::ReadAnnotations)?;

    info!(feature_count = features.len(), "read features");
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

    let stdout = io::stdout().lock();
    let mut writer = normalization::Writer::new(stdout);

    let mut feature_names: Vec<_> = features.keys().collect();
    feature_names.sort();

    writer
        .write_values(&feature_names, &values)
        .map_err(NormalizeError::Io)?;

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
