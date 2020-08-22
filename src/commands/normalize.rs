use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

use anyhow::Context;
use log::info;

use crate::{
    count,
    normalization::{self, calculate_fpkms, calculate_tpms},
    read_features,
};

pub fn normalize<P, Q>(
    counts_src: P,
    annotations_src: Q,
    feature_type: &str,
    id: &str,
    method: normalization::Method,
) -> anyhow::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let mut reader = File::open(counts_src.as_ref())
        .map(BufReader::new)
        .map(count::Reader::new)
        .with_context(|| format!("Could not open {}", counts_src.as_ref().display()))?;

    let count_map = reader
        .read_counts()
        .with_context(|| format!("Could not read {}", counts_src.as_ref().display()))?;

    let mut gff_reader = crate::gff::open(annotations_src.as_ref())
        .with_context(|| format!("Could not open {}", annotations_src.as_ref().display()))?;

    let feature_map = read_features(&mut gff_reader, feature_type, id)
        .with_context(|| format!("Could not read {}", annotations_src.as_ref().display()))?;

    let feature_ids: Vec<_> = feature_map.keys().map(|id| id.into()).collect();

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = normalization::Writer::new(handle);

    let values = match method {
        normalization::Method::Fpkm => {
            info!("calculating fpkms");

            calculate_fpkms(&count_map, &feature_map)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                .context("Could not calculate FPKM")?
        }
        normalization::Method::Tpm => {
            info!("calculating tpms");

            calculate_tpms(&count_map, &feature_map)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                .context("Could not calculate TPM")?
        }
    };

    writer
        .write_values(&feature_ids, &values)
        .context("Could not write to stdout")?;

    Ok(())
}
