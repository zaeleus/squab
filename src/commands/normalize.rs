use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

use log::info;

use crate::{
    count,
    normalization::{self, calculate_fpkms, calculate_tpms},
    read_features,
    writer::write_normalized_count_values,
};

pub fn normalize<P, Q>(
    counts_src: P,
    annotations_src: Q,
    feature_type: &str,
    id: &str,
    method: normalization::Method,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let mut reader = File::open(counts_src).map(|f| count::Reader::new(BufReader::new(f)))?;
    let count_map = reader.read_counts()?;

    let feature_map = read_features(annotations_src, feature_type, id)?;
    let feature_ids: Vec<_> = feature_map.keys().map(|id| id.into()).collect();

    let stdout = io::stdout();
    let mut handle = stdout.lock();

    let values = match method {
        normalization::Method::Fpkm => {
            info!("calculating fpkms");
            calculate_fpkms(&count_map, &feature_map)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
        }
        normalization::Method::Tpm => {
            info!("calculating tpms");
            calculate_tpms(&count_map, &feature_map)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
        }
    };

    info!("writing values");
    write_normalized_count_values(&mut handle, &values, &feature_ids)?;

    Ok(())
}
