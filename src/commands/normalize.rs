use std::{fs::File, io, path::Path};

use log::info;

use crate::{
    normalization::{self, calculate_fpkms, calculate_tpms},
    read_features,
    reader::read_counts,
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
    let mut file = File::open(counts_src)?;
    let count_map = read_counts(&mut file)?;

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
