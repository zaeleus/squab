use log::info;
use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

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
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let mut reader = File::open(counts_src).map(|f| count::Reader::new(BufReader::new(f)))?;
    let count_map = reader.read_counts()?;

    let mut gff_reader = crate::gff::open(annotations_src)?;
    let feature_map = read_features(&mut gff_reader, feature_type, id)?;
    let feature_ids: Vec<_> = feature_map.keys().map(|id| id.into()).collect();

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = normalization::Writer::new(handle);

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

    writer.write_values(&feature_ids, &values)?;

    Ok(())
}
