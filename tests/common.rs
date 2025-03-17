use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

use noodles::{
    bam,
    sam::{self, alignment::io::Write},
};

pub fn convert_sam_to_bam<P, Q>(src: P, dst: Q) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;

    let header = reader.read_header()?;

    let mut writer = File::create(dst).map(bam::io::Writer::new)?;
    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    Ok(())
}
