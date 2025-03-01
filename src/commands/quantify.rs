use std::{
    fs::File,
    io::{self, BufWriter, Write},
    num::NonZeroUsize,
    path::Path,
};

use anyhow::Context as AnyhowContext;
use noodles::{bam, bgzf};
use tracing::{info, warn};

use crate::{
    StrandSpecification, StrandSpecificationOption, build_interval_trees,
    count::{self, Filter, count_paired_end_records, count_single_end_records},
    detect::{LibraryLayout, detect_specification},
    read_features,
};

#[allow(clippy::too_many_arguments)]
pub fn quantify<P, Q, R>(
    bam_src: P,
    annotations_src: Q,
    feature_type: &str,
    id: &str,
    filter: Filter,
    strand_specification_option: StrandSpecificationOption,
    worker_count: NonZeroUsize,
    dst: Option<R>,
) -> anyhow::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: AsRef<Path>,
{
    let mut gff_reader = crate::gff::open(annotations_src.as_ref())
        .with_context(|| format!("Could not open {}", annotations_src.as_ref().display()))?;

    let mut reader = bam::io::reader::Builder.build_from_path(bam_src.as_ref())?;
    let header = reader.read_header()?;

    info!(src = ?annotations_src.as_ref(), feature_type, feature_id = id, "reading features");

    let (reference_sequence_names, features) = read_features(&mut gff_reader, feature_type, id)?;

    info!(feature_count = features.len(), "read features");

    let interval_trees = build_interval_trees(&header, &reference_sequence_names, &features);

    let decoder: Box<dyn bgzf::io::Read + Send> = if worker_count.get() > 1 {
        File::open(bam_src.as_ref())
            .map(|f| bgzf::MultithreadedReader::with_worker_count(worker_count, f))
            .map(Box::new)
            .with_context(|| format!("Could not open {}", bam_src.as_ref().display()))?
    } else {
        bgzf::reader::Builder
            .build_from_path(bam_src.as_ref())
            .map(Box::new)
            .with_context(|| format!("Could not open {}", bam_src.as_ref().display()))?
    };

    let mut reader = bam::io::Reader::from(decoder);
    reader.read_header()?;

    info!("detecting library type");

    let (library_layout, detected_strand_specification, strandedness_confidence) =
        detect_specification(&bam_src, &interval_trees)?;

    info!("detected library layout: {library_layout}");
    info!(
        "strand specification: {detected_strand_specification} (confidence: {strandedness_confidence:.2})"
    );

    let strand_specification = match strand_specification_option {
        StrandSpecificationOption::None => StrandSpecification::None,
        StrandSpecificationOption::Forward => StrandSpecification::Forward,
        StrandSpecificationOption::Reverse => StrandSpecification::Reverse,
        StrandSpecificationOption::Auto => detected_strand_specification,
    };

    if strand_specification != detected_strand_specification {
        warn!(
            expected = ?strand_specification,
            actual = ?detected_strand_specification,
            "strand specification mismatch",
        );
    }

    info!("counting features");

    let ctx = match library_layout {
        LibraryLayout::SingleEnd => count_single_end_records(
            reader,
            &interval_trees,
            &filter,
            strand_specification,
            worker_count,
        )?,
        LibraryLayout::PairedEnd => count_paired_end_records(
            reader,
            &interval_trees,
            &filter,
            strand_specification,
            worker_count,
        )?,
    };

    let writer: Box<dyn Write> = if let Some(dst) = dst {
        File::create(dst.as_ref())
            .map(BufWriter::new)
            .map(Box::new)
            .with_context(|| format!("Could not open {}", dst.as_ref().display()))?
    } else {
        let stdout = io::stdout().lock();
        Box::new(BufWriter::new(stdout))
    };

    info!("writing counts");

    let mut count_writer = count::Writer::new(writer);

    let mut feature_names: Vec<_> = features.keys().collect();
    feature_names.sort();

    count_writer.write_counts(&feature_names, &ctx.counts)?;
    count_writer.write_stats(&ctx)?;

    Ok(())
}
