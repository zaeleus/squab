use std::{fs::File, io::BufWriter, num::NonZeroUsize, path::Path};

use anyhow::Context as AnyhowContext;
use noodles::{bam, bgzf};
use tracing::{info, warn};

use crate::{
    build_interval_trees,
    count::{self, count_paired_end_records, count_single_end_records, Filter},
    detect::{detect_specification, LibraryLayout},
    read_features, StrandSpecification, StrandSpecificationOption,
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
    results_dst: R,
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

    let (reference_sequence_names, feature_map) = read_features(&mut gff_reader, feature_type, id)?;
    let interval_trees = build_interval_trees(&header, &reference_sequence_names, &feature_map);

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

    let mut feature_ids: Vec<_> = feature_map.keys().collect();
    feature_ids.sort();

    info!("detecting library type");

    let (library_layout, detected_strand_specification, strandedness_confidence) =
        detect_specification(&bam_src, &interval_trees)?;

    info!("detected library layout: {library_layout}");
    info!("strand specification: {detected_strand_specification} (confidence: {strandedness_confidence:.2})");

    let strand_specification = match strand_specification_option {
        StrandSpecificationOption::None => StrandSpecification::None,
        StrandSpecificationOption::Forward => StrandSpecification::Forward,
        StrandSpecificationOption::Reverse => StrandSpecification::Reverse,
        StrandSpecificationOption::Auto => detected_strand_specification,
    };

    if strand_specification != detected_strand_specification {
        warn!(
            "input strand specification ({:?}) does not match detected strandedness ({:?})",
            strand_specification, detected_strand_specification,
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

    let writer = File::create(results_dst.as_ref())
        .map(BufWriter::new)
        .with_context(|| format!("Could not open {}", results_dst.as_ref().display()))?;

    info!("writing counts");

    let mut count_writer = count::Writer::new(writer);
    count_writer.write_counts(&feature_ids, &ctx.counts)?;
    count_writer.write_stats(&ctx)?;

    Ok(())
}
