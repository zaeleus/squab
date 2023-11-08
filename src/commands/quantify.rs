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

    let feature_map = read_features(&mut gff_reader, feature_type, id)?;
    let (features, names) = build_interval_trees(&feature_map);

    let mut reader = bgzf::reader::Builder::default()
        .set_worker_count(worker_count)
        .build_from_path(bam_src.as_ref())
        .map(bam::Reader::from)
        .with_context(|| format!("Could not open {}", bam_src.as_ref().display()))?;

    let header = reader.read_header()?;
    let reference_sequences = header.reference_sequences().clone();

    let mut feature_ids = Vec::with_capacity(names.len());
    feature_ids.extend(names);
    feature_ids.sort();

    info!("detecting library type");

    let (library_layout, detected_strand_specification, strandedness_confidence) =
        detect_specification(&bam_src, &reference_sequences, &features)?;

    info!("detected library layout: {library_layout}");

    match detected_strand_specification {
        StrandSpecification::None => info!(
            "strand specification: none (confidence: {:.2})",
            strandedness_confidence
        ),
        StrandSpecification::Forward => info!(
            "strand specification: forward (confidence: {:.2})",
            strandedness_confidence
        ),
        StrandSpecification::Reverse => info!(
            "strand specification: reverse (confidence: {:.2})",
            strandedness_confidence
        ),
    }

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
            &features,
            &reference_sequences,
            &filter,
            strand_specification,
            worker_count,
        )?,
        LibraryLayout::PairedEnd => count_paired_end_records(
            reader,
            &features,
            &reference_sequences,
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
