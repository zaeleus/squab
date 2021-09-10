use std::{
    fs::File,
    io::{self, BufWriter},
    path::Path,
    sync::Arc,
};

use anyhow::Context as AnyhowContext;
use noodles::{bam, sam};
use tokio::io::AsyncRead;
use tracing::{info, warn};

use crate::{
    build_interval_trees,
    count::{self, count_paired_end_records, count_single_end_records, Filter},
    detect::{detect_specification, LibraryLayout},
    read_features, StrandSpecification, StrandSpecificationOption,
};

#[allow(clippy::too_many_arguments)]
pub async fn quantify<P, Q, R>(
    bam_src: P,
    annotations_src: Q,
    feature_type: &str,
    id: &str,
    filter: Filter,
    strand_specification_option: StrandSpecificationOption,
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

    let mut reader = tokio::fs::File::open(bam_src.as_ref())
        .await
        .map(bam::AsyncReader::new)
        .with_context(|| format!("Could not open {}", bam_src.as_ref().display()))?;

    let header = read_header(&mut reader).await?;
    let reference_sequences = header.reference_sequences().clone();

    let mut feature_ids = Vec::with_capacity(names.len());
    feature_ids.extend(names.into_iter());
    feature_ids.sort();

    info!("detecting library type");

    let (library_layout, detected_strand_specification, strandedness_confidence) =
        detect_specification(&bam_src, &reference_sequences, &features).await?;

    match library_layout {
        LibraryLayout::SingleEnd => info!("library layout: single end"),
        LibraryLayout::PairedEnd => info!("library layout: paired end"),
    }

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

    let reference_sequences = Arc::new(reference_sequences);
    let features = Arc::new(features);

    let ctx = match library_layout {
        LibraryLayout::SingleEnd => {
            count_single_end_records(
                reader,
                features.clone(),
                reference_sequences.clone(),
                filter.clone(),
                strand_specification,
            )
            .await?
        }
        LibraryLayout::PairedEnd => {
            count_paired_end_records(
                reader,
                features.clone(),
                reference_sequences.clone(),
                filter.clone(),
                strand_specification,
            )
            .await?
        }
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

async fn read_header<R>(reader: &mut bam::AsyncReader<R>) -> anyhow::Result<sam::Header>
where
    R: AsyncRead + Unpin,
{
    let mut header: sam::Header = reader
        .read_header()
        .await?
        .parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .context("Could not parse BAM header")?;

    let reference_sequences = reader.read_reference_sequences().await?;

    if header.reference_sequences().is_empty() {
        *header.reference_sequences_mut() = reference_sequences;
    }

    Ok(header)
}
