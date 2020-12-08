use std::{
    fs::File,
    io::{self, BufWriter},
    path::Path,
    sync::Arc,
};

use anyhow::Context as AnyhowContext;
use log::{info, warn};
use noodles::Region;
use noodles_bam::{self as bam, bai};
use noodles_sam::{self as sam, header::ReferenceSequences};

use crate::{
    build_interval_trees,
    count::{
        self, count_paired_end_record_singletons, count_paired_end_records,
        count_single_end_records, Filter,
    },
    detect::{detect_specification, LibraryLayout},
    normalization::{self, calculate_fpkms, calculate_tpms},
    read_features, Context, Features, StrandSpecification, StrandSpecificationOption,
};

#[allow(clippy::too_many_arguments)]
pub fn quantify<P, Q, R>(
    bam_src: P,
    annotations_src: Q,
    feature_type: &str,
    id: &str,
    filter: Filter,
    strand_specification_option: StrandSpecificationOption,
    threads: usize,
    normalize: Option<normalization::Method>,
    results_dst: R,
) -> anyhow::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: AsRef<Path>,
{
    let mut gff_reader = crate::gff::open(annotations_src)?;
    let feature_map = read_features(&mut gff_reader, feature_type, id)?;
    let (features, names) = build_interval_trees(&feature_map);

    let mut reader = File::open(bam_src.as_ref())
        .map(bam::Reader::new)
        .with_context(|| format!("Could not open {}", bam_src.as_ref().display()))?;

    let header: sam::Header = reader
        .read_header()?
        .parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .context("Could not parse BAM header")?;

    let bai_src = bam_src.as_ref().with_extension("bam.bai");
    let index =
        bai::read(&bai_src).with_context(|| format!("Could not read {}", bai_src.display()))?;

    let reference_sequences = header.reference_sequences().clone();

    let mut feature_ids = Vec::with_capacity(names.len());
    feature_ids.extend(names.into_iter());
    feature_ids.sort();

    info!("detecting library type");

    let (library_layout, detected_strand_specification, strandedness_confidence) =
        detect_specification(&bam_src, &reference_sequences, &features)?;

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

    info!("using {} thread(s)", threads);

    let runtime = tokio::runtime::Builder::new_multi_thread()
        .worker_threads(threads)
        .build()?;

    let index = Arc::new(index);
    let reference_sequences = Arc::new(reference_sequences);
    let features = Arc::new(features);

    let ctx = runtime.block_on(async {
        match library_layout {
            LibraryLayout::SingleEnd => {
                let tasks: Vec<_> = reference_sequences
                    .values()
                    .map(|reference_sequence| {
                        tokio::spawn(count_single_end_records_by_region(
                            bam_src.as_ref().to_path_buf(),
                            index.clone(),
                            reference_sequences.clone(),
                            reference_sequence.name().into(),
                            features.clone(),
                            filter.clone(),
                            strand_specification,
                        ))
                    })
                    .collect();

                let mut ctx = Context::default();

                for task in tasks {
                    let region_ctx = task.await??;
                    ctx.add(&region_ctx);
                }

                Ok::<Context, anyhow::Error>(ctx)
            }
            LibraryLayout::PairedEnd => {
                let tasks: Vec<_> = reference_sequences
                    .values()
                    .map(|reference_sequence| {
                        tokio::spawn(count_paired_end_records_by_region(
                            bam_src.as_ref().to_path_buf(),
                            index.clone(),
                            reference_sequences.clone(),
                            reference_sequence.name().into(),
                            features.clone(),
                            filter.clone(),
                            strand_specification,
                        ))
                    })
                    .collect();

                let mut ctx1 = Context::default();
                let mut pairs = Vec::with_capacity(reference_sequences.len());

                for task in tasks {
                    let (region_ctx, region_pairs) = task.await??;
                    ctx1.add(&region_ctx);
                    pairs.push(region_pairs);
                }

                let records = pairs.into_iter().flat_map(|r| r.into_iter()).map(Ok);
                let (ctx2, mut pairs) = count_paired_end_records(
                    records,
                    &features,
                    &reference_sequences,
                    &filter,
                    strand_specification,
                )?;

                let singletons = pairs.singletons().map(Ok);
                let ctx3 = count_paired_end_record_singletons(
                    singletons,
                    &features,
                    &reference_sequences,
                    &filter,
                    strand_specification,
                )?;

                ctx1.add(&ctx2);
                ctx1.add(&ctx3);

                Ok::<Context, anyhow::Error>(ctx1)
            }
        }
    })?;

    let writer = File::create(results_dst.as_ref())
        .map(BufWriter::new)
        .with_context(|| format!("Could not open {}", results_dst.as_ref().display()))?;

    if let Some(normalization_method) = normalize {
        let mut value_writer = normalization::Writer::new(writer);

        match normalization_method {
            normalization::Method::Fpkm => {
                info!("calculating fpkms");
                let fpkms = calculate_fpkms(&ctx.counts, &feature_map)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                info!("writing fpkms");
                value_writer.write_values(&feature_ids, &fpkms)?;
            }
            normalization::Method::Tpm => {
                info!("calculating tpms");
                let tpms = calculate_tpms(&ctx.counts, &feature_map)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                info!("writing tpms");
                value_writer.write_values(&feature_ids, &tpms)?;
            }
        }
    } else {
        info!("writing counts");

        let mut count_writer = count::Writer::new(writer);
        count_writer.write_counts(&feature_ids, &ctx.counts)?;
        count_writer.write_stats(&ctx)?;
    }

    Ok(())
}

async fn count_single_end_records_by_region<P>(
    bam_src: P,
    index: Arc<bai::Index>,
    reference_sequences: Arc<ReferenceSequences>,
    reference_sequence_name: String,
    features: Arc<Features>,
    filter: Filter,
    strand_specification: StrandSpecification,
) -> anyhow::Result<Context>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(bam_src.as_ref())
        .map(bam::Reader::new)
        .with_context(|| format!("Could not open {}", bam_src.as_ref().display()))?;

    let reference_sequence_len = reference_sequences
        .get(&reference_sequence_name)
        .map(|rs| rs.len())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid reference sequence name",
            )
        })?;

    let region = Region::mapped(reference_sequence_name, 1, reference_sequence_len);
    let query = reader.query(&reference_sequences, &index, &region)?;

    let ctx = count_single_end_records(
        query,
        &features,
        &reference_sequences,
        &filter,
        strand_specification,
    )?;

    Ok(ctx)
}

async fn count_paired_end_records_by_region<P>(
    bam_src: P,
    index: Arc<bai::Index>,
    reference_sequences: Arc<ReferenceSequences>,
    reference_sequence_name: String,
    features: Arc<Features>,
    filter: Filter,
    strand_specification: StrandSpecification,
) -> anyhow::Result<(Context, Vec<bam::Record>)>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(bam_src.as_ref())
        .map(bam::Reader::new)
        .with_context(|| format!("Could not open {}", bam_src.as_ref().display()))?;

    let reference_sequence_len = reference_sequences
        .get(&reference_sequence_name)
        .map(|rs| rs.len())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid reference sequence name",
            )
        })?;

    let region = Region::mapped(reference_sequence_name, 1, reference_sequence_len);
    let query = reader.query(&reference_sequences, &index, &region)?;

    let (ctx, mut pairs) = count_paired_end_records(
        query,
        &features,
        &reference_sequences,
        &filter,
        strand_specification,
    )?;

    Ok((ctx, pairs.singletons().collect()))
}
