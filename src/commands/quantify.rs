use std::{
    fs::File,
    io::{self, BufWriter},
    path::Path,
    sync::Arc,
};

use log::{info, warn};
use noodles::Region;
use noodles_bam::{self as bam, bai};
use noodles_sam as sam;

use crate::{
    build_interval_trees,
    count::{
        count_paired_end_record_singletons, count_paired_end_records, count_single_end_records,
        Filter,
    },
    detect::{detect_specification, LibraryLayout},
    normalization::{self, calculate_fpkms, calculate_tpms},
    read_features,
    writer::{write_counts, write_normalized_count_values, write_stats},
    Context, Features, StrandSpecification, StrandSpecificationOption,
};

#[allow(clippy::too_many_arguments)]
pub fn quantify<P, Q, R>(
    bam_src: P,
    annotations_src: Q,
    feature_type: &str,
    id: &str,
    min_mapq: u8,
    with_secondary_records: bool,
    with_supplementary_records: bool,
    with_nonunique_records: bool,
    strand_specification_option: StrandSpecificationOption,
    threads: usize,
    normalize: Option<normalization::Method>,
    results_dst: R,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: AsRef<Path>,
{
    let feature_map = read_features(annotations_src, feature_type, id)?;
    let (features, names) = build_interval_trees(&feature_map);

    let mut reader = File::open(bam_src.as_ref()).map(bam::Reader::new)?;

    let header: sam::Header = reader
        .read_header()?
        .parse()
        // TODO: Pass `sam::header::ParseError` as error.
        .map_err(|_| io::Error::from(io::ErrorKind::InvalidData))?;

    let reference_sequences = header.reference_sequences();

    let mut feature_ids = Vec::with_capacity(names.len());
    feature_ids.extend(names.into_iter());
    feature_ids.sort();

    info!("detecting library type");

    let (library_layout, detected_strand_specification, strandedness_confidence) =
        detect_specification(&bam_src, reference_sequences, &features)?;

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

    let filter = Filter::new(
        min_mapq,
        with_secondary_records,
        with_supplementary_records,
        with_nonunique_records,
    );

    info!("counting features");

    info!("using {} thread(s)", threads);

    let mut runtime = tokio::runtime::Builder::new()
        .threaded_scheduler()
        .core_threads(threads)
        .build()?;

    let features = Arc::new(features);

    let ctx = match library_layout {
        LibraryLayout::SingleEnd => runtime.block_on(async {
            let tasks: Vec<_> = reference_sequences
                .values()
                .map(|reference_sequence| {
                    tokio::spawn(count_single_end_records_by_region(
                        bam_src.as_ref().to_path_buf(),
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

            Ok::<Context, io::Error>(ctx)
        })?,
        LibraryLayout::PairedEnd => runtime.block_on(async {
            let tasks: Vec<_> = reference_sequences
                .values()
                .map(|reference_sequence| {
                    tokio::spawn(count_paired_end_records_by_region(
                        bam_src.as_ref().to_path_buf(),
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
                reference_sequences,
                &filter,
                strand_specification,
            )?;

            let singletons = pairs.singletons().map(Ok);
            let ctx3 = count_paired_end_record_singletons(
                singletons,
                &features,
                reference_sequences,
                &filter,
                strand_specification,
            )?;

            ctx1.add(&ctx2);
            ctx1.add(&ctx3);

            Ok::<Context, io::Error>(ctx1)
        })?,
    };

    let mut writer = File::create(results_dst).map(BufWriter::new)?;

    if let Some(normalization_method) = normalize {
        match normalization_method {
            normalization::Method::Fpkm => {
                info!("calculating fpkms");
                let fpkms = calculate_fpkms(&ctx.counts, &feature_map)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                info!("writing fpkms");
                write_normalized_count_values(&mut writer, &fpkms, &feature_ids)?;
            }
            normalization::Method::Tpm => {
                info!("calculating tpms");
                let tpms = calculate_tpms(&ctx.counts, &feature_map)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                info!("writing tpms");
                write_normalized_count_values(&mut writer, &tpms, &feature_ids)?;
            }
        }
    } else {
        info!("writing counts");
        write_counts(&mut writer, &ctx.counts, &feature_ids)?;
        write_stats(&mut writer, &ctx)?;
    }

    Ok(())
}

async fn count_single_end_records_by_region<P>(
    bam_src: P,
    reference_sequence_name: String,
    features: Arc<Features>,
    filter: Filter,
    strand_specification: StrandSpecification,
) -> io::Result<Context>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(bam_src.as_ref()).map(bam::Reader::new)?;

    let header: sam::Header = reader
        .read_header()?
        .parse()
        // TODO: Pass `sam::header::ParseError` as error.
        .map_err(|_| io::Error::from(io::ErrorKind::InvalidData))?;

    let reference_sequences = header.reference_sequences();

    let bai_src = bam_src.as_ref().with_extension("bam.bai");
    let index = bai::read(bai_src)?;

    let region = Region::mapped(reference_sequence_name, 0, None);
    let query = reader.query(reference_sequences, &index, &region)?;

    count_single_end_records(
        query,
        &features,
        reference_sequences,
        &filter,
        strand_specification,
    )
}

async fn count_paired_end_records_by_region<P>(
    bam_src: P,
    reference_sequence_name: String,
    features: Arc<Features>,
    filter: Filter,
    strand_specification: StrandSpecification,
) -> io::Result<(Context, Vec<bam::Record>)>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(bam_src.as_ref()).map(bam::Reader::new)?;

    let header: sam::Header = reader
        .read_header()?
        .parse()
        // TODO: Pass `sam::header::ParseError` as error.
        .map_err(|_| io::Error::from(io::ErrorKind::InvalidData))?;

    let reference_sequences = header.reference_sequences();

    let bai_src = bam_src.as_ref().with_extension("bam.bai");
    let index = bai::read(bai_src)?;

    let region = Region::mapped(reference_sequence_name, 1, None);
    let query = reader.query(reference_sequences, &index, &region)?;

    let (ctx, mut pairs) = count_paired_end_records(
        query,
        &features,
        reference_sequences,
        &filter,
        strand_specification,
    )?;

    Ok((ctx, pairs.singletons().collect()))
}
