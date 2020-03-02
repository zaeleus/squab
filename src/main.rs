use std::{convert::TryFrom, fs::File, io::BufWriter, path::Path, sync::Arc};

use clap::{crate_name, value_t, App, Arg};
use git_testament::{git_testament, render_testament};
use log::{info, warn, LevelFilter};
use noodles::formats::bai;
use noodles_bam as bam;
use noodles_squab::{
    build_interval_trees,
    count::{count_paired_end_record_singletons, count_paired_end_records, Filter},
    count_single_end_records,
    detect::{detect_specification, LibraryLayout},
    normalization::{self, calculate_fpkms, calculate_tpms},
    read_features,
    writer::{write_counts, write_normalized_count_values, write_stats},
    Context, Features, StrandSpecification, StrandSpecificationOption,
};

git_testament!(TESTAMENT);

async fn count_single_end_records_by_region<P>(
    bam_src: P,
    ref_id: usize,
    reference: bam::Reference,
    features: Arc<Features>,
    references: Arc<Vec<bam::Reference>>,
    filter: Filter,
    strand_specification: StrandSpecification,
) -> Context
where
    P: AsRef<Path>,
{
    let file = File::open(bam_src.as_ref()).unwrap();
    let mut reader = bam::Reader::new(file);

    let bai_src = bam_src.as_ref().with_extension("bam.bai");
    let bai_file = File::open(bai_src).unwrap();
    let mut bai_reader = bai::Reader::new(bai_file);
    bai_reader.header().unwrap();
    let index = bai_reader.read_index().unwrap();

    let index_ref = &index.references[ref_id];

    let start = 0;
    let end = reference.len();

    let query = reader.query(index_ref, start, end);

    count_single_end_records(query, &features, &references, &filter, strand_specification).unwrap()
}

fn match_args_from_env() -> clap::ArgMatches<'static> {
    App::new(crate_name!())
        .version(render_testament!(TESTAMENT).as_str())
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .long("verbose")
                .help("Use verbose logging"),
        )
        .arg(
            Arg::with_name("with-secondary-records")
                .long("with-secondary-records")
                .help("Count secondary records (BAM flag 0x100)"),
        )
        .arg(
            Arg::with_name("with-supplementary-records")
                .long("with-supplementary-records")
                .help("Count supplementary records (BAM flag 0x800)"),
        )
        .arg(
            Arg::with_name("with-nonunique-records")
                .long("with-nonunique-records")
                .help("Count nonunique records (BAM data tag NH > 1)"),
        )
        .arg(
            Arg::with_name("strand-specification")
                .long("strand-specification")
                .value_name("str")
                .help("Strand specification")
                .possible_values(&["none", "forward", "reverse", "auto"])
                .default_value("auto"),
        )
        .arg(
            Arg::with_name("type")
                .short("t")
                .long("type")
                .value_name("str")
                .help("Feature type to count")
                .default_value("exon"),
        )
        .arg(
            Arg::with_name("id")
                .short("i")
                .long("id")
                .value_name("str")
                .help("Feature attribute to use as the feature identity")
                .default_value("gene_id"),
        )
        .arg(
            Arg::with_name("min-mapq")
                .long("min-mapq")
                .value_name("u8")
                .help("Minimum mapping quality to consider an alignment")
                .default_value("10"),
        )
        .arg(
            Arg::with_name("normalize")
                .long("normalize")
                .value_name("str")
                .help("Quantification normalization method")
                .possible_values(&["fpkm", "tpm"]),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .value_name("file")
                .help("Output destination for feature counts")
                .required(true),
        )
        .arg(
            Arg::with_name("annotations")
                .short("a")
                .long("annotations")
                .value_name("file")
                .help("Input annotations file (GTF/GFFv2)")
                .required(true),
        )
        .arg(
            Arg::with_name("threads")
                .long("threads")
                .value_name("uint")
                .help("Force a specific number of threads"),
        )
        .arg(
            Arg::with_name("bam")
                .help("Input alignment file")
                .required(true)
                .index(1),
        )
        .get_matches()
}

#[allow(clippy::cognitive_complexity)]
fn main() {
    let matches = match_args_from_env();

    if matches.is_present("verbose") {
        env_logger::Builder::from_default_env()
            .filter(Some("noodles_squab"), LevelFilter::Info)
            .init();
    } else {
        env_logger::init();
    }

    let bam_src = matches.value_of("bam").unwrap();
    let annotations_src = matches.value_of("annotations").unwrap();

    let normalize = matches.value_of("normalize").map(|_| {
        value_t!(matches, "normalize", normalization::Method).unwrap_or_else(|e| e.exit())
    });

    let results_dst = matches.value_of("output").unwrap();

    let feature_type = matches.value_of("type").unwrap();
    let id = matches.value_of("id").unwrap();
    let min_mapq = value_t!(matches, "min-mapq", u8).unwrap_or_else(|e| e.exit());

    let with_secondary_records = matches.is_present("with-secondary-records");
    let with_supplementary_records = matches.is_present("with-supplementary-records");
    let with_nonunique_records = matches.is_present("with-nonunique-records");

    let threads = value_t!(matches, "threads", usize).unwrap_or_else(|_| num_cpus::get());

    let strand_specification_option =
        value_t!(matches, "strand-specification", StrandSpecificationOption)
            .unwrap_or_else(|e| e.exit());

    let feature_map = read_features(annotations_src, feature_type, id).unwrap();
    let (features, names) = build_interval_trees(&feature_map);

    let file = File::open(&bam_src).unwrap();
    let mut reader = bam::Reader::new(file);
    let (_, references) = reader.header().expect("failed to read bam header");

    let mut feature_ids = Vec::with_capacity(names.len());
    feature_ids.extend(names.into_iter());
    feature_ids.sort();

    info!("detecting library type");

    let (library_layout, detected_strand_specification, strandedness_confidence) =
        detect_specification(&bam_src, &references, &features).unwrap();

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
        StrandSpecificationOption::Auto => detected_strand_specification,
        _ => StrandSpecification::try_from(strand_specification_option).unwrap(),
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

    let ctx = match library_layout {
        LibraryLayout::SingleEnd => {
            info!("using {} thread(s)", threads);

            let mut runtime = tokio::runtime::Builder::new()
                .threaded_scheduler()
                .core_threads(threads)
                .build()
                .unwrap();

            let features = Arc::new(features);
            let references = Arc::new(references);

            runtime.block_on(async {
                let tasks: Vec<_> = references
                    .iter()
                    .enumerate()
                    .map(|(ref_id, reference)| {
                        tokio::spawn(count_single_end_records_by_region(
                            bam_src.to_string(),
                            ref_id,
                            reference.clone(),
                            features.clone(),
                            references.clone(),
                            filter.clone(),
                            strand_specification,
                        ))
                    })
                    .collect();

                let mut ctx = Context::default();

                for task in tasks {
                    let region_ctx = task.await.unwrap();
                    ctx.add(&region_ctx);
                }

                ctx
            })
        }
        LibraryLayout::PairedEnd => {
            let records = reader.records();
            let (mut ctx1, mut pairs) = count_paired_end_records(
                records,
                &features,
                &references,
                &filter,
                strand_specification,
            )
            .unwrap();

            let singletons = pairs.singletons().map(Ok);
            let ctx2 = count_paired_end_record_singletons(
                singletons,
                &features,
                &references,
                &filter,
                strand_specification,
            )
            .unwrap();

            ctx1.add(&ctx2);

            ctx1
        }
    };

    let file = File::create(results_dst).unwrap();
    let mut writer = BufWriter::new(file);

    if let Some(normalization_method) = normalize {
        match normalization_method {
            normalization::Method::Fpkm => {
                info!("calculating fpkms");
                let fpkms = calculate_fpkms(&ctx.counts, &feature_map).unwrap();
                info!("writing fpkms");
                write_normalized_count_values(&mut writer, &fpkms, &feature_ids).unwrap();
            }
            normalization::Method::Tpm => {
                info!("calculating tpms");
                let tpms = calculate_tpms(&ctx.counts, &feature_map).unwrap();
                info!("writing tpms");
                write_normalized_count_values(&mut writer, &tpms, &feature_ids).unwrap();
            }
        }
    } else {
        info!("writing counts");
        write_counts(&mut writer, &ctx.counts, &feature_ids).unwrap();
        write_stats(&mut writer, &ctx).unwrap();
    }
}
