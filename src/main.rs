use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use clap::{crate_name, App, AppSettings, Arg, ArgMatches};
use git_testament::{git_testament, render_testament};
use noodles_squab::{commands, count::Filter};
use tracing::{info, warn};

git_testament!(TESTAMENT);

fn match_args_from_env() -> clap::ArgMatches {
    let quantify_cmd = App::new("quantify")
        .about("Gene expression quantification")
        .arg(
            Arg::new("with-secondary-records")
                .long("with-secondary-records")
                .help("Count secondary records (BAM flag 0x100)"),
        )
        .arg(
            Arg::new("with-supplementary-records")
                .long("with-supplementary-records")
                .help("Count supplementary records (BAM flag 0x800)"),
        )
        .arg(
            Arg::new("with-nonunique-records")
                .long("with-nonunique-records")
                .help("Count nonunique records (BAM data tag NH > 1)"),
        )
        .arg(
            Arg::new("strand-specification")
                .long("strand-specification")
                .value_name("str")
                .help("Strand specification")
                .possible_values(&["none", "forward", "reverse", "auto"])
                .default_value("auto"),
        )
        .arg(
            Arg::new("feature-type")
                .short('t')
                .long("feature-type")
                .value_name("str")
                .help("Feature type to count")
                .default_value("exon"),
        )
        .arg(
            Arg::new("id")
                .short('i')
                .long("id")
                .value_name("str")
                .help("Feature attribute to use as the feature identity")
                .default_value("gene_id"),
        )
        .arg(
            Arg::new("min-mapping-quality")
                .long("min-mapping-quality")
                .value_name("u8")
                .help("Minimum mapping quality to consider an alignment")
                .default_value("10"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("file")
                .help("Output destination for feature counts")
                .required(true),
        )
        .arg(
            Arg::new("annotations")
                .short('a')
                .long("annotations")
                .value_name("file")
                .help("Input annotations file (GFF3)")
                .required(true),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .value_name("uint")
                .help("Force a specific number of threads"),
        )
        .arg(
            Arg::new("bam")
                .help("Input alignment file")
                .required(true)
                .index(1),
        );

    let normalize_cmd = App::new("normalize")
        .about("Normalize counts")
        .arg(
            Arg::new("feature-type")
                .short('t')
                .long("feature-type")
                .value_name("str")
                .help("Feature type to count")
                .default_value("exon"),
        )
        .arg(
            Arg::new("id")
                .short('i')
                .long("id")
                .value_name("str")
                .help("Feature attribute to use as the feature identity")
                .default_value("gene_id"),
        )
        .arg(
            Arg::new("annotations")
                .short('a')
                .long("annotations")
                .value_name("file")
                .help("Input annotations file (GFF3)")
                .required(true),
        )
        .arg(
            Arg::new("method")
                .long("method")
                .value_name("str")
                .help("Quantification normalization method")
                .possible_values(&["fpkm", "tpm"])
                .default_value("tpm"),
        )
        .arg(
            Arg::new("counts")
                .help("Input counts file")
                .required(true)
                .index(1),
        );

    App::new(crate_name!())
        .version(render_testament!(TESTAMENT).as_str())
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .help("Use verbose logging")
                .hide(true),
        )
        .subcommand(quantify_cmd)
        .subcommand(normalize_cmd)
        .get_matches()
}

fn quantify(matches: &ArgMatches) -> anyhow::Result<()> {
    let bam_src = matches.value_of("bam").unwrap();
    let annotations_src = matches.value_of("annotations").unwrap();

    let results_dst = matches.value_of("output").unwrap();

    let feature_type = matches.value_of("feature-type").unwrap();
    let id = matches.value_of("id").unwrap();
    let min_mapping_quality = matches
        .value_of_t("min-mapping-quality")
        .unwrap_or_else(|e| e.exit());

    let with_secondary_records = matches.is_present("with-secondary-records");
    let with_supplementary_records = matches.is_present("with-supplementary-records");
    let with_nonunique_records = matches.is_present("with-nonunique-records");

    let threads = matches
        .value_of_t("threads")
        .unwrap_or_else(|_| num_cpus::get());

    let strand_specification_option = matches
        .value_of_t("strand-specification")
        .unwrap_or_else(|e| e.exit());

    let filter = Filter::new(
        min_mapping_quality,
        with_secondary_records,
        with_supplementary_records,
        with_nonunique_records,
    );

    info!("using {} thread(s)", threads);

    let runtime = tokio::runtime::Builder::new_multi_thread()
        .worker_threads(threads)
        .build()?;

    runtime.block_on(commands::quantify(
        bam_src,
        annotations_src,
        feature_type,
        id,
        filter,
        strand_specification_option,
        threads,
        results_dst,
    ))
}

fn normalize(matches: &ArgMatches) -> anyhow::Result<()> {
    let counts_src = matches.value_of("counts").unwrap();
    let annotations_src = matches.value_of("annotations").unwrap();

    let feature_type = matches.value_of("feature-type").unwrap();
    let id = matches.value_of("id").unwrap();

    let method = matches.value_of_t("method").unwrap_or_else(|e| e.exit());

    commands::normalize(counts_src, annotations_src, feature_type, id, method)
}

fn main() -> anyhow::Result<()> {
    let matches = match_args_from_env();

    tracing_subscriber::fmt::init();

    if matches.is_present("verbose") {
        warn!("`-v`/`--verbose` is deprecated and will be removed in a future version. Logging is now always enabled.");
    }

    if let Some(submatches) = matches.subcommand_matches("quantify") {
        quantify(submatches)
    } else if let Some(submatches) = matches.subcommand_matches("normalize") {
        normalize(submatches)
    } else {
        unreachable!()
    }
}
