use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::{path::PathBuf, thread};

use clap::{crate_name, value_parser, Arg, ArgMatches, Command};
use git_testament::{git_testament, render_testament};
use squab::{commands, count::Filter, normalization::Method, StrandSpecificationOption};
use tracing::{info, warn};

git_testament!(TESTAMENT);

fn match_args_from_env() -> clap::ArgMatches {
    let quantify_cmd = Command::new("quantify")
        .about("Gene expression quantification")
        .arg(
            Arg::new("with-secondary-records")
                .long("with-secondary-records")
                .help("Count secondary records (BAM flag 0x100)")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("with-supplementary-records")
                .long("with-supplementary-records")
                .help("Count supplementary records (BAM flag 0x800)")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("with-nonunique-records")
                .long("with-nonunique-records")
                .help("Count nonunique records (BAM data tag NH > 1)")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("strand-specification")
                .long("strand-specification")
                .value_name("str")
                .help("Strand specification")
                .default_value("auto")
                .value_parser(value_parser!(StrandSpecificationOption)),
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
                .default_value("10")
                .value_parser(value_parser!(u8)),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("path")
                .help("Output destination for feature counts")
                .value_parser(value_parser!(PathBuf))
                .required(true),
        )
        .arg(
            Arg::new("annotations")
                .short('a')
                .long("annotations")
                .value_name("path")
                .help("Input annotations file (GFF3)")
                .value_parser(value_parser!(PathBuf))
                .required(true),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .value_name("uint")
                .help("Force a specific number of threads")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("bam")
                .index(1)
                .help("Input alignment file")
                .value_parser(value_parser!(PathBuf))
                .required(true),
        );

    let normalize_cmd = Command::new("normalize")
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
                .value_parser(value_parser!(PathBuf))
                .required(true),
        )
        .arg(
            Arg::new("method")
                .long("method")
                .value_name("str")
                .help("Quantification normalization method")
                .default_value("tpm")
                .value_parser(value_parser!(Method)),
        )
        .arg(
            Arg::new("counts")
                .index(1)
                .help("Input counts file")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        );

    Command::new(crate_name!())
        .version(render_testament!(TESTAMENT))
        .subcommand_required(true)
        .arg_required_else_help(true)
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .help("Use verbose logging")
                .action(clap::ArgAction::SetTrue)
                .hide(true),
        )
        .subcommand(quantify_cmd)
        .subcommand(normalize_cmd)
        .get_matches()
}

fn quantify(matches: &ArgMatches) -> anyhow::Result<()> {
    let bam_src: &PathBuf = matches.get_one("bam").unwrap();
    let annotations_src: &PathBuf = matches.get_one("annotations").unwrap();

    let results_dst: &PathBuf = matches.get_one("output").unwrap();

    let feature_type: &String = matches.get_one("feature-type").unwrap();
    let id: &String = matches.get_one("id").unwrap();
    let min_mapping_quality = *matches.get_one("min-mapping-quality").unwrap();

    let with_secondary_records = *matches.get_one("with-secondary-records").unwrap();
    let with_supplementary_records = *matches.get_one("with-supplementary-records").unwrap();
    let with_nonunique_records = *matches.get_one("with-nonunique-records").unwrap();

    let threads = match matches.get_one("threads") {
        Some(n) => *n,
        None => thread::available_parallelism().map(usize::from)?,
    };

    let strand_specification_option = *matches.get_one("strand-specification").unwrap();

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
    let counts_src: &PathBuf = matches.get_one("counts").unwrap();
    let annotations_src: &PathBuf = matches.get_one("annotations").unwrap();

    let feature_type: &String = matches.get_one("feature-type").unwrap();
    let id: &String = matches.get_one("id").unwrap();

    let method = *matches.get_one("method").unwrap();

    commands::normalize(counts_src, annotations_src, feature_type, id, method)
}

fn main() -> anyhow::Result<()> {
    let matches = match_args_from_env();

    tracing_subscriber::fmt::init();

    if let Some(true) = matches.get_one("verbose").copied() {
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
