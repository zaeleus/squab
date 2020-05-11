use std::io;

use clap::{crate_name, value_t, App, AppSettings, Arg, ArgMatches, SubCommand};
use git_testament::{git_testament, render_testament};
use log::LevelFilter;
use noodles_squab::{commands, normalization, StrandSpecificationOption};

git_testament!(TESTAMENT);

fn match_args_from_env() -> clap::ArgMatches<'static> {
    let quantify_cmd = SubCommand::with_name("quantify")
        .about("Gene expression quantification")
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
        );

    let normalize_cmd = SubCommand::with_name("normalize")
        .about("Normalize counts")
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
            Arg::with_name("annotations")
                .short("a")
                .long("annotations")
                .value_name("file")
                .help("Input annotations file (GTF/GFFv2)")
                .required(true),
        )
        .arg(
            Arg::with_name("method")
                .long("method")
                .value_name("str")
                .help("Quantification normalization method")
                .possible_values(&["fpkm", "tpm"])
                .default_value("tpm"),
        )
        .arg(
            Arg::with_name("counts")
                .help("Input counts file")
                .required(true)
                .index(1),
        );

    App::new(crate_name!())
        .version(render_testament!(TESTAMENT).as_str())
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .long("verbose")
                .help("Use verbose logging"),
        )
        .subcommand(quantify_cmd)
        .subcommand(normalize_cmd)
        .get_matches()
}

fn quantify(matches: &ArgMatches<'_>) -> io::Result<()> {
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

    commands::quantify(
        bam_src,
        annotations_src,
        feature_type,
        id,
        min_mapq,
        with_secondary_records,
        with_supplementary_records,
        with_nonunique_records,
        strand_specification_option,
        threads,
        normalize,
        results_dst,
    )
}

fn normalize(matches: &ArgMatches<'_>) -> io::Result<()> {
    let counts_src = matches.value_of("counts").unwrap();
    let annotations_src = matches.value_of("annotations").unwrap();

    let feature_type = matches.value_of("type").unwrap();
    let id = matches.value_of("id").unwrap();

    let method = value_t!(matches, "method", normalization::Method).unwrap_or_else(|e| e.exit());

    commands::normalize(counts_src, annotations_src, feature_type, id, method)
}

fn main() {
    let matches = match_args_from_env();

    if matches.is_present("verbose") {
        env_logger::Builder::from_default_env()
            .filter(Some("noodles_squab"), LevelFilter::Info)
            .init();
    } else {
        env_logger::init();
    }

    if let Some(submatches) = matches.subcommand_matches("quantify") {
        quantify(submatches).unwrap();
    } else if let Some(submatches) = matches.subcommand_matches("normalize") {
        normalize(submatches).unwrap();
    }
}
