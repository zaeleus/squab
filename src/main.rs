use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::thread;

use clap::Parser;
use squab::{
    cli::{self, Command},
    commands,
    count::Filter,
    Cli,
};
use tracing::info;

fn quantify(options: cli::Quantify) -> anyhow::Result<()> {
    let bam_src = options.src;
    let annotations_src = options.annotations;

    let results_dst = options.output;

    let threads = match options.threads {
        Some(n) => n,
        None => thread::available_parallelism()?,
    };

    let strand_specification_option = options.strand_specification;

    let filter = Filter::new(
        options.min_mapping_quality,
        options.with_secondary_records,
        options.with_supplementary_records,
        options.with_nonunique_records,
    );

    info!("using {} thread(s)", threads);

    commands::quantify(
        bam_src,
        annotations_src,
        &options.feature_type,
        &options.id,
        filter,
        strand_specification_option,
        threads,
        results_dst,
    )
}

fn normalize(options: cli::Normalize) -> anyhow::Result<()> {
    let counts_src = options.src;
    let annotations_src = options.annotations;

    commands::normalize(
        counts_src,
        annotations_src,
        &options.feature_type,
        &options.id,
        options.method,
    )
}

fn main() -> anyhow::Result<()> {
    tracing_subscriber::fmt::init();

    let cli = Cli::parse();

    match cli.command {
        Command::Quantify(options) => quantify(options),
        Command::Normalize(options) => normalize(options),
    }
}
