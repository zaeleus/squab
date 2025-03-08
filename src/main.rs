use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::{io, num::NonZero, thread};

use clap::Parser;
use squab::{
    Cli,
    cli::{self, Command},
    commands,
    count::Filter,
};
use tracing::{info, warn};

const DEFAULT_FEATURE_ID: &str = "gene_id";

fn quantify(mut options: cli::Quantify) -> anyhow::Result<()> {
    let bam_src = options.src;
    let annotations_src = options.annotations;

    if options.id != DEFAULT_FEATURE_ID {
        warn!("The --id option is deprecated. Use `--feature-id` instead.");
        options.feature_id = options.id.clone();
    }

    if let Some(threads) = options.threads {
        warn!("The --threads option is deprecated. Use `--worker-count` instead.");
        options.worker_count = Some(threads);
    }

    let worker_count = options
        .worker_count
        .unwrap_or_else(|| thread::available_parallelism().unwrap_or(NonZero::<usize>::MIN));

    let strand_specification_option = options.strand_specification;

    let filter = Filter::new(
        options.min_mapping_quality,
        options.with_secondary_records,
        options.with_supplementary_records,
        options.with_nonunique_records,
    );

    info!("using {} thread(s)", worker_count);

    commands::quantify(
        bam_src,
        annotations_src,
        &options.feature_type,
        &options.feature_id,
        filter,
        strand_specification_option,
        worker_count,
        options.output,
    )?;

    Ok(())
}

fn normalize(mut options: cli::Normalize) -> anyhow::Result<()> {
    let counts_src = options.src;
    let annotations_src = options.annotations;

    if options.id != DEFAULT_FEATURE_ID {
        warn!("The --id option is deprecated. Use `--feature-id` instead.");
        options.feature_id = options.id.clone();
    }

    commands::normalize(
        counts_src,
        annotations_src,
        &options.feature_type,
        &options.feature_id,
        options.method,
        options.output,
    )?;

    Ok(())
}

fn main() -> anyhow::Result<()> {
    tracing_subscriber::fmt().with_writer(io::stderr).init();

    let cli = Cli::parse();

    match cli.command {
        Command::Quantify(options) => quantify(options),
        Command::Normalize(options) => normalize(options),
    }
}
