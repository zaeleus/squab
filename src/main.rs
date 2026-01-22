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
use tracing::info;

fn quantify(options: cli::Quantify) -> anyhow::Result<()> {
    let bam_src = options.src;
    let annotations_src = options.annotations;

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

fn main() -> anyhow::Result<()> {
    tracing_subscriber::fmt().with_writer(io::stderr).init();

    let cli = Cli::parse();

    match cli.command {
        Command::Quantify(options) => quantify(options),
    }
}
