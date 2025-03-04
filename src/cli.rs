use std::{num::NonZero, path::PathBuf};

use clap::{Parser, Subcommand};
use git_testament::{git_testament, render_testament};
use noodles::sam::alignment::record::MappingQuality;

use crate::{StrandSpecificationOption, normalization::Method};

git_testament!(TESTAMENT);

#[derive(Parser)]
#[command(version = render_testament!(TESTAMENT))]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand)]
pub enum Command {
    /// Normalize features counts.
    Normalize(Normalize),
    /// Gene expression quantification.
    Quantify(Quantify),
}

#[derive(Parser)]
pub struct Normalize {
    /// Feature type to count.
    #[arg(short = 't', long, default_value = "exon")]
    pub feature_type: String,

    /// Feature attribute to use as the feature identity.
    #[arg(short = 'i', long, default_value = "gene_id")]
    pub id: String,

    /// Input annotations file (GFF3).
    #[arg(short = 'a', long)]
    pub annotations: PathBuf,

    /// Quantification normalization method.
    #[arg(long, value_enum, default_value_t = Method::Tpm)]
    pub method: Method,

    /// Output destination.
    ///
    /// If not set, output is written to stdout.
    #[arg(short = 'o', long)]
    pub output: Option<PathBuf>,

    /// Input counts file.
    pub src: PathBuf,
}

#[derive(Parser)]
pub struct Quantify {
    /// Count secondary records (BAM flag 0x100).
    #[arg(long)]
    pub with_secondary_records: bool,

    /// Count supplementary records (BAM flag 0x800).
    #[arg(long)]
    pub with_supplementary_records: bool,

    /// Count nonunique records (BAM data tag NH > 1).
    #[arg(long)]
    pub with_nonunique_records: bool,

    /// Strand specification.
    #[arg(long, value_enum, default_value_t = StrandSpecificationOption::Auto)]
    pub strand_specification: StrandSpecificationOption,

    /// Feature type to count.
    #[arg(short = 't', long, default_value = "exon")]
    pub feature_type: String,

    /// Feature attribute to use as the feature identity.
    #[arg(short = 'i', long, default_value = "gene_id")]
    pub id: String,

    #[arg(long, value_parser = parse_mapping_quality, default_value = "10")]
    pub min_mapping_quality: MappingQuality,

    /// Output destination.
    ///
    /// If not set, output is written to stdout.
    #[arg(short = 'o', long)]
    pub output: Option<PathBuf>,

    /// Input annotations file (GFF3).
    #[arg(short = 'a', long)]
    pub annotations: PathBuf,

    /// Force a specific number of threads.
    #[arg(long)]
    pub threads: Option<NonZero<usize>>,

    /// The number of workers to spawn.
    ///
    /// By default, this (usually) uses the number of available CPUs.
    #[arg(long)]
    pub worker_count: Option<NonZero<usize>>,

    /// Input alignment file.
    pub src: PathBuf,
}

fn parse_mapping_quality(s: &str) -> Result<MappingQuality, &'static str> {
    s.parse::<u8>()
        .map_err(|_| "invalid input")
        .and_then(|n| MappingQuality::new(n).ok_or("missing mapping quality"))
}
