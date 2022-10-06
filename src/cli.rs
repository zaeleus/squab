use std::path::PathBuf;

use clap::{Parser, Subcommand};
use git_testament::{git_testament, render_testament};
use noodles::sam::record::MappingQuality;

use crate::{normalization::Method, StrandSpecificationOption};

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
    #[arg(short = 't', long, default_value_t = String::from("exon"))]
    pub feature_type: String,

    /// Feature attribute to use as the feature identity.
    #[arg(short = 'i', long, default_value_t = String::from("gene_id"))]
    pub id: String,

    /// Input annotations file (GFF3).
    #[arg(short = 'a', long)]
    pub annotations: PathBuf,

    /// Quantification normalization method.
    #[arg(long, default_value_t = Method::Tpm)]
    pub method: Method,

    /// Input counts file.
    pub counts: PathBuf,
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
    #[arg(long, default_value = "auto")]
    pub strand_specification: StrandSpecificationOption,

    /// Feature type to count.
    #[arg(short = 't', long, default_value_t = String::from("exon"))]
    pub feature_type: String,

    /// Feature attribute to use as the feature identity.
    #[arg(short = 'i', long, default_value_t = String::from("gene_id"))]
    pub id: String,

    #[arg(long, default_value = "10")]
    pub min_mapping_quality: MappingQuality,

    /// Output destination for feature counts.
    #[arg(short = 'o', long)]
    pub output: PathBuf,

    /// Input annotations file (GFF3).
    #[arg(short = 'a', long)]
    pub annotations: PathBuf,

    /// Force a specific number of threads.
    #[arg(long)]
    pub threads: Option<usize>,

    /// Input alignment file.
    pub bam: PathBuf,
}
