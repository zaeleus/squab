# noodles count-features

**noodles count-features** counts the number of aligned records that intersects
a set of features.

The goal of this project is to provide a faster alternative to [htseq-count].
It uses the same counting rules and outputs a compatible data table.

[htseq-count]: https://htseq.readthedocs.io/en/master/count.html

## Installation

Use [Cargo] to install `noodles-count-features`.

```
$ cargo install --git https://github.com/zaeleus/noodles-count-features.git
```

[Cargo]: https://doc.rust-lang.org/cargo/getting-started/installation.html

## Usage

```
noodles-count-features 0.1.0

USAGE:
    noodles-count-features [FLAGS] [OPTIONS] <bam> --annotations <file> --output <file>

FLAGS:
    -h, --help                          Prints help information
    -V, --version                       Prints version information
    -v, --verbose                       Use verbose logging
        --with-secondary-records        Count secondary records (BAM flag 0x100)
        --with-supplementary-records    Count supplementary records (BAM flag 0x800)

OPTIONS:
    -a, --annotations <file>    Input annotations file (GTF/GFFv2)
    -i, --id <str>              Feature attribute to use as the feature identity [default: gene_id]
        --min-mapq <u8>         Minimum mapping quality to consider an alignment [default: 10]
    -o, --output <file>         Output destination for feature counts
    -t, --type <str>            Feature type to count [default: exon]

ARGS:
    <bam>    Input alignment file
```

The output is a tab-delimited text file with two columns: the feature
identifier and the number of reads from the input alignment that overlap it.
This file is compatible as output from htseq-count, meaning it includes
statistics in the trailer.

## Limitations

  * Counts are taken only as the union of matched feature sets, i.e., reads that
    overlap any part of the feature is considered once.
  * Nonunique reads (BAM data tag `NH` > 1) are always counted.
  * For paired end alignments, a read that matches itself before a mate is
    found replaces the previously known record.
