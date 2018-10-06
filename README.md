# noodles count-features

**noodles count-features** counts how many aligned reads map to a list of features.

## Usage

By default, the input BAM is assumed to be single-end reads. Use `--paired-end`
if the input is paired-end reads.

## Limitations

  * Counts are taken only as the union of matched feature sets, i.e., reads that
    overlap any part of the feature is considered once.
  * Secondary and supplementary reads (BAM flags 0x100 and 0x800) are always
    counted.
  * Nonunique reads (BAM data tag `NH` > 1) are always counted.
  * For paired end alignments, a read that matches itself before a mate is
    found replaces the previously known record.
