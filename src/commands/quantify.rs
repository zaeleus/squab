use std::{
    fs::File,
    io::{self, BufWriter, Write},
    num::NonZeroUsize,
    path::Path,
};

use anyhow::Context as _;
use noodles::{bam, bgzf};
use tracing::{info, warn};

use crate::{
    Context, StrandSpecification, StrandSpecificationOption, build_interval_trees,
    count::{Filter, context::Counts, count_paired_end_records, count_single_end_records},
    detect::{LibraryLayout, detect_specification},
    read_features,
};

#[allow(clippy::too_many_arguments)]
pub fn quantify<P, Q, R>(
    bam_src: P,
    annotations_src: Q,
    feature_type: &str,
    id: &str,
    filter: Filter,
    strand_specification_option: StrandSpecificationOption,
    worker_count: NonZeroUsize,
    dst: Option<R>,
) -> anyhow::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: AsRef<Path>,
{
    let mut gff_reader = crate::gff::open(annotations_src.as_ref())
        .with_context(|| format!("Could not open {}", annotations_src.as_ref().display()))?;

    let mut reader = bam::io::reader::Builder.build_from_path(bam_src.as_ref())?;
    let header = reader.read_header()?;

    info!(src = ?annotations_src.as_ref(), feature_type, feature_id = id, "reading features");

    let (reference_sequence_names, features) = read_features(&mut gff_reader, feature_type, id)?;

    info!(feature_count = features.len(), "read features");

    let interval_trees = build_interval_trees(&header, &reference_sequence_names, &features);

    let decoder: Box<dyn bgzf::io::Read + Send> = if worker_count.get() > 1 {
        File::open(bam_src.as_ref())
            .map(|f| bgzf::MultithreadedReader::with_worker_count(worker_count, f))
            .map(Box::new)
            .with_context(|| format!("Could not open {}", bam_src.as_ref().display()))?
    } else {
        bgzf::reader::Builder
            .build_from_path(bam_src.as_ref())
            .map(Box::new)
            .with_context(|| format!("Could not open {}", bam_src.as_ref().display()))?
    };

    let mut reader = bam::io::Reader::from(decoder);
    reader.read_header()?;

    info!("detecting library type");

    let (library_layout, detected_strand_specification, strandedness_confidence) =
        detect_specification(&bam_src, &interval_trees)?;

    info!("detected library layout: {library_layout}");
    info!(
        "strand specification: {detected_strand_specification} (confidence: {strandedness_confidence:.2})"
    );

    let strand_specification = strand_specification_from_option_or(
        strand_specification_option,
        detected_strand_specification,
    );

    if strand_specification != detected_strand_specification {
        warn!(
            expected = ?strand_specification,
            actual = ?detected_strand_specification,
            "strand specification mismatch",
        );
    }

    info!("counting features");

    let ctx = match library_layout {
        LibraryLayout::SingleEnd => count_single_end_records(
            reader,
            &interval_trees,
            &filter,
            strand_specification,
            worker_count,
        )?,
        LibraryLayout::PairedEnd => count_paired_end_records(
            reader,
            &interval_trees,
            &filter,
            strand_specification,
            worker_count,
        )?,
    };

    let mut writer: Box<dyn Write> = if let Some(dst) = dst {
        File::create(dst.as_ref())
            .map(BufWriter::new)
            .map(Box::new)
            .with_context(|| format!("Could not open {}", dst.as_ref().display()))?
    } else {
        let stdout = io::stdout().lock();
        Box::new(BufWriter::new(stdout))
    };

    let mut feature_names: Vec<_> = features.keys().collect();
    feature_names.sort();

    info!("writing counts");

    write_counts(&mut writer, &feature_names, &ctx.counts)?;
    write_metadata(&mut writer, &ctx)?;

    Ok(())
}

fn strand_specification_from_option_or(
    strand_specification_option: StrandSpecificationOption,
    detected_strand_specification: StrandSpecification,
) -> StrandSpecification {
    match strand_specification_option {
        StrandSpecificationOption::None => StrandSpecification::None,
        StrandSpecificationOption::Forward => StrandSpecification::Forward,
        StrandSpecificationOption::Reverse => StrandSpecification::Reverse,
        StrandSpecificationOption::Auto => detected_strand_specification,
    }
}

const DELIMITER: char = '\t';

fn write_counts<W>(writer: &mut W, feature_names: &[&String], counts: &Counts) -> io::Result<()>
where
    W: Write,
{
    const MISSING: u64 = 0;

    for name in feature_names {
        let count = counts.get(name.as_str()).copied().unwrap_or(MISSING);
        writeln!(writer, "{name}{DELIMITER}{count}")?;
    }

    Ok(())
}

fn write_metadata<W>(writer: &mut W, ctx: &Context) -> io::Result<()>
where
    W: Write,
{
    writeln!(writer, "__no_feature{DELIMITER}{}", ctx.miss)?;
    writeln!(writer, "__ambiguous{DELIMITER}{}", ctx.ambiguous)?;
    writeln!(writer, "__too_low_aQual{DELIMITER}{}", ctx.low_quality)?;
    writeln!(writer, "__not_aligned{DELIMITER}{}", ctx.unmapped)?;
    writeln!(writer, "__alignment_not_unique{DELIMITER}{}", ctx.nonunique)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strand_specification_from_option_or() {
        assert_eq!(
            strand_specification_from_option_or(
                StrandSpecificationOption::None,
                StrandSpecification::Forward
            ),
            StrandSpecification::None
        );

        assert_eq!(
            strand_specification_from_option_or(
                StrandSpecificationOption::Forward,
                StrandSpecification::Reverse
            ),
            StrandSpecification::Forward
        );

        assert_eq!(
            strand_specification_from_option_or(
                StrandSpecificationOption::Reverse,
                StrandSpecification::None
            ),
            StrandSpecification::Reverse
        );

        assert_eq!(
            strand_specification_from_option_or(
                StrandSpecificationOption::Auto,
                StrandSpecification::None
            ),
            StrandSpecification::None
        );
    }

    #[test]
    fn test_write_counts() -> io::Result<()> {
        let counts = [("AADAT", 302), ("CLN3", 37), ("PAK4", 145)]
            .into_iter()
            .collect();

        let names = [
            &String::from("AADAT"),
            &String::from("CLN3"),
            &String::from("NEO1"),
            &String::from("PAK4"),
        ];

        let mut buf = Vec::new();
        write_counts(&mut buf, &names, &counts)?;

        let expected = b"AADAT\t302\nCLN3\t37\nNEO1\t0\nPAK4\t145\n";
        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_metadata() -> io::Result<()> {
        let ctx = Context {
            miss: 735,
            ambiguous: 5,
            low_quality: 60,
            unmapped: 8,
            nonunique: 13,
            ..Default::default()
        };

        let mut buf = Vec::new();
        write_metadata(&mut buf, &ctx)?;

        let expected = b"__no_feature\t735
__ambiguous\t5
__too_low_aQual\t60
__not_aligned\t8
__alignment_not_unique\t13
";

        assert_eq!(buf, expected);

        Ok(())
    }
}
