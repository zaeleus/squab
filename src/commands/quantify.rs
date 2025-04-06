use std::{
    fs::File,
    io::{self, BufWriter, Write},
    num::NonZero,
    path::{Path, PathBuf},
};

use bstr::{BString, ByteSlice};
use noodles::{bam, bgzf};
use thiserror::Error;
use tracing::{info, warn};

use crate::{
    Context, ReadFeaturesError, StrandSpecification, StrandSpecificationOption,
    build_interval_trees,
    count::{Filter, context::Counts, count_paired_end_records, count_single_end_records},
    detect::{LibraryLayout, detect_specification},
    read_features,
};

#[allow(clippy::too_many_arguments)]
pub fn quantify<P, Q, R>(
    src: P,
    annotations_src: Q,
    feature_type: &str,
    feature_id: &str,
    filter: Filter,
    strand_specification_option: StrandSpecificationOption,
    worker_count: NonZero<usize>,
    dst: Option<R>,
) -> Result<(), QuantifyError>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: AsRef<Path>,
{
    let src = src.as_ref();
    let annotations_src = annotations_src.as_ref();

    let mut gff_reader = crate::gff::open(annotations_src)
        .map_err(|e| QuantifyError::OpenFile(e, annotations_src.into()))?;

    let mut reader = bam::io::reader::Builder
        .build_from_path(src)
        .map_err(|e| QuantifyError::OpenFile(e, src.into()))?;

    let header = reader.read_header()?;

    info!(src = ?annotations_src, feature_type, feature_id, "reading features");

    let (reference_sequence_names, features) =
        read_features(&mut gff_reader, feature_type, feature_id)?;

    info!(feature_count = features.len(), "read features");

    let interval_trees = build_interval_trees(&header, &reference_sequence_names, &features);

    let decoder: Box<dyn bgzf::io::Read + Send> = if worker_count.get() > 1 {
        File::open(src)
            .map(|f| bgzf::MultithreadedReader::with_worker_count(worker_count, f))
            .map(Box::new)
            .map_err(|e| QuantifyError::OpenFile(e, src.into()))?
    } else {
        bgzf::io::reader::Builder
            .build_from_path(src)
            .map(Box::new)
            .map_err(|e| QuantifyError::OpenFile(e, src.into()))?
    };

    let mut reader = bam::io::Reader::from(decoder);
    reader.read_header()?;

    info!("detecting library type");

    let (library_layout, detected_strand_specification, strandedness_confidence) =
        detect_specification(src, &interval_trees)?;

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
        let dst = dst.as_ref();

        File::create(dst)
            .map(BufWriter::new)
            .map(Box::new)
            .map_err(|e| QuantifyError::CreateFile(e, dst.into()))?
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

#[derive(Debug, Error)]
pub enum QuantifyError {
    #[error("I/O error")]
    Io(#[from] io::Error),
    #[error("could not open file: {1}")]
    OpenFile(#[source] io::Error, PathBuf),
    #[error("could not create file: {1}")]
    CreateFile(#[source] io::Error, PathBuf),
    #[error("invalid annotations")]
    ReadAnnotations(#[from] ReadFeaturesError),
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

fn write_counts<W>(writer: &mut W, feature_names: &[&BString], counts: &Counts) -> io::Result<()>
where
    W: Write,
{
    const MISSING: u64 = 0;

    for name in feature_names {
        let count = counts.get(name.as_bstr()).copied().unwrap_or(MISSING);
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
    use bstr::BStr;

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
        let counts = [
            (BStr::new("AADAT"), 302),
            (BStr::new("CLN3"), 37),
            (BStr::new("PAK4"), 145),
        ]
        .into_iter()
        .collect();

        let names = [
            &BString::from("AADAT"),
            &BString::from("CLN3"),
            &BString::from("NEO1"),
            &BString::from("PAK4"),
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
