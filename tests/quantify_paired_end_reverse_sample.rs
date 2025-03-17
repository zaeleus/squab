mod common;

use std::{env, fs, num::NonZero};

use noodles::sam::alignment::record::MappingQuality;
use squab::{StrandSpecificationOption, commands::quantify, count::Filter};

use self::common::convert_sam_to_bam;

const MIN_MAPPING_QUALITY: MappingQuality = MappingQuality::new(10).unwrap();

#[test]
fn test_quantify_paired_end_reverse_sample() -> anyhow::Result<()> {
    let feature_type = "exon";
    let feature_id = "gene_id";
    let filter = Filter::new(MIN_MAPPING_QUALITY, false, false, false);
    let strand_specification_option = StrandSpecificationOption::Reverse;
    let worker_count = NonZero::<usize>::MIN;

    let working_prefix = env::temp_dir();
    fs::create_dir_all(&working_prefix)?;

    let src = working_prefix.join("in.bam");
    convert_sam_to_bam("tests/fixtures/sample.paired-end.reverse-strand.sam", &src)?;

    let dst = working_prefix.join("out.txt");

    quantify(
        src,
        "tests/fixtures/annotations.gff3",
        feature_type,
        feature_id,
        filter,
        strand_specification_option,
        worker_count,
        Some(&dst),
    )?;

    let actual = fs::read_to_string(dst)?;

    let expected = "\
sq0g0\t1
sq1g0\t1
sq2g0\t0
sq2g1\t1
__no_feature\t2
__ambiguous\t1
__too_low_aQual\t1
__not_aligned\t1
__alignment_not_unique\t1
";

    assert_eq!(actual, expected);

    Ok(())
}
