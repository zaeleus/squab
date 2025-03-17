mod common;

use std::{env, fs, num::NonZero};

use noodles::sam::alignment::record::MappingQuality;
use squab::{StrandSpecificationOption, commands::quantify, count::Filter};

use self::common::convert_sam_to_bam;

const MIN_MAPPING_QUALITY: MappingQuality = MappingQuality::new(10).unwrap();

#[test]
fn test_quantify_with_single_end_forward_sample() -> anyhow::Result<()> {
    let feature_type = "exon";
    let id = "gene_id";
    let filter = Filter::new(MIN_MAPPING_QUALITY, false, false, false);
    let strand_specification_option = StrandSpecificationOption::Auto;
    let worker_count = NonZero::<usize>::MIN;

    let working_prefix = env::temp_dir();
    fs::create_dir_all(&working_prefix)?;

    let src = working_prefix.join("in.bam");
    convert_sam_to_bam("tests/fixtures/sample.single-end.forward-strand.sam", &src)?;

    let dst = working_prefix.join("out.txt");

    quantify(
        src,
        "tests/fixtures/annotations.gff3",
        feature_type,
        id,
        filter,
        strand_specification_option,
        worker_count,
        Some(&dst),
    )?;

    let actual = fs::read_to_string(dst)?;

    let expected = "\
sq0g0\t4
sq1g0\t2
sq2g0\t0
sq2g1\t0
__no_feature\t5
__ambiguous\t1
__too_low_aQual\t1
__not_aligned\t1
__alignment_not_unique\t1
";

    assert_eq!(actual, expected);

    Ok(())
}
