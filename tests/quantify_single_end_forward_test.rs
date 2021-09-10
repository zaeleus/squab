use std::{env, fs};

use noodles_squab::{commands::quantify, count::Filter, StrandSpecificationOption};

#[tokio::test]
async fn test_quantify_with_single_end_forward_sample() -> anyhow::Result<()> {
    let feature_type = "exon";
    let id = "gene_id";
    let filter = Filter::new(10, false, false, false);
    let strand_specification_option = StrandSpecificationOption::Auto;

    let working_prefix = env::temp_dir();
    fs::create_dir_all(&working_prefix)?;

    let results_dst = working_prefix.join("out.txt");

    quantify(
        "tests/fixtures/sample.single-end.forward-strand.bam",
        "tests/fixtures/annotations.gff3",
        feature_type,
        id,
        filter,
        strand_specification_option,
        &results_dst,
    )
    .await?;

    let actual = fs::read_to_string(results_dst)?;

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
