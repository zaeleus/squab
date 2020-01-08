use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
    sync::Arc,
};

use clap::{crate_name, value_t, App, Arg};
use git_testament::{git_testament, render_testament};
use log::{info, LevelFilter};
use noodles::formats::bai;
use noodles_bam::{self as bam, Record};
use noodles_count_features::{
    count::count_paired_end_record_singletons, count::Filter, count_paired_end_records,
    count_single_end_records, read_features, Context, Features,
};

git_testament!(TESTAMENT);

fn write_counts<W>(
    writer: &mut W,
    counts: &HashMap<String, u64>,
    feature_ids: &[String],
) -> io::Result<()>
where
    W: Write,
{
    for id in feature_ids {
        let count = counts.get(id).unwrap_or(&0);
        writeln!(writer, "{}\t{}", id, count)?;
    }

    Ok(())
}

fn write_stats<W>(writer: &mut W, ctx: &Context) -> io::Result<()>
where
    W: Write,
{
    writeln!(writer, "__no_feature\t{}", ctx.no_feature)?;
    writeln!(writer, "__ambiguous\t{}", ctx.ambiguous)?;
    writeln!(writer, "__too_low_aQual\t{}", ctx.low_quality)?;
    writeln!(writer, "__not_aligned\t{}", ctx.unmapped)?;
    writeln!(writer, "__alignment_not_unique\t{}", ctx.nonunique)?;
    Ok(())
}

fn is_paired_end<P>(src: P) -> io::Result<bool>
where
    P: AsRef<Path>,
{
    let file = File::open(src)?;
    let mut reader = bam::Reader::new(file);
    reader.meta()?;

    let mut record = Record::new();
    reader.read_record(&mut record)?;

    Ok(record.flag().is_paired())
}

async fn count_single_end_records_by_region<P>(
    bam_src: P,
    ref_id: usize,
    reference: bam::Reference,
    features: Arc<Features>,
    references: Arc<Vec<bam::Reference>>,
    filter: Filter,
    strand_irrelevant: bool,
) -> Context
where
    P: AsRef<Path>,
{
    let file = File::open(bam_src.as_ref()).unwrap();
    let mut reader = bam::Reader::new(file);

    let bai_src = bam_src.as_ref().with_extension("bam.bai");
    let bai_file = File::open(bai_src).unwrap();
    let mut bai_reader = bai::Reader::new(bai_file);
    bai_reader.header().unwrap();
    let index = bai_reader.read_index().unwrap();

    let index_ref = &index.references[ref_id];

    let start = 0;
    let end = reference.len();

    let query = reader.query(index_ref, start, end);

    count_single_end_records(query, &features, &references, &filter, strand_irrelevant).unwrap()
}

#[tokio::main]
async fn main() {
    let matches = App::new(crate_name!())
        .version(render_testament!(TESTAMENT).as_str())
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .long("verbose")
                .help("Use verbose logging"),
        )
        .arg(
            Arg::with_name("with-secondary-records")
                .long("with-secondary-records")
                .help("Count secondary records (BAM flag 0x100)"),
        )
        .arg(
            Arg::with_name("with-supplementary-records")
                .long("with-supplementary-records")
                .help("Count supplementary records (BAM flag 0x800)"),
        )
        .arg(
            Arg::with_name("with-nonunique-records")
                .long("with-nonunique-records")
                .help("Count nonunique records (BAM data tag NH > 1)"),
        )
        .arg(
            Arg::with_name("strand-irrelevant")
                .long("strand-irrelevant")
                .help("Whether the sequencing protocol lacks strandedness"),
        )
        .arg(
            Arg::with_name("type")
                .short("t")
                .long("type")
                .value_name("str")
                .help("Feature type to count")
                .default_value("exon"),
        )
        .arg(
            Arg::with_name("id")
                .short("i")
                .long("id")
                .value_name("str")
                .help("Feature attribute to use as the feature identity")
                .default_value("gene_id"),
        )
        .arg(
            Arg::with_name("min-mapq")
                .long("min-mapq")
                .value_name("u8")
                .help("Minimum mapping quality to consider an alignment")
                .default_value("10"),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .value_name("file")
                .help("Output destination for feature counts")
                .required(true),
        )
        .arg(
            Arg::with_name("annotations")
                .short("a")
                .long("annotations")
                .value_name("file")
                .help("Input annotations file (GTF/GFFv2)")
                .required(true),
        )
        .arg(
            Arg::with_name("bam")
                .help("Input alignment file")
                .required(true)
                .index(1),
        )
        .get_matches();

    if matches.is_present("verbose") {
        env_logger::Builder::from_default_env()
            .filter(Some("noodles_count_features"), LevelFilter::Info)
            .init();
    } else {
        env_logger::init();
    }

    let bam_src = matches.value_of("bam").unwrap();
    let annotations_src = matches.value_of("annotations").unwrap();

    let results_dst = matches.value_of("output").unwrap();

    let feature_type = matches.value_of("type").unwrap();
    let id = matches.value_of("id").unwrap();
    let min_mapq = value_t!(matches, "min-mapq", u8).unwrap_or_else(|e| e.exit());

    let with_secondary_records = matches.is_present("with-secondary-records");
    let with_supplementary_records = matches.is_present("with-supplementary-records");
    let with_nonunique_records = matches.is_present("with-nonunique-records");
    let strand_irrelevant = matches.is_present("strand-irrelevant");

    let (features, names) = read_features(annotations_src, feature_type, id).unwrap();

    let file = File::open(&bam_src).unwrap();
    let mut reader = bam::Reader::new(file);
    let (_, references) = reader.meta().expect("failed to read bam header");

    let mut feature_ids = Vec::with_capacity(names.len());
    feature_ids.extend(names.into_iter());
    feature_ids.sort();

    let is_paired_end = is_paired_end(&bam_src).unwrap();

    let filter = Filter::new(
        min_mapq,
        with_secondary_records,
        with_supplementary_records,
        with_nonunique_records,
    );

    let ctx = if is_paired_end {
        info!("counting features for paired end records");

        let records = reader.records();
        let (mut ctx1, mut pairs) =
            count_paired_end_records(records, &features, &references, &filter, strand_irrelevant)
                .unwrap();

        let singletons = pairs.singletons().map(Ok);
        let ctx2 = count_paired_end_record_singletons(
            singletons,
            &features,
            &references,
            &filter,
            strand_irrelevant,
        )
        .unwrap();

        ctx1.add(&ctx2);

        ctx1
    } else {
        info!("counting features for single end records");

        let features = Arc::new(features);
        let references = Arc::new(references);

        let tasks: Vec<_> = references
            .iter()
            .enumerate()
            .map(|(ref_id, reference)| {
                tokio::spawn(count_single_end_records_by_region(
                    bam_src.to_string(),
                    ref_id,
                    reference.clone(),
                    features.clone(),
                    references.clone(),
                    filter.clone(),
                    strand_irrelevant,
                ))
            })
            .collect();

        let mut ctx = Context::default();

        for task in tasks {
            let region_ctx = task.await.unwrap();
            ctx.add(&region_ctx);
        }

        ctx
    };

    info!("writing counts");

    let file = File::create(results_dst).unwrap();
    let mut writer = BufWriter::new(file);
    write_counts(&mut writer, &ctx.counts, &feature_ids).unwrap();
    write_stats(&mut writer, &ctx).unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_counts() {
        let counts: HashMap<String, u64> = [
            (String::from("AADAT"), 302),
            (String::from("CLN3"), 37),
            (String::from("PAK4"), 145),
        ]
        .iter()
        .cloned()
        .collect();

        let ids = vec![
            String::from("AADAT"),
            String::from("CLN3"),
            String::from("NEO1"),
            String::from("PAK4"),
        ];

        let mut buf = Vec::new();

        write_counts(&mut buf, &counts, &ids).unwrap();

        let actual = String::from_utf8(buf).unwrap();

        let expected = "\
AADAT\t302
CLN3\t37
NEO1\t0
PAK4\t145
";

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_write_stats() {
        let mut buf = Vec::new();

        let mut ctx = Context::default();
        ctx.no_feature = 735;
        ctx.ambiguous = 5;
        ctx.low_quality = 60;
        ctx.unmapped = 8;
        ctx.nonunique = 13;

        write_stats(&mut buf, &ctx).unwrap();

        let actual = String::from_utf8(buf).unwrap();

        let expected = "\
__no_feature\t735
__ambiguous\t5
__too_low_aQual\t60
__not_aligned\t8
__alignment_not_unique\t13
";

        assert_eq!(actual, expected);
    }
}
