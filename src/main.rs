#[macro_use] extern crate clap;
extern crate env_logger;
extern crate noodles_count_features;
extern crate interval_tree;
#[macro_use] extern crate log;
extern crate noodles;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::ops::Range;
use std::path::Path;

use clap::{App, Arg};
use interval_tree::IntervalTree;
use log::LevelFilter;
use noodles::formats::bam::{self, ByteRecord, Flag, Reference};
use noodles_count_features::{Entry, Features, RecordPairs, Strand, cigar_to_intervals, read_features};

#[derive(Default)]
struct Context {
    counts: HashMap<String, u64>,
    no_feature: u64,
    ambiguous: u64,
    low_quality: u64,
    unmapped: u64,
}

fn count_single_end_records<R>(
    mut reader: bam::Reader<R>,
    features: Features,
    references: Vec<Reference>,
    min_mapq: u8,
) -> io::Result<Context>
where
    R: Read,
{
    let mut ctx = Context::default();
    let mut record = ByteRecord::new();

    loop {
        let bytes_read = reader.read_byte_record(&mut record)?;

        if bytes_read == 0 {
            break;
        }

        let flag = Flag::new(record.flag());

        if flag.is_unmapped() {
            ctx.unmapped += 1;
            continue;
        }

        if record.mapq() < min_mapq {
            ctx.low_quality += 1;
            continue;
        }

        if record.ref_id() < 0 {
            panic!("no ref id");
        }

        let ref_id = record.ref_id() as usize;
        let reference = &references[ref_id];

        let intervals = cigar_to_intervals(&record, false);

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            },
        };

        let set = find(tree, &intervals);

        if set.len() == 0 {
            ctx.no_feature += 1;
        } else if set.len() == 1 {
            for gene_name in set {
                let count = ctx.counts.entry(gene_name.clone()).or_insert(0);
                *count += 1;
            }
        } else if set.len() > 1 {
            ctx.ambiguous += 1;
        }
    }

    Ok(ctx)
}

fn count_paired_end_records<R>(
    reader: bam::Reader<R>,
    features: Features,
    references: Vec<Reference>,
    min_mapq: u8,
) -> io::Result<Context>
where
    R: Read,
{
    let mut ctx = Context::default();
    let mut pairs = RecordPairs::new(reader);

    for pair in &mut pairs {
        let (r1, r2) = pair?;

        let f1 = Flag::new(r1.flag());
        let f2 = Flag::new(r2.flag());

        if f1.is_unmapped() && f2.is_unmapped() {
            ctx.unmapped += 1;
            continue;
        }

        if r1.mapq() < min_mapq || r2.mapq() < min_mapq {
            ctx.low_quality += 1;
            continue;
        }

        let ref_id = r1.ref_id();

        if ref_id < 0 {
            panic!("no ref id");
        }

        let reference = &references[ref_id as usize];

        let intervals = cigar_to_intervals(&r1, false);

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            },
        };

        let mut set = find(tree, &intervals);

        let ref_id = r2.ref_id();

        if ref_id < 0 {
            panic!("no ref id");
        }

        let reference = &references[ref_id as usize];

        let intervals = cigar_to_intervals(&r2, true);

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            },
        };

        let set2 = find(tree, &intervals);

        set.extend(set2.into_iter());

        if set.len() == 0 {
            ctx.no_feature += 1;
        } else if set.len() == 1 {
            for gene_name in set {
                let count = ctx.counts.entry(gene_name.clone()).or_insert(0);
                *count += 1;
            }
        } else if set.len() > 1 {
            ctx.ambiguous += 1;
        }
    }

    let rest = pairs.orphan_pairs();

    for (r1, r2) in rest {
        let mut intervals = if let Some(ref r1) = r1 {
            let f1 = Flag::new(r1.flag());

            if f1.is_unmapped() {
                ctx.unmapped += 1;
                continue;
            }

            if r1.mapq() < min_mapq {
                ctx.low_quality += 1;
                continue;
            }

            cigar_to_intervals(&r1, false)
        } else {
            vec![]
        };

        if let Some(ref r2) = r2 {
            let f2 = Flag::new(r2.flag());

            if f2.is_unmapped() {
                ctx.unmapped += 1;
                continue;
            }

            if r2.mapq() < min_mapq {
                ctx.low_quality += 1;
                continue;
            }

            intervals.extend(cigar_to_intervals(&r2, true));
        }

        let ref_id = if let Some(ref r1) = r1 {
            r1.ref_id()
        } else if let Some(ref r2) = r2 {
            r2.ref_id()
        } else {
            unreachable!("one orphan record must be set");
        };

        if ref_id < 0 {
            panic!("no ref id");
        }

        let reference = &references[ref_id as usize];

        let name = reference.name();
        let tree = match features.get(name) {
            Some(t) => t,
            None => {
                ctx.no_feature += 1;
                continue;
            },
        };

        let set = find(tree, &intervals);

        if set.len() == 0 {
            ctx.no_feature += 1;
        } else if set.len() == 1 {
            for gene_name in set {
                let count = ctx.counts.entry(gene_name.clone()).or_insert(0);
                *count += 1;
            }
        } else if set.len() > 1 {
            ctx.ambiguous += 1;
        }
    }

    Ok(ctx)
}

fn find(
    tree: &IntervalTree<u64, Entry>,
    intervals: &[(Range<u64>, bool)],
) -> HashSet<String> {
    let mut set = HashSet::new();

    for (interval, is_reverse) in intervals {
        for entry in tree.find(interval.clone()) {
            let gene_name = &entry.value.0;
            let strand = &entry.value.1;

            if (strand == &Strand::Reverse && *is_reverse)
                    || (strand == &Strand::Forward && !*is_reverse) {
                set.insert(gene_name.to_string());
            }
        }
    }

    set
}

fn write_counts<P>(
    dst: P,
    counts: &HashMap<String, u64>,
    feature_ids: &[String],
) -> io::Result<()>
where
    P: AsRef<Path>,
{
    let file = File::create(dst)?;
    let mut writer = BufWriter::new(file);

    for id in feature_ids {
        let count = counts.get(id).unwrap_or(&0);
        writeln!(&mut writer, "{}\t{}", id, count)?;
    }

    Ok(())
}

fn main() {
    let matches = App::new(crate_name!())
        .version(crate_version!())
        .arg(Arg::with_name("verbose")
             .short("v")
             .long("verbose")
             .help("Use verbose logging"))
        .arg(Arg::with_name("type")
             .short("t")
             .long("type")
             .value_name("str")
             .help("Feature type to count")
             .default_value("exon"))
        .arg(Arg::with_name("id")
             .short("i")
             .long("id")
             .value_name("str")
             .help("Feature attribute to use as the feature identity")
             .default_value("gene_id"))
        .arg(Arg::with_name("min-mapq")
             .long("min-mapq")
             .value_name("u8")
             .help("Minimum mapping quality to consider an alignment")
             .default_value("10"))
        .arg(Arg::with_name("paired-end")
             .long("paired-end")
             .help("Read records as mate pairs"))
        .arg(Arg::with_name("output")
             .short("o")
             .long("output")
             .value_name("file")
             .help("Output destination for feature counts")
             .required(true))
        .arg(Arg::with_name("annotations")
             .short("a")
             .long("annotations")
             .value_name("file")
             .help("Input annotations file (GTF/GFFv2)")
             .required(true))
        .arg(Arg::with_name("bam")
             .help("Input alignment file")
             .required(true)
             .index(1))
        .get_matches();

    if matches.is_present("verbose") {
        env_logger::Builder::from_default_env()
            .filter(Some("feature_count"), LevelFilter::Info)
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

    let is_paired_end = matches.is_present("paired-end");

    let (features, names) = read_features(annotations_src, feature_type, id).unwrap();

    let mut reader = bam::Reader::<File>::open(bam_src).unwrap();

    let _header = reader.header().unwrap();
    let references: Vec<Reference> = reader.references()
        .unwrap()
        .filter_map(Result::ok)
        .collect();

    let mut feature_ids = Vec::with_capacity(names.len());
    feature_ids.extend(names.into_iter());
    feature_ids.sort();

    let ctx = if is_paired_end {
        info!("counting features for paired end records");
        count_paired_end_records(reader, features, references, min_mapq).unwrap()
    } else {
        info!("counting features for single end records");
        count_single_end_records(reader, features, references, min_mapq).unwrap()
    };

    info!("writing counts");

    write_counts(results_dst, &ctx.counts, &feature_ids).unwrap();

    eprintln!("__no_feature: {}", ctx.no_feature);
    eprintln!("__ambiguous: {}", ctx.ambiguous);
    eprintln!("__too_low_aQual: {}", ctx.low_quality);
    eprintln!("__not_aligned: {}", ctx.unmapped);
    eprintln!("__alignment_not_unique: -");
}
