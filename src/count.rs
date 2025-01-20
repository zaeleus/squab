pub mod context;
mod filter;
pub(crate) mod reader;
mod writer;

pub use self::{context::Context, filter::Filter, reader::Reader, writer::Writer};

use std::{
    collections::HashSet,
    io::{self, Read},
    num::NonZeroUsize,
    thread,
};

use crate::{
    Entry, IntervalTrees, MatchIntervals, RecordPairs, SegmentPosition, StrandSpecification,
};
use interval_tree::IntervalTree;
use noodles::{bam, core::Position, gff};

use self::context::Event;

const CHUNK_SIZE: usize = 8192;

pub fn count_single_end_records<R>(
    mut reader: bam::io::Reader<R>,
    interval_trees: &IntervalTrees,
    filter: &Filter,
    strand_specification: StrandSpecification,
    worker_count: NonZeroUsize,
) -> io::Result<Context>
where
    R: Read + Send,
{
    thread::scope(move |scope| {
        let (tx, rx) = crossbeam_channel::bounded(worker_count.get());

        scope.spawn(move || {
            let mut records = reader.records();

            loop {
                let mut chunk = Vec::with_capacity(CHUNK_SIZE);

                for result in records.by_ref().take(CHUNK_SIZE) {
                    let record = result?;
                    chunk.push(record);
                }

                if chunk.is_empty() {
                    drop(tx);
                    break;
                } else {
                    tx.send(chunk).expect("reader channel unexpectedly closed");
                }
            }

            Ok::<_, io::Error>(())
        });

        let handles: Vec<_> = (0..worker_count.get())
            .map(|_| {
                let rx = rx.clone();

                scope.spawn(move || {
                    let mut ctx = Context::default();

                    while let Ok(chunk) = rx.recv() {
                        for record in chunk {
                            let event = count_single_end_record(
                                interval_trees,
                                filter,
                                strand_specification,
                                &record,
                            )?;

                            ctx.add_event(event);
                        }
                    }

                    Ok::<_, io::Error>(ctx)
                })
            })
            .collect();

        let mut ctx = Context::default();

        for handle in handles {
            let c = handle.join().unwrap()?;
            ctx.add(&c);
        }

        Ok(ctx)
    })
}

pub fn count_single_end_record(
    interval_trees: &IntervalTrees,
    filter: &Filter,
    strand_specification: StrandSpecification,
    record: &bam::Record,
) -> io::Result<Event> {
    if let Some(event) = filter.filter(record)? {
        return Ok(event);
    }

    let reference_sequence_id = record
        .reference_sequence_id()
        .transpose()?
        .expect("missing reference sequence ID");

    let Some(interval_tree) = interval_trees.get(reference_sequence_id) else {
        return Ok(Event::NoFeature);
    };

    let cigar = record.cigar();
    let start = record.alignment_start().expect("missing alignment start")?;
    let intervals = MatchIntervals::new(&cigar, start);

    let flags = record.flags();
    let is_reverse = match strand_specification {
        StrandSpecification::Reverse => !flags.is_reverse_complemented(),
        _ => flags.is_reverse_complemented(),
    };

    let set = find(interval_tree, intervals, strand_specification, is_reverse)?;

    Ok(update_intersections(set))
}

pub fn count_paired_end_records<R>(
    reader: bam::io::Reader<R>,
    interval_trees: &IntervalTrees,
    filter: &Filter,
    strand_specification: StrandSpecification,
    worker_count: NonZeroUsize,
) -> io::Result<Context>
where
    R: Read + Send,
{
    let primary_only = !filter.with_secondary_records() && !filter.with_supplementary_records();
    let mut record_pairs = RecordPairs::new(reader, primary_only);

    let (mut ctx, mut record_pairs) = thread::scope(move |scope| {
        let (tx, rx) = crossbeam_channel::bounded(worker_count.get());

        let reader_handle = scope.spawn(move || {
            loop {
                let mut chunk = Vec::with_capacity(CHUNK_SIZE);

                for _ in 0..CHUNK_SIZE {
                    match record_pairs.next_pair()? {
                        Some(pair) => chunk.push(pair),
                        None => break,
                    }
                }

                if chunk.is_empty() {
                    drop(tx);
                    break;
                } else {
                    tx.send(chunk).expect("reader channel unexpectedly closed");
                }
            }

            Ok::<_, io::Error>(record_pairs)
        });

        let handles: Vec<_> = (0..worker_count.get())
            .map(|_| {
                let rx = rx.clone();

                scope.spawn(move || {
                    let mut ctx = Context::default();

                    while let Ok(chunk) = rx.recv() {
                        for (r1, r2) in chunk {
                            let event = count_paired_end_record_pair(
                                interval_trees,
                                filter,
                                strand_specification,
                                &r1,
                                &r2,
                            )?;

                            ctx.add_event(event);
                        }
                    }

                    Ok::<_, io::Error>(ctx)
                })
            })
            .collect();

        let record_pairs = reader_handle.join().unwrap()?;

        let mut ctx = Context::default();

        for handle in handles {
            let c = handle.join().unwrap()?;
            ctx.add(&c);
        }

        Ok::<_, io::Error>((ctx, record_pairs))
    })?;

    for record in record_pairs.singletons() {
        let event = count_paired_end_singleton_record(
            interval_trees,
            filter,
            strand_specification,
            &record,
        )?;

        ctx.add_event(event);
    }

    Ok(ctx)
}

pub fn count_paired_end_record_pair(
    interval_trees: &IntervalTrees,
    filter: &Filter,
    strand_specification: StrandSpecification,
    r1: &bam::Record,
    r2: &bam::Record,
) -> io::Result<Event> {
    if let Some(event) = filter.filter_pair(r1, r2)? {
        return Ok(event);
    }

    let reference_sequence_id = r1
        .reference_sequence_id()
        .transpose()?
        .expect("missing reference sequence ID");

    let Some(interval_tree) = interval_trees.get(reference_sequence_id) else {
        return Ok(Event::NoFeature);
    };

    let cigar = r1.cigar();
    let start = r1.alignment_start().expect("missing alignment start")?;
    let intervals = MatchIntervals::new(&cigar, start);

    let f1 = r1.flags();
    let is_reverse = match strand_specification {
        StrandSpecification::Reverse => !f1.is_reverse_complemented(),
        _ => f1.is_reverse_complemented(),
    };

    let mut set = find(interval_tree, intervals, strand_specification, is_reverse)?;

    let reference_sequence_id = r2
        .reference_sequence_id()
        .transpose()?
        .expect("missing reference sequence ID");

    let Some(interval_tree) = interval_trees.get(reference_sequence_id) else {
        return Ok(Event::NoFeature);
    };

    let cigar = r2.cigar();
    let start = r2.alignment_start().expect("missing alignment start")?;
    let intervals = MatchIntervals::new(&cigar, start);

    let f2 = r2.flags();
    let is_reverse = match strand_specification {
        StrandSpecification::Reverse => f2.is_reverse_complemented(),
        _ => !f2.is_reverse_complemented(),
    };

    let set2 = find(interval_tree, intervals, strand_specification, is_reverse)?;

    set.extend(set2);

    Ok(update_intersections(set))
}

pub fn count_paired_end_singleton_record(
    interval_trees: &IntervalTrees,
    filter: &Filter,
    strand_specification: StrandSpecification,
    record: &bam::Record,
) -> io::Result<Event> {
    if let Some(event) = filter.filter(record)? {
        return Ok(event);
    }

    let reference_sequence_id = record
        .reference_sequence_id()
        .transpose()?
        .expect("missing reference sequence ID");

    let Some(interval_tree) = interval_trees.get(reference_sequence_id) else {
        return Ok(Event::NoFeature);
    };

    let cigar = record.cigar();
    let start = record.alignment_start().expect("missing alignment start")?;
    let intervals = MatchIntervals::new(&cigar, start);

    let flags = record.flags();
    let is_reverse = match SegmentPosition::try_from(flags) {
        Ok(SegmentPosition::First) => match strand_specification {
            StrandSpecification::Reverse => !flags.is_reverse_complemented(),
            _ => flags.is_reverse_complemented(),
        },
        Ok(SegmentPosition::Last) => match strand_specification {
            StrandSpecification::Reverse => flags.is_reverse_complemented(),
            _ => !flags.is_reverse_complemented(),
        },
        Err(e) => return Err(io::Error::new(io::ErrorKind::InvalidData, e)),
    };

    let set = find(interval_tree, intervals, strand_specification, is_reverse)?;

    Ok(update_intersections(set))
}

fn find(
    tree: &IntervalTree<Position, Entry>,
    intervals: MatchIntervals,
    strand_specification: StrandSpecification,
    is_reverse: bool,
) -> io::Result<HashSet<String>> {
    let mut set = HashSet::new();

    for result in intervals {
        let interval = result?;

        for entry in tree.find(interval.clone()) {
            let (gene_name, strand) = entry.get();

            match strand_specification {
                StrandSpecification::None => {
                    set.insert(gene_name.to_string());
                }
                StrandSpecification::Forward | StrandSpecification::Reverse => {
                    if (strand == &gff::record::Strand::Reverse && is_reverse)
                        || (strand == &gff::record::Strand::Forward && !is_reverse)
                    {
                        set.insert(gene_name.to_string());
                    }
                }
            }
        }
    }

    Ok(set)
}

fn update_intersections(mut intersections: HashSet<String>) -> Event {
    if intersections.is_empty() {
        Event::NoFeature
    } else if intersections.len() == 1 {
        intersections.drain().next().map(Event::Hit).unwrap()
    } else {
        Event::Ambiguous
    }
}
