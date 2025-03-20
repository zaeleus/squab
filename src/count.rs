pub mod context;
mod filter;

pub use self::{context::Context, filter::Filter};

use std::{
    collections::HashSet,
    io::{self, Read},
    num::NonZero,
    thread,
};

use crate::{Entry, IntervalTrees, MatchIntervals, RecordPairs, StrandSpecification};
use interval_tree::IntervalTree;
use noodles::{bam, core::Position};
use tracing::warn;

use self::context::Event;

const CHUNK_SIZE: usize = 8192;

pub fn count_single_end_records<'f, R>(
    mut reader: bam::io::Reader<R>,
    interval_trees: &IntervalTrees<'f>,
    filter: &'f Filter,
    strand_specification: StrandSpecification,
    worker_count: NonZero<usize>,
) -> io::Result<Context<'f>>
where
    R: Read + Send,
{
    thread::scope(move |scope| {
        let (tx, rx) = crossbeam_channel::bounded(worker_count.get());

        scope.spawn(move || {
            let mut records = reader.records();

            loop {
                let chunk: Vec<_> = records
                    .by_ref()
                    .take(CHUNK_SIZE)
                    .collect::<io::Result<_>>()?;

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

pub fn count_single_end_record<'f>(
    interval_trees: &IntervalTrees<'f>,
    filter: &'f Filter,
    strand_specification: StrandSpecification,
    record: &bam::Record,
) -> io::Result<Event<'f>> {
    if let Some(event) = filter.filter(record)? {
        return Ok(event);
    }

    let mut intersections = HashSet::new();

    let is_reverse_complemented = resolve_is_reverse_complemented(
        record.flags().is_reverse_complemented(),
        strand_specification,
    );

    if let Some(event) = count_record(
        interval_trees,
        strand_specification,
        is_reverse_complemented,
        record,
        &mut intersections,
    )? {
        return Ok(event);
    }

    Ok(resolve_intersections(&intersections))
}

pub fn count_paired_end_records<'f, R>(
    reader: bam::io::Reader<R>,
    interval_trees: &IntervalTrees<'f>,
    filter: &'f Filter,
    strand_specification: StrandSpecification,
    worker_count: NonZero<usize>,
) -> io::Result<Context<'f>>
where
    R: Read + Send,
{
    let primary_only = !filter.with_secondary_records() && !filter.with_supplementary_records();
    let mut record_pairs = RecordPairs::new(reader, primary_only);

    let (mut ctx, record_pairs) = thread::scope(move |scope| {
        let (tx, rx) = crossbeam_channel::bounded(worker_count.get());

        let reader_handle = scope.spawn(move || {
            loop {
                let chunk: Vec<_> = record_pairs
                    .by_ref()
                    .take(CHUNK_SIZE)
                    .collect::<io::Result<_>>()?;

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

    let unmatched_records = record_pairs.unmatched_records();
    let unmatched_record_count = unmatched_records.len();

    if unmatched_record_count > 0 {
        warn!(unmatched_record_count, "found unmatched record");

        for record in unmatched_records {
            let event =
                count_single_end_record(interval_trees, filter, strand_specification, &record)?;

            ctx.add_event(event);
        }
    }

    Ok(ctx)
}

pub fn count_paired_end_record_pair<'f>(
    interval_trees: &IntervalTrees<'f>,
    filter: &'f Filter,
    strand_specification: StrandSpecification,
    r1: &bam::Record,
    r2: &bam::Record,
) -> io::Result<Event<'f>> {
    if let Some(event) = filter.filter_pair(r1, r2)? {
        return Ok(event);
    }

    let mut intersections = HashSet::new();

    let r1_is_reverse_complemented =
        resolve_is_reverse_complemented(r1.flags().is_reverse_complemented(), strand_specification);

    if let Some(event) = count_record(
        interval_trees,
        strand_specification,
        r1_is_reverse_complemented,
        r1,
        &mut intersections,
    )? {
        return Ok(event);
    }

    let r2_is_reverse_complemented = !resolve_is_reverse_complemented(
        r2.flags().is_reverse_complemented(),
        strand_specification,
    );

    if let Some(event) = count_record(
        interval_trees,
        strand_specification,
        r2_is_reverse_complemented,
        r2,
        &mut intersections,
    )? {
        return Ok(event);
    }

    Ok(resolve_intersections(&intersections))
}

fn count_record<'f>(
    interval_trees: &IntervalTrees<'f>,
    strand_specification: StrandSpecification,
    is_reverse_complemented: bool,
    record: &bam::Record,
    intersections: &mut HashSet<&'f str>,
) -> io::Result<Option<Event<'f>>> {
    let reference_sequence_id = record
        .reference_sequence_id()
        .transpose()?
        .expect("missing reference sequence ID");

    let Some(interval_tree) = interval_trees.get(reference_sequence_id) else {
        return Ok(Some(Event::Miss));
    };

    let cigar = record.cigar();
    let mut ops = cigar.iter();

    let alignment_start = record
        .alignment_start()
        .transpose()?
        .expect("missing alignment start");

    let intervals = MatchIntervals::new(&mut ops, alignment_start);

    intersect(
        intersections,
        interval_tree,
        intervals,
        strand_specification,
        is_reverse_complemented,
    )?;

    Ok(None)
}

fn intersect<'f>(
    intersections: &mut HashSet<&'f str>,
    interval_tree: &IntervalTree<Position, Entry<'f>>,
    intervals: MatchIntervals,
    strand_specification: StrandSpecification,
    is_reverse_complemented: bool,
) -> io::Result<()> {
    use noodles::gff::feature::record::Strand;

    for result in intervals {
        let interval = result?;

        for entry in interval_tree.find(interval.clone()) {
            let (name, strand) = entry.get();

            if strand_specification == StrandSpecification::None
                || (*strand == Strand::Reverse && is_reverse_complemented)
                || (*strand == Strand::Forward && !is_reverse_complemented)
            {
                intersections.insert(name);
            }
        }
    }

    Ok(())
}

fn resolve_intersections<'f>(intersections: &HashSet<&'f str>) -> Event<'f> {
    if intersections.is_empty() {
        Event::Miss
    } else if intersections.len() == 1 {
        // SAFETY: `intersections` is nonempty.
        let name = intersections.iter().next().unwrap();
        Event::Hit(name)
    } else {
        Event::Ambiguous
    }
}

fn resolve_is_reverse_complemented(
    is_reverse_complemented: bool,
    strand_specification: StrandSpecification,
) -> bool {
    if strand_specification == StrandSpecification::Reverse {
        !is_reverse_complemented
    } else {
        is_reverse_complemented
    }
}
