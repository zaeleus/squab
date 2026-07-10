pub mod context;
mod filter;
mod intersections;
mod try_buffered_chunks;

use std::{
    io::{self, Read},
    num::NonZero,
};

use noodles::{bam, core::Position};
use rayon::iter::{ParallelBridge, ParallelIterator};

pub use self::{context::Context, filter::Filter};
use self::{context::Event, intersections::Intersections, try_buffered_chunks::TryBufferedChunks};
use crate::{
    Entry, IntervalTrees, MatchIntervals, RecordPairs, StrandSpecification,
    collections::IntervalTree,
};

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
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(worker_count.get())
        .build()
        .unwrap();

    pool.install(|| {
        TryBufferedChunks::new(reader.records(), CHUNK_SIZE)
            .par_bridge()
            .try_fold(Context::default, |mut ctx, result| {
                let chunk = result?;

                for record in chunk {
                    let event = count_single_end_record(
                        interval_trees,
                        filter,
                        strand_specification,
                        &record,
                    )?;

                    ctx.add_event(event);
                }

                Ok(ctx)
            })
            .try_reduce(Context::default, |mut ctx, c| {
                ctx.add(&c);
                Ok(ctx)
            })
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

    let mut intersections = Intersections::default();

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
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(worker_count.get())
        .build()
        .unwrap();

    let primary_only = !filter.with_secondary_records() && !filter.with_supplementary_records();
    let mut record_pairs = RecordPairs::new(reader, primary_only);

    let mut ctx = pool.install(|| {
        TryBufferedChunks::new(&mut record_pairs, CHUNK_SIZE)
            .par_bridge()
            .try_fold(Context::default, |mut ctx, result| {
                let chunk = result?;

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

                Ok::<_, io::Error>(ctx)
            })
            .try_reduce(Context::default, |mut ctx, c| {
                ctx.add(&c);
                Ok(ctx)
            })
    })?;

    for record in record_pairs.unmatched_records() {
        let event = count_single_end_record(interval_trees, filter, strand_specification, &record)?;
        ctx.add_event(event);
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

    let mut intersections = Intersections::default();

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
    intersections: &mut Intersections<'f>,
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
    intersections: &mut Intersections<'f>,
    interval_tree: &IntervalTree<Position, Entry<'f>>,
    intervals: MatchIntervals,
    strand_specification: StrandSpecification,
    is_reverse_complemented: bool,
) -> io::Result<()> {
    use noodles::gff::feature::record::Strand;

    for result in intervals {
        let interval = result?;

        for (_, (name, strand)) in interval_tree.find(interval.clone()) {
            if strand_specification == StrandSpecification::None
                || (*strand == Strand::Reverse && is_reverse_complemented)
                || (*strand == Strand::Forward && !is_reverse_complemented)
            {
                intersections.insert(name);

                if matches!(intersections, Intersections::Many) {
                    return Ok(());
                }
            }
        }
    }

    Ok(())
}

fn resolve_intersections<'f>(intersections: &Intersections<'f>) -> Event<'f> {
    match intersections {
        Intersections::Empty => Event::Miss,
        Intersections::One(name) => Event::Hit(name),
        Intersections::Many => Event::Ambiguous,
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
