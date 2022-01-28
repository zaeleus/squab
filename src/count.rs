pub mod context;
mod filter;
mod reader;
mod writer;

pub use self::{context::Context, filter::Filter, reader::Reader, writer::Writer};

use std::{collections::HashSet, io, sync::Arc};

use futures::{StreamExt, TryFutureExt, TryStreamExt};
use interval_tree::IntervalTree;
use noodles::{
    bam, gff,
    sam::{self, header::ReferenceSequences},
};
use tokio::io::AsyncRead;

use crate::{Entry, Features, MatchIntervals, PairPosition, RecordPairs, StrandSpecification};

use self::context::Event;

const CHUNK_SIZE: usize = 8192;

pub async fn count_single_end_records<R>(
    mut reader: bam::AsyncReader<R>,
    features: Arc<Features>,
    reference_sequences: Arc<ReferenceSequences>,
    filter: Filter,
    strand_specification: StrandSpecification,
    worker_count: usize,
) -> io::Result<Context>
where
    R: AsyncRead + Unpin,
{
    use tokio::task::spawn_blocking;

    let ctx = reader
        .records()
        .try_chunks(CHUNK_SIZE)
        .map(move |result| {
            let features = features.clone();
            let reference_sequences = reference_sequences.clone();
            let filter = filter.clone();

            result
                .map(|records| {
                    spawn_blocking(move || {
                        let mut ctx = Context::default();

                        for record in records {
                            let event = count_single_end_record(
                                &features,
                                &reference_sequences,
                                &filter,
                                strand_specification,
                                &record,
                            )?;

                            ctx.add_event(event);
                        }

                        Ok(ctx)
                    })
                    .map_err(io::Error::from)
                })
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
        .try_buffer_unordered(worker_count)
        .try_fold(Context::default(), |mut ctx, result| async move {
            result.map(|c| {
                ctx.add(&c);
                ctx
            })
        })
        .await?;

    Ok(ctx)
}

pub fn count_single_end_record(
    features: &Features,
    reference_sequences: &ReferenceSequences,
    filter: &Filter,
    strand_specification: StrandSpecification,
    record: &bam::Record,
) -> io::Result<Event> {
    if let Some(event) = filter.filter(record)? {
        return Ok(event);
    }

    let cigar = record.cigar();
    let start = record
        .position()
        .map(i32::from)
        .map(|n| n as u64)
        .expect("record cannot be unmapped");
    let flags = record.flags();

    let is_reverse = match strand_specification {
        StrandSpecification::Reverse => !flags.is_reverse_complemented(),
        _ => flags.is_reverse_complemented(),
    };

    let tree = match get_tree(
        features,
        reference_sequences,
        record.reference_sequence_id(),
    )? {
        Some(t) => t,
        None => return Ok(Event::NoFeature),
    };

    let intervals = MatchIntervals::new(cigar, start);

    let set = find(tree, intervals, strand_specification, is_reverse)?;

    Ok(update_intersections(set))
}

pub async fn count_paired_end_records<R>(
    reader: bam::AsyncReader<R>,
    features: Arc<Features>,
    reference_sequences: Arc<ReferenceSequences>,
    filter: Filter,
    strand_specification: StrandSpecification,
    worker_count: usize,
) -> io::Result<Context>
where
    R: AsyncRead + Unpin,
{
    use futures::stream;
    use tokio::task::spawn_blocking;

    let primary_only = !filter.with_secondary_records() && !filter.with_supplementary_records();
    let mut record_pairs = RecordPairs::new(reader, primary_only);

    let features_clone = features.clone();
    let reference_sequences_clone = reference_sequences.clone();
    let filter_clone = filter.clone();

    let mut ctx = stream::unfold(&mut record_pairs, |pairs| async {
        match pairs.next_pair().await {
            Ok(Some(pair)) => Some((Ok(pair), pairs)),
            Ok(None) => None,
            Err(e) => Some((Err(e), pairs)),
        }
    })
    .try_chunks(CHUNK_SIZE)
    .map(move |result| {
        let features = features.clone();
        let reference_sequences = reference_sequences.clone();
        let filter = filter.clone();

        result
            .map(|pairs| {
                spawn_blocking(move || {
                    let mut ctx = Context::default();

                    for (r1, r2) in pairs {
                        let event = count_paired_end_record_pair(
                            &features,
                            &reference_sequences,
                            &filter,
                            strand_specification,
                            &r1,
                            &r2,
                        )?;

                        ctx.add_event(event);
                    }

                    Ok(ctx)
                })
                .map_err(io::Error::from)
            })
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
    .try_buffer_unordered(worker_count)
    .try_fold(Context::default(), |mut ctx, result| async move {
        result.map(|c| {
            ctx.add(&c);
            ctx
        })
    })
    .await?;

    for record in record_pairs.singletons() {
        let event = count_paired_end_singleton_record(
            &features_clone,
            &reference_sequences_clone,
            &filter_clone,
            strand_specification,
            &record,
        )?;

        ctx.add_event(event);
    }

    Ok(ctx)
}

pub fn count_paired_end_record_pair(
    features: &Features,
    reference_sequences: &ReferenceSequences,
    filter: &Filter,
    strand_specification: StrandSpecification,
    r1: &bam::Record,
    r2: &bam::Record,
) -> io::Result<Event> {
    if let Some(event) = filter.filter_pair(r1, r2)? {
        return Ok(event);
    }

    let cigar = r1.cigar();
    let start = r1
        .position()
        .map(i32::from)
        .map(|n| n as u64)
        .expect("record cannot be unmapped");
    let f1 = r1.flags();

    let is_reverse = match strand_specification {
        StrandSpecification::Reverse => !f1.is_reverse_complemented(),
        _ => f1.is_reverse_complemented(),
    };

    let tree = match get_tree(features, reference_sequences, r1.reference_sequence_id())? {
        Some(t) => t,
        None => return Ok(Event::NoFeature),
    };

    let intervals = MatchIntervals::new(cigar, start);

    let mut set = find(tree, intervals, strand_specification, is_reverse)?;

    let cigar = r2.cigar();
    let start = r2
        .position()
        .map(i32::from)
        .map(|n| n as u64)
        .expect("record cannot be unmapped");
    let f2 = r2.flags();

    let is_reverse = match strand_specification {
        StrandSpecification::Reverse => f2.is_reverse_complemented(),
        _ => !f2.is_reverse_complemented(),
    };

    let tree = match get_tree(features, reference_sequences, r2.reference_sequence_id())? {
        Some(t) => t,
        None => return Ok(Event::NoFeature),
    };

    let intervals = MatchIntervals::new(cigar, start);

    let set2 = find(tree, intervals, strand_specification, is_reverse)?;

    set.extend(set2.into_iter());

    Ok(update_intersections(set))
}

pub fn count_paired_end_singleton_record(
    features: &Features,
    reference_sequences: &ReferenceSequences,
    filter: &Filter,
    strand_specification: StrandSpecification,
    record: &bam::Record,
) -> io::Result<Event> {
    if let Some(event) = filter.filter(record)? {
        return Ok(event);
    }

    let cigar = record.cigar();
    let start = record
        .position()
        .map(i32::from)
        .map(|n| n as u64)
        .expect("record cannot be unmapped");
    let flags = record.flags();

    let is_reverse = match PairPosition::try_from(record) {
        Ok(PairPosition::First) => match strand_specification {
            StrandSpecification::Reverse => !flags.is_reverse_complemented(),
            _ => flags.is_reverse_complemented(),
        },
        Ok(PairPosition::Second) => match strand_specification {
            StrandSpecification::Reverse => flags.is_reverse_complemented(),
            _ => !flags.is_reverse_complemented(),
        },
        Err(e) => return Err(io::Error::new(io::ErrorKind::InvalidData, e)),
    };

    let tree = match get_tree(
        features,
        reference_sequences,
        record.reference_sequence_id(),
    )? {
        Some(t) => t,
        None => return Ok(Event::NoFeature),
    };

    let intervals = MatchIntervals::new(cigar, start);

    let set = find(tree, intervals, strand_specification, is_reverse)?;

    Ok(update_intersections(set))
}

fn find(
    tree: &IntervalTree<u64, Entry>,
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

fn get_reference_sequence(
    reference_sequences: &ReferenceSequences,
    reference_sequence_id: Option<bam::record::ReferenceSequenceId>,
) -> io::Result<&sam::header::ReferenceSequence> {
    reference_sequence_id
        .and_then(|id| {
            reference_sequences
                .get_index(i32::from(id) as usize)
                .map(|(_, rs)| rs)
        })
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid reference sequence ID: {:?}", reference_sequence_id),
            )
        })
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

pub fn get_tree<'t>(
    features: &'t Features,
    reference_sequences: &ReferenceSequences,
    reference_sequence_id: Option<bam::record::ReferenceSequenceId>,
) -> io::Result<Option<&'t IntervalTree<u64, Entry>>> {
    let reference_sequence = get_reference_sequence(reference_sequences, reference_sequence_id)?;
    let reference_sequence_name = reference_sequence.name();
    Ok(features.get(reference_sequence_name.as_str()))
}

#[cfg(test)]
mod tests {
    use sam::header::ReferenceSequence;

    use super::*;

    fn build_reference_sequences() -> Result<ReferenceSequences, Box<dyn std::error::Error>> {
        let reference_sequences = [
            ("sq0".parse()?, 8),
            ("sq1".parse()?, 13),
            ("sq2".parse()?, 21),
        ]
        .into_iter()
        .map(
            |(name, len): (sam::header::reference_sequence::Name, i32)| {
                let sn = name.to_string();
                ReferenceSequence::new(name, len).map(|rs| (sn, rs))
            },
        )
        .collect::<Result<_, _>>()?;

        Ok(reference_sequences)
    }

    #[test]
    fn test_get_reference() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequences = build_reference_sequences()?;

        let reference_sequence_id = Some(bam::record::ReferenceSequenceId::try_from(1)?);
        let reference_sequence =
            get_reference_sequence(&reference_sequences, reference_sequence_id)?;
        assert_eq!(reference_sequence.name().as_str(), "sq1");
        assert_eq!(reference_sequence.len(), 13);

        let reference_sequence_id = None;
        let reference_sequence =
            get_reference_sequence(&reference_sequences, reference_sequence_id);
        assert!(reference_sequence.is_err());

        let reference_sequence_id = Some(bam::record::ReferenceSequenceId::try_from(5)?);
        let reference_sequence =
            get_reference_sequence(&reference_sequences, reference_sequence_id);
        assert!(reference_sequence.is_err());

        Ok(())
    }
}
