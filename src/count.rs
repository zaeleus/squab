mod context;
mod filter;
mod reader;
mod writer;

pub use self::{context::Context, filter::Filter, reader::Reader, writer::Writer};

use std::{collections::HashSet, convert::TryFrom, io};

use interval_tree::IntervalTree;
use noodles::{
    bam, gff,
    sam::{self, header::ReferenceSequences},
};

use crate::{Entry, Features, MatchIntervals, PairPosition, RecordPairs, StrandSpecification};

use self::context::Event;

pub fn count_single_end_records<I>(
    records: I,
    features: &Features,
    references: &ReferenceSequences,
    filter: &Filter,
    strand_specification: StrandSpecification,
) -> io::Result<Context>
where
    I: Iterator<Item = io::Result<bam::Record>>,
{
    let mut ctx = Context::default();

    for result in records {
        let record = result?;

        count_single_end_record(
            &mut ctx,
            features,
            references,
            filter,
            strand_specification,
            &record,
        )?;
    }

    Ok(ctx)
}

pub fn count_single_end_record(
    ctx: &mut Context,
    features: &Features,
    reference_sequences: &ReferenceSequences,
    filter: &Filter,
    strand_specification: StrandSpecification,
    record: &bam::Record,
) -> io::Result<()> {
    if let Some(event) = filter.filter(record)? {
        ctx.add_event(event);
        return Ok(());
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

    let intervals = MatchIntervals::new(&cigar, start);

    let tree = match get_tree(
        features,
        reference_sequences,
        record.reference_sequence_id(),
    )? {
        Some(t) => t,
        None => {
            ctx.add_event(Event::NoFeature);
            return Ok(());
        }
    };

    let set = find(tree, intervals, strand_specification, is_reverse)?;
    let event = update_intersections(set);

    ctx.add_event(event);

    Ok(())
}

pub fn count_paired_end_records<I>(
    records: I,
    features: &Features,
    reference_sequences: &ReferenceSequences,
    filter: &Filter,
    strand_specification: StrandSpecification,
) -> io::Result<(Context, RecordPairs<I>)>
where
    I: Iterator<Item = io::Result<bam::Record>>,
{
    let mut ctx = Context::default();

    let primary_only = !filter.with_secondary_records() && !filter.with_supplementary_records();
    let mut pairs = RecordPairs::new(records, primary_only);

    for pair in &mut pairs {
        let (r1, r2) = pair?;

        if let Some(event) = filter.filter_pair(&r1, &r2)? {
            ctx.add_event(event);
            continue;
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

        let intervals = MatchIntervals::new(&cigar, start);

        let tree = match get_tree(features, reference_sequences, r1.reference_sequence_id())? {
            Some(t) => t,
            None => {
                ctx.add_event(Event::NoFeature);
                continue;
            }
        };

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

        let intervals = MatchIntervals::new(&cigar, start);

        let tree = match get_tree(features, reference_sequences, r2.reference_sequence_id())? {
            Some(t) => t,
            None => {
                ctx.add_event(Event::NoFeature);
                continue;
            }
        };

        let set2 = find(tree, intervals, strand_specification, is_reverse)?;

        set.extend(set2.into_iter());

        let event = update_intersections(set);

        ctx.add_event(event);
    }

    Ok((ctx, pairs))
}

pub fn count_paired_end_record_singletons<I>(
    records: I,
    features: &Features,
    reference_sequences: &ReferenceSequences,
    filter: &Filter,
    strand_specification: StrandSpecification,
) -> io::Result<Context>
where
    I: Iterator<Item = io::Result<bam::Record>>,
{
    let mut ctx = Context::default();

    for result in records {
        let record = result?;

        if let Some(event) = filter.filter(&record)? {
            ctx.add_event(event);
            continue;
        }

        let cigar = record.cigar();
        let start = record
            .position()
            .map(i32::from)
            .map(|n| n as u64)
            .expect("record cannot be unmapped");
        let flags = record.flags();

        let is_reverse = match PairPosition::try_from(&record) {
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

        let intervals = MatchIntervals::new(&cigar, start);

        let tree = match get_tree(
            features,
            reference_sequences,
            record.reference_sequence_id(),
        )? {
            Some(t) => t,
            None => {
                ctx.add_event(Event::NoFeature);
                continue;
            }
        };

        let set = find(tree, intervals, strand_specification, is_reverse)?;
        let event = update_intersections(set);

        ctx.add_event(event);
    }

    Ok(ctx)
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
    let name = reference_sequence.name();
    Ok(features.get(name))
}

#[cfg(test)]
mod tests {
    use sam::header::ReferenceSequence;

    use super::*;

    fn build_reference_sequences(
    ) -> Result<ReferenceSequences, sam::header::reference_sequence::NewError> {
        vec![("sq0", 8), ("sq1", 13), ("sq2", 21)]
            .into_iter()
            .map(|(name, len)| ReferenceSequence::new(name, len).map(|rs| (name.into(), rs)))
            .collect()
    }

    #[test]
    fn test_get_reference() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequences = build_reference_sequences()?;

        let reference_sequence_id = Some(bam::record::ReferenceSequenceId::try_from(1)?);
        let reference_sequence =
            get_reference_sequence(&reference_sequences, reference_sequence_id)?;
        assert_eq!(reference_sequence.name(), "sq1");
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
