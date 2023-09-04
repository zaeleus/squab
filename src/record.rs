use std::io;

use noodles::{bam::lazy::record::Cigar, core::Position};

pub fn alignment_end(
    alignment_start: Option<Position>,
    cigar: &Cigar<'_>,
) -> Option<io::Result<Position>> {
    let start = alignment_start?;

    let span = match alignment_span(cigar) {
        Ok(n) => n,
        Err(e) => return Some(Err(e)),
    };

    let end = usize::from(start) + span - 1;

    Position::new(end).map(Ok)
}

fn alignment_span(cigar: &Cigar<'_>) -> io::Result<usize> {
    let mut span = 0;

    for result in cigar.iter() {
        let op = result?;

        if op.kind().consumes_reference() {
            span += op.len();
        }
    }

    Ok(span)
}
