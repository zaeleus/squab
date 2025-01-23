#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Event<'f> {
    Hit(&'f str),
    Miss,
    Ambiguous,
    LowQuality,
    Unmapped,
    Nonunique,
    Skip,
}
