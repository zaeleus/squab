#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Event {
    NoFeature,
    Ambiguous,
    LowQuality,
    Unmapped,
    Nonunique,
}
