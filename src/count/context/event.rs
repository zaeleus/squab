#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Event {
    Hit(String),
    NoFeature,
    Ambiguous,
    LowQuality,
    Unmapped,
    Nonunique,
    Skip,
}
