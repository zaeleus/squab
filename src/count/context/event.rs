#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Event {
    Hit(String),
    Miss,
    Ambiguous,
    LowQuality,
    Unmapped,
    Nonunique,
    Skip,
}
