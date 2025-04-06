use bstr::BStr;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Event<'f> {
    Hit(&'f BStr),
    Miss,
    Ambiguous,
    LowQuality,
    Unmapped,
    Nonunique,
    Skip,
}
