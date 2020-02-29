use std::str::FromStr;

#[derive(Clone, Copy, Debug)]
pub enum QuantificationMethod {
    Count,
    Fpkm,
    Tpm,
}

impl FromStr for QuantificationMethod {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "count" => Ok(Self::Count),
            "fpkm" => Ok(Self::Fpkm),
            "tpm" => Ok(Self::Tpm),
            _ => Err(()),
        }
    }
}
