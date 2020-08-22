mod normalize;
mod quantify;

pub use self::{normalize::normalize, quantify::quantify};

use std::str::FromStr;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum StrandSpecificationOption {
    None,
    Forward,
    Reverse,
    Auto,
}

impl FromStr for StrandSpecificationOption {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "none" => Ok(Self::None),
            "forward" => Ok(Self::Forward),
            "reverse" => Ok(Self::Reverse),
            "auto" => Ok(Self::Auto),
            _ => Err(()),
        }
    }
}
