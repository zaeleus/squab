mod normalize;
mod quantify;

pub use self::{normalize::normalize, quantify::quantify};

use std::str::FromStr;

#[derive(Clone, Copy, Debug, Eq, PartialEq, clap::ValueEnum)]
pub enum StrandSpecificationOption {
    None,
    Forward,
    Reverse,
    Auto,
}

impl FromStr for StrandSpecificationOption {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "none" => Ok(Self::None),
            "forward" => Ok(Self::Forward),
            "reverse" => Ok(Self::Reverse),
            "auto" => Ok(Self::Auto),
            _ => Err(format!("invalid strand specification option: {}", s)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!("none".parse(), Ok(StrandSpecificationOption::None));
        assert_eq!("forward".parse(), Ok(StrandSpecificationOption::Forward));
        assert_eq!("reverse".parse(), Ok(StrandSpecificationOption::Reverse));
        assert_eq!("auto".parse(), Ok(StrandSpecificationOption::Auto));

        assert!("".parse::<StrandSpecificationOption>().is_err());
        assert!("noodles".parse::<StrandSpecificationOption>().is_err());
    }
}
