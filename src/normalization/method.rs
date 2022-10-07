use std::{error, fmt, str::FromStr};

/// Normalization method
#[derive(Clone, Copy, Debug, Eq, PartialEq, clap::ValueEnum)]
pub enum Method {
    /// fragments per kilobase per million mapped reads
    Fpkm,
    /// transcripts per million
    Tpm,
}

#[derive(Debug, Eq, PartialEq)]
pub struct ParseError(String);

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid method: {}", self.0)
    }
}

impl error::Error for ParseError {}

impl FromStr for Method {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "fpkm" => Ok(Self::Fpkm),
            "tpm" => Ok(Self::Tpm),
            _ => Err(ParseError(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("fpkm".parse::<Method>()?, Method::Fpkm);
        assert_eq!("tpm".parse::<Method>()?, Method::Tpm);

        assert!("".parse::<Method>().is_err());
        assert!("squab".parse::<Method>().is_err());
        assert!("FPKM".parse::<Method>().is_err());

        Ok(())
    }
}
