use bstr::BStr;

#[derive(Default)]
pub(super) enum Intersections<'f> {
    #[default]
    Empty,
    One(&'f BStr),
    Many,
}

impl<'f> Intersections<'f> {
    pub(super) fn insert(&mut self, name: &'f BStr) {
        match self {
            Self::Empty => *self = Self::One(name),
            Self::One(other) if name != *other => *self = Self::Many,
            _ => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;

    use super::*;

    #[test]
    fn test_default() {
        assert!(matches!(Intersections::default(), Intersections::Empty));
    }

    #[test]
    fn test_insert() {
        let name_0 = b"f0".as_bstr();
        let name_1 = b"f1".as_bstr();

        let mut intersections = Intersections::default();

        intersections.insert(name_0);
        assert!(matches!(intersections, Intersections::One(s) if s == name_0));

        intersections.insert(name_0);
        assert!(matches!(intersections, Intersections::One(s) if s == name_0));

        intersections.insert(name_1);
        assert!(matches!(intersections, Intersections::Many));
    }
}
