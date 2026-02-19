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
