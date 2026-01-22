mod event;

use bstr::BStr;
use rapidhash::RapidHashMap;

pub use self::event::Event;

pub type Counts<'f> = RapidHashMap<&'f BStr, u64>;

#[derive(Default)]
pub struct Context<'f> {
    pub counts: Counts<'f>,
    pub miss: u64,
    pub ambiguous: u64,
    pub low_quality: u64,
    pub unmapped: u64,
    pub nonunique: u64,
}

impl<'f> Context<'f> {
    pub fn add(&mut self, other: &Context<'f>) {
        for (name, count) in other.counts.iter() {
            let entry = self.counts.entry(name).or_insert(0);
            *entry += count;
        }

        self.miss += other.miss;
        self.ambiguous += other.ambiguous;
        self.low_quality += other.low_quality;
        self.unmapped += other.unmapped;
        self.nonunique += other.nonunique;
    }

    pub fn add_event(&mut self, event: Event<'f>) {
        match event {
            Event::Hit(name) => {
                let count = self.counts.entry(name).or_insert(0);
                *count += 1;
            }
            Event::Miss => self.miss += 1,
            Event::Ambiguous => self.ambiguous += 1,
            Event::LowQuality => self.low_quality += 1,
            Event::Unmapped => self.unmapped += 1,
            Event::Nonunique => self.nonunique += 1,
            Event::Skip => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        let mut ctx_a = Context::default();

        ctx_a.counts.insert(BStr::new("AADAT"), 2);
        ctx_a.miss = 3;
        ctx_a.ambiguous = 5;
        ctx_a.low_quality = 8;
        ctx_a.unmapped = 13;
        ctx_a.nonunique = 21;

        let mut ctx_b = Context::default();

        ctx_b.counts.insert(BStr::new("AADAT"), 2);
        ctx_b.counts.insert(BStr::new("CLN3"), 3);
        ctx_b.miss = 5;
        ctx_b.ambiguous = 8;
        ctx_b.low_quality = 13;
        ctx_b.unmapped = 21;
        ctx_b.nonunique = 34;

        ctx_a.add(&ctx_b);

        assert_eq!(ctx_a.counts.len(), 2);
        assert_eq!(ctx_a.counts[BStr::new("AADAT")], 4);
        assert_eq!(ctx_a.counts[BStr::new("CLN3")], 3);

        assert_eq!(ctx_a.miss, 8);
        assert_eq!(ctx_a.ambiguous, 13);
        assert_eq!(ctx_a.low_quality, 21);
        assert_eq!(ctx_a.unmapped, 34);
        assert_eq!(ctx_a.nonunique, 55);
    }

    #[test]
    fn test_add_event() {
        let mut ctx = Context::default();
        ctx.add_event(Event::Hit(BStr::new("AADAT")));
        ctx.add_event(Event::Miss);
        ctx.add_event(Event::Ambiguous);
        ctx.add_event(Event::LowQuality);
        ctx.add_event(Event::Unmapped);
        ctx.add_event(Event::Nonunique);

        assert_eq!(ctx.counts.len(), 1);
        assert_eq!(ctx.counts[BStr::new("AADAT")], 1);

        assert_eq!(ctx.miss, 1);
        assert_eq!(ctx.ambiguous, 1);
        assert_eq!(ctx.low_quality, 1);
        assert_eq!(ctx.unmapped, 1);
        assert_eq!(ctx.nonunique, 1);
    }
}
