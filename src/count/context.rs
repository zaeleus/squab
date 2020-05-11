use std::collections::HashMap;

#[derive(Default)]
pub struct Context {
    pub counts: HashMap<String, u64>,
    pub no_feature: u64,
    pub ambiguous: u64,
    pub low_quality: u64,
    pub unmapped: u64,
    pub nonunique: u64,
}

impl Context {
    pub fn add(&mut self, other: &Context) {
        for (name, count) in other.counts.iter() {
            let entry = self.counts.entry(name.to_string()).or_insert(0);
            *entry += count;
        }

        self.no_feature += other.no_feature;
        self.ambiguous += other.ambiguous;
        self.low_quality += other.low_quality;
        self.unmapped += other.unmapped;
        self.nonunique += other.nonunique;
    }
}
