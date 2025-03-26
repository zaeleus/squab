use std::{iter::FusedIterator, ops::RangeInclusive, slice};

const MIN_LEVEL: usize = 3;

/// An augmented interval tree.
///
/// This is a port of [`lh3/cgranges/cpp/IITree.h`][IITree.h].
///
/// [IITree.h]: https://github.com/lh3/cgranges/blob/b3d5e2c5b9a0a379f54592ab85f6cff5d58c387e/cpp/IITree.h
#[derive(Debug, Default)]
pub struct IntervalTree<K, V> {
    nodes: Vec<Node<K, V>>,
    height: usize,
}

impl<K, V> IntervalTree<K, V>
where
    K: Copy + Ord,
{
    pub fn find(&self, key: RangeInclusive<K>) -> Find<'_, K, V> {
        Find::new(&self.nodes, self.height, key)
    }

    fn index(&mut self) {
        let nodes = &mut self.nodes;
        let len = nodes.len();

        if len == 0 {
            return;
        }

        nodes.sort_unstable_by_key(|i| i.start());

        let mut last_id = if len % 2 == 0 { len - 2 } else { len - 1 };
        let mut last_max = nodes[last_id].max;

        let mut level = 1;
        let mut step = 1 << level;

        while step <= len {
            let first_id = first_id_of_level(level);

            for id in (first_id..len).step_by(step) {
                let l = nodes[left_child_id(id, level)].max;

                let r = nodes
                    .get(right_child_id(id, level))
                    .map(|node| node.max)
                    .unwrap_or(last_max);

                nodes[id].max = nodes[id].end().max(l).max(r);
            }

            last_id = parent_id(last_id, level - 1);

            if let Some(node) = nodes.get(last_id) {
                last_max = last_max.max(node.max);
            }

            level += 1;
            step = 1 << level;
        }

        self.height = level - 1;
    }
}

impl<K, V> FromIterator<(RangeInclusive<K>, V)> for IntervalTree<K, V>
where
    K: Copy + Ord,
{
    fn from_iter<T: IntoIterator<Item = (RangeInclusive<K>, V)>>(iter: T) -> Self {
        let nodes = iter
            .into_iter()
            .map(|(key, value)| Node::new(key, value))
            .collect();

        let mut granges = Self { nodes, height: 0 };
        granges.index();
        granges
    }
}

#[derive(Debug)]
struct Node<K, V> {
    key: RangeInclusive<K>,
    value: V,
    max: K,
}

impl<K, V> Node<K, V>
where
    K: Copy + Ord,
{
    fn new(key: RangeInclusive<K>, value: V) -> Self {
        let max = *key.end();
        Self { key, value, max }
    }

    fn start(&self) -> K {
        *self.key.start()
    }

    fn end(&self) -> K {
        *self.key.end()
    }
}

enum State<'t, K, V> {
    Traverse,
    Scan(slice::Iter<'t, Node<K, V>>),
}

pub struct Find<'t, K, V> {
    nodes: &'t [Node<K, V>],
    key: RangeInclusive<K>,
    stack: Vec<(usize, usize)>,
    state: State<'t, K, V>,
}

impl<'t, K, V> Find<'t, K, V> {
    fn new(nodes: &'t [Node<K, V>], height: usize, key: RangeInclusive<K>) -> Self {
        let mut stack = Vec::with_capacity(1 << 4);
        stack.push((root_id(height), height));

        Self {
            nodes,
            key,
            stack,
            state: State::Traverse,
        }
    }
}

impl<'t, K, V> Iterator for Find<'t, K, V>
where
    K: Copy + Ord,
{
    type Item = (&'t RangeInclusive<K>, &'t V);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.state {
                State::Traverse => {
                    let (id, level) = self.stack.pop()?;

                    if level <= MIN_LEVEL {
                        let start_id = id >> level << level;
                        let node_count = perfect_binary_tree_node_count(level);
                        let end_id = (start_id + node_count).min(self.nodes.len());
                        self.state = State::Scan(self.nodes[start_id..end_id].iter());
                    } else {
                        let i = left_child_id(id, level);

                        if self
                            .nodes
                            .get(i)
                            .map(|node| node.max >= *self.key.start())
                            .unwrap_or(true)
                        {
                            self.stack.push((i, level - 1));
                        }

                        if let Some(node) = self.nodes.get(id) {
                            if node.start() <= *self.key.end() {
                                let i = right_child_id(id, level);
                                self.stack.push((i, level - 1));

                                if *self.key.start() <= node.end() {
                                    return Some((&node.key, &node.value));
                                }
                            }
                        }
                    }
                }
                State::Scan(ref mut iter) => match iter.next() {
                    Some(node) if node.start() <= *self.key.end() => {
                        if *self.key.start() <= node.end() {
                            return Some((&node.key, &node.value));
                        }
                    }
                    _ => self.state = State::Traverse,
                },
            }
        }
    }
}

impl<K, V> FusedIterator for Find<'_, K, V> where K: Copy + Ord {}

fn perfect_binary_tree_node_count(height: usize) -> usize {
    (1 << (height + 1)) - 1
}

fn root_id(height: usize) -> usize {
    first_id_of_level(height)
}

fn first_id_of_level(level: usize) -> usize {
    (1 << level) - 1
}

fn parent_id(id: usize, level: usize) -> usize {
    if (id >> (level + 1)) % 2 == 0 {
        id + (1 << level)
    } else {
        id - (1 << level)
    }
}

fn left_child_id(id: usize, level: usize) -> usize {
    id - (1 << (level - 1))
}

fn right_child_id(id: usize, level: usize) -> usize {
    id + (1 << (level - 1))
}

#[cfg(test)]
mod tests {
    use super::*;

    // 000000000111111111122222222223333333333444444444
    // 123456789012345678901234567890123456789012345678
    // ─────────┼─────────┼─────────┼─────────┼────────
    //  aaaa eeee g  iilll mmmm pp  qqqqqqqqqq    x
    //   bbbbbbf  hhhhkkkkkkk oooo  rr u  vv         yy
    //    ccccc             nnnnnn   ssssswwww        z
    //     dddd                      ttttt
    const NODES: [(RangeInclusive<usize>, char); 26] = [
        (2..=5, 'a'),
        (3..=8, 'b'),
        (4..=8, 'c'),
        (5..=8, 'd'),
        (7..=10, 'e'),
        (9..=9, 'f'),
        (12..=12, 'g'),
        (12..=15, 'h'),
        (15..=16, 'i'),
        (15..=18, 'j'),
        (16..=22, 'k'),
        (17..=19, 'l'),
        (21..=24, 'm'),
        (22..=27, 'n'),
        (24..=27, 'o'),
        (26..=27, 'p'),
        (30..=39, 'q'),
        (30..=31, 'r'),
        (31..=35, 's'),
        (31..=35, 't'),
        (33..=33, 'u'),
        (36..=37, 'v'),
        (36..=39, 'w'),
        (44..=44, 'x'),
        (47..=48, 'y'),
        (48..=48, 'z'),
    ];

    // level
    //   4 ┬                                                      15
    //     │                                                      2-5
    //     │                                  /                                        \
    //   3 ┼                           7                                                     23
    //     │                         12-15                                                  44-44
    //     │                 /                   \                                  /                   \
    //   2 ┼             3                          11                         19                           27
    //     │            5-8                        17-19                      31-35                         ()
    //     │       /           \               /           \              /           \               /           \
    //   1 ┼      1             5             9            13            17            21            25            29
    //     │     3-8           9-9          15-18         22-27         30-31         36-37         48-48          ()
    //     │   /     \       /     \       /     \       /     \       /     \       /     \       /     \       /     \
    //   0 ┴  0      2      4      6      8     10     12     14     16     18     20     22     24     26     28     30
    //       2-5    4-8    7-10  12-12  15-16  16-22  21-24  24-27  30-39  31-35  33-33  36-39  47-48   ()     ()     ()

    #[test]
    fn test_index() {
        let mut interval_tree = IntervalTree {
            nodes: NODES
                .into_iter()
                .map(|(key, value)| Node::new(key, value))
                .collect(),
            height: 0,
        };

        interval_tree.index();

        let actual: Vec<_> = interval_tree.nodes.iter().map(|node| node.max).collect();
        let expected = [
            5, 8, 8, 12, 10, 12, 12, 27, 16, 22, 22, 27, 24, 27, 27, 48, 39, 39, 35, 39, 33, 39,
            39, 48, 48, 48,
        ];

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_find() {
        fn t(
            interval_tree: &IntervalTree<usize, char>,
            key: RangeInclusive<usize>,
            expected: &[(RangeInclusive<usize>, char)],
        ) {
            assert_eq!(
                interval_tree.find(key).collect::<Vec<_>>(),
                expected
                    .iter()
                    .map(|(key, value)| (key, value))
                    .collect::<Vec<_>>()
            );
        }

        let interval_tree: IntervalTree<_, _> = NODES.into_iter().collect();

        t(&interval_tree, 2..=3, &[(2..=5, 'a'), (3..=8, 'b')]);
        t(
            &interval_tree,
            22..=25,
            &[
                (16..=22, 'k'),
                (21..=24, 'm'),
                (22..=27, 'n'),
                (24..=27, 'o'),
            ],
        );
        t(&interval_tree, 40..=43, &[]);
        t(&interval_tree, 44..=44, &[(44..=44, 'x')]);
        t(&interval_tree, 47..=48, &[(47..=48, 'y'), (48..=48, 'z')]);
        t(&interval_tree, 50..=56, &[]);
    }

    #[test]
    fn test_perfect_binary_tree_node_count() {
        assert_eq!(perfect_binary_tree_node_count(0), 1);
        assert_eq!(perfect_binary_tree_node_count(1), 3);
        assert_eq!(perfect_binary_tree_node_count(2), 7);
        assert_eq!(perfect_binary_tree_node_count(3), 15);
    }

    #[test]
    fn test_root_id() {
        assert_eq!(root_id(0), 0);
        assert_eq!(root_id(1), 1);
        assert_eq!(root_id(2), 3);
        assert_eq!(root_id(3), 7);
    }

    #[test]
    fn test_first_id_of_level() {
        assert_eq!(first_id_of_level(0), 0);
        assert_eq!(first_id_of_level(1), 1);
        assert_eq!(first_id_of_level(2), 3);
        assert_eq!(first_id_of_level(3), 7);
    }

    #[test]
    fn test_parent_id() {
        assert_eq!(parent_id(0, 0), 1);
        assert_eq!(parent_id(2, 0), 1);
        assert_eq!(parent_id(4, 0), 5);
        assert_eq!(parent_id(6, 0), 5);
        assert_eq!(parent_id(8, 0), 9);
        assert_eq!(parent_id(10, 0), 9);
        assert_eq!(parent_id(12, 0), 13);
        assert_eq!(parent_id(14, 0), 13);

        assert_eq!(parent_id(1, 1), 3);
        assert_eq!(parent_id(5, 1), 3);
        assert_eq!(parent_id(9, 1), 11);
        assert_eq!(parent_id(13, 1), 11);

        assert_eq!(parent_id(3, 2), 7);
        assert_eq!(parent_id(11, 2), 7);
    }

    #[test]
    fn test_left_child_id() {
        assert_eq!(left_child_id(7, 3), 3);

        assert_eq!(left_child_id(3, 2), 1);
        assert_eq!(left_child_id(11, 2), 9);

        assert_eq!(left_child_id(1, 1), 0);
        assert_eq!(left_child_id(5, 1), 4);
        assert_eq!(left_child_id(9, 1), 8);
        assert_eq!(left_child_id(13, 1), 12);
    }

    #[test]
    fn test_right_child_id() {
        assert_eq!(right_child_id(7, 3), 11);

        assert_eq!(right_child_id(3, 2), 5);
        assert_eq!(right_child_id(11, 2), 13);

        assert_eq!(right_child_id(1, 1), 2);
        assert_eq!(right_child_id(5, 1), 6);
        assert_eq!(right_child_id(9, 1), 10);
        assert_eq!(right_child_id(13, 1), 14);
    }
}
