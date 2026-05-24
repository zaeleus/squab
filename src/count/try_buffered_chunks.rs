use std::io;

pub(super) struct TryBufferedChunks<I> {
    iter: I,
    chunk_size: usize,
}

impl<I> TryBufferedChunks<I> {
    pub(super) fn new(iter: I, chunk_size: usize) -> Self {
        Self { iter, chunk_size }
    }
}

impl<I, T> Iterator for TryBufferedChunks<I>
where
    I: Iterator<Item = io::Result<T>>,
{
    type Item = io::Result<Vec<T>>;

    fn next(&mut self) -> Option<Self::Item> {
        let result: Self::Item = self.iter.by_ref().take(self.chunk_size).collect();

        match result {
            Ok(items) if items.is_empty() => None,
            Ok(items) => Some(Ok(items)),
            Err(e) => Some(Err(e)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_next() -> io::Result<()> {
        let list = [Ok(0), Ok(1), Ok(2)];
        let buffered_chunks = TryBufferedChunks::new(list.into_iter(), 2);

        let actual: Vec<_> = buffered_chunks.collect::<io::Result<_>>()?;
        let expected = [vec![0, 1], vec![2]];

        assert_eq!(actual, expected);

        Ok(())
    }
}
