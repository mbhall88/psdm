use anyhow::{anyhow, Context, Result};
use noodles_fasta as fasta;
use std::collections::HashSet;
use std::io::BufRead;
use std::iter::FromIterator;
use structopt::StructOpt;

trait SortExt<T> {
    fn argsort(&self) -> Vec<usize>;
    fn sort_by_indices(&mut self, indices: &mut Vec<usize>);
}

impl<T: Ord + Clone> SortExt<T> for Vec<T> {
    fn argsort(&self) -> Vec<usize> {
        let mut indices = (0..self.len()).collect::<Vec<_>>();
        indices.sort_by_key(|&i| &self[i]);
        indices
    }

    fn sort_by_indices(&mut self, indices: &mut Vec<usize>) {
        for idx in 0..self.len() {
            if indices[idx] != usize::MAX {
                let mut current_idx = idx;
                loop {
                    let target_idx = indices[current_idx];
                    indices[current_idx] = usize::MAX;
                    if indices[target_idx] == usize::MAX {
                        break;
                    }
                    self.swap(current_idx, target_idx);
                    current_idx = target_idx;
                }
            }
        }
    }
}

fn parse_ignored_chars(s: &str) -> HashSet<u8> {
    HashSet::from_iter(s.as_bytes().to_vec())
}

/// A struct to hold all of the options for the transforming sequences
#[derive(StructOpt, Debug, Default)]
pub struct Transformer {
    /// Ignore case - i.e., dist(a, A) = 0
    #[structopt(short, long)]
    ignore_case: bool,
    /// Sort the alignment(s) by ID
    #[structopt(short, long)]
    sort: bool,
    /// String of characters to ignore - e.g., `-e NX` -> dist(A, N) = 0 and dist(A, X) = 0
    #[structopt(short = "e", long, default_value="N", parse(from_str=parse_ignored_chars))]
    ignored_chars: HashSet<u8>,
}

impl Transformer {
    pub fn load_alignment<R: BufRead>(
        &self,
        reader: &mut fasta::Reader<R>,
        starting_seqlen: usize,
    ) -> Result<(Vec<String>, Vec<Vec<u8>>), anyhow::Error> {
        let mut seqlen: usize = starting_seqlen;
        let mut names: Vec<String> = vec![];
        let mut seqs: Vec<Vec<u8>> = vec![];

        for result in reader.records() {
            let record = result.context("Failed to parse record")?;
            names.push(record.name().to_owned());
            if seqlen > 0 && seqlen != record.sequence().len() {
                return Err(anyhow!(format!(
                    "Alignment sequences must all be the same length [id: {}]",
                    record.name()
                )));
            } else if seqlen == 0 {
                seqlen = record.sequence().len();
            }
            seqs.push(record.sequence().to_vec());
        }

        if self.sort {
            let mut indices = names.argsort();
            names.sort();
            seqs.sort_by_indices(&mut indices);
        }

        Ok((names, seqs))
    }
}
#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn parse_chars() {
        let s = "N-X";

        let actual = parse_ignored_chars(s);
        let expected: HashSet<u8> = HashSet::from_iter([b'N', b'-', b'X']);

        assert_eq!(actual, expected)
    }

    #[test]
    fn alignments_all_have_same_length() {
        let data = b">s1\nACGT\n>s0\nCCCC\n";
        let mut reader = fasta::Reader::new(&data[..]);
        let t: Transformer = Default::default();

        let actual = t.load_alignment(&mut reader, 0).unwrap();
        let expected = (
            vec!["s1".to_string(), "s0".to_string()],
            vec![b"ACGT".to_vec(), b"CCCC".to_vec()],
        );

        assert_eq!(actual, expected)
    }

    #[test]
    fn alignments_do_not_all_have_same_length() {
        let data = b">s0\nACGT\n>s1\nCCCCC\n";
        let mut reader = fasta::Reader::new(&data[..]);
        let t: Transformer = Default::default();

        let actual = t.load_alignment(&mut reader, 0).unwrap_err();
        assert!(actual.to_string().contains("[id: s1]"))
    }

    #[test]
    fn alignments_do_not_have_same_length_as_starting_seqlen() {
        let data = b">s0\nACGT\n>s1\nCCCC\n";
        let mut reader = fasta::Reader::new(&data[..]);
        let t: Transformer = Default::default();

        let actual = t.load_alignment(&mut reader, 1).unwrap_err();
        assert!(actual.to_string().contains("[id: s0]"))
    }

    #[test]
    fn alignments_sorted_by_id() {
        let data = b">s10\nACGT\n>s51\nCCCC\n>s0\nGGCC\n";
        let mut reader = fasta::Reader::new(&data[..]);
        let t: Transformer = Transformer {
            sort: true,
            ..Default::default()
        };

        let actual = t.load_alignment(&mut reader, 0).unwrap();
        let expected = (
            vec!["s0".to_string(), "s10".to_string(), "s51".to_string()],
            vec![b"GGCC".to_vec(), b"ACGT".to_vec(), b"CCCC".to_vec()],
        );
        assert_eq!(actual, expected)
    }

    #[test]
    fn argsort() {
        let v = vec![1, 7, 4, 2];
        let i = v.argsort();
        assert_eq!(i, &[0, 3, 2, 1]);
    }

    #[test]
    fn argsort_with_str() {
        let v = vec!["a", "c", "B", "-"];
        let i = v.argsort();
        assert_eq!(i, &[3, 2, 0, 1]);
    }

    #[test]
    fn argsort_on_empty() {
        let v: Vec<u8> = vec![];
        let i = v.argsort();
        assert!(i.is_empty());
    }

    #[test]
    fn sort_by_index() {
        let mut i = vec![5, 3, 2, 1, 6, 4, 0];
        let mut v = vec!["a", "b", "c", "d", "e", "f", "g"];
        v.sort_by_indices(&mut i);

        assert_eq!(v, &["f", "d", "c", "b", "g", "e", "a"]);
    }
}
