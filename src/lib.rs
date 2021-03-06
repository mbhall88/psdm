use anyhow::{anyhow, Context, Result};
use clap::Parser;
use itertools::iproduct;
use ndarray::{ArrayBase, Ix2, OwnedRepr};
use noodles_fasta as fasta;
use std::collections::HashSet;
use std::io::{BufRead, Error, Write};
use std::iter::FromIterator;

const IGNORE: u8 = b'.';

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

// A struct to hold all of the options for the transforming sequences
#[derive(Parser, Debug, Default)]
pub struct Transformer {
    /// Case matters - i.e., dist(a, A) = 1
    #[clap(short, long)]
    case_sensitive: bool,
    /// Sort the alignment(s) by ID
    #[clap(short, long)]
    sort: bool,
    /// String of characters to ignore - e.g., `-e N-` -> dist(A, N) = 0 and dist(A, -) = 0
    ///
    /// Note, if using `--case-sensitive` the upper- and lower-case form of a character is needed.
    /// To not ignore any characters, use `-e ''` or `-e ""`
    #[clap(short = 'e', long, default_value="N-", parse(from_str=parse_ignored_chars), allow_hyphen_values = true)]
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

        let skip_transform = self.ignored_chars.is_empty() && !self.case_sensitive;
        if !skip_transform {
            for seq in seqs.iter_mut() {
                self.transform(seq);
            }
        }

        Ok((names, seqs))
    }

    fn transform(&self, seq: &mut Vec<u8>) {
        for b in seq {
            if !self.case_sensitive {
                b.make_ascii_uppercase();
            }
            if self.ignored_chars.contains(b) {
                *b = IGNORE.to_owned();
            }
        }
    }
}

fn dist(a: u8, b: u8) -> u64 {
    (a != b && a != IGNORE && b != IGNORE) as u64
}

pub fn hamming_distance(a: &[u8], b: &[u8]) -> u64 {
    a.iter().zip(b).fold(0, |acc, (x, y)| acc + dist(*x, *y))
}

pub trait ToTable {
    fn to_csv(
        &self,
        ostream: &mut Box<dyn Write>,
        delimiter: char,
        column_names: &[String],
        row_names: &[String],
    ) -> Result<(), Error>;
    fn to_long(
        &self,
        ostream: &mut Box<dyn Write>,
        delimiter: char,
        column_names: &[String],
        row_names: &[String],
    ) -> Result<(), Error>;
}

impl ToTable for ArrayBase<OwnedRepr<u64>, Ix2> {
    fn to_csv(
        &self,
        ostream: &mut Box<dyn Write>,
        delimiter: char,
        column_names: &[String],
        row_names: &[String],
    ) -> Result<(), Error> {
        // write empty top-left corner cell
        write!(ostream, "{}", delimiter)?;
        let header: String = column_names.to_vec().join(&delimiter.to_string());
        writeln!(ostream, "{}", header)?;

        for (row_idx, row_name) in row_names.iter().enumerate() {
            write!(ostream, "{}", row_name)?;
            let row = self.row(row_idx);
            let s = row
                .iter()
                .map(|x| format!("{}{}", delimiter, x))
                .collect::<String>();
            writeln!(ostream, "{}", s)?;
        }
        Ok(())
    }

    fn to_long(
        &self,
        ostream: &mut Box<dyn Write>,
        delimiter: char,
        column_names: &[String],
        row_names: &[String],
    ) -> Result<(), Error> {
        for (i, j) in iproduct!(0..column_names.len(), 0..row_names.len()) {
            let dist = self[[j, i]];
            let c_name = &column_names[i];
            let r_name = &row_names[j];
            writeln!(ostream, "{}{d}{}{d}{}", c_name, r_name, dist, d = delimiter)?;
        }
        Ok(())
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

    #[test]
    fn transform_preserve_case() {
        let mut s = b"aC-t".to_vec();
        let expected = s.clone();
        let t: Transformer = Transformer {
            case_sensitive: true,
            ..Default::default()
        };

        t.transform(&mut s);

        assert_eq!(s, expected)
    }

    #[test]
    fn transform_ignore_case() {
        let mut s = b"aC-t".to_vec();
        let t = Transformer {
            case_sensitive: false,
            ..Default::default()
        };

        t.transform(&mut s);
        let expected = b"AC-T".to_vec();

        assert_eq!(s, expected)
    }

    #[test]
    fn transform_ignore_chars() {
        let ignore = HashSet::from_iter(b"N-x".to_vec());
        let t = Transformer {
            ignored_chars: ignore,
            case_sensitive: true,
            ..Default::default()
        };
        let mut s = b"AxC-GNt".to_vec();

        t.transform(&mut s);
        let expected = vec![b'A', IGNORE, b'C', IGNORE, b'G', IGNORE, b't'];

        assert_eq!(s, expected)
    }

    #[test]
    fn transform_ignore_chars_case_sensitive() {
        let ignore = HashSet::from_iter(b"N".to_vec());
        let t = Transformer {
            ignored_chars: ignore,
            case_sensitive: true,
            ..Default::default()
        };
        let mut s = b"ACGnt".to_vec();

        t.transform(&mut s);
        let expected = vec![b'A', b'C', b'G', b'n', b't'];

        assert_eq!(s, expected)
    }

    #[test]
    fn transform_ignore_chars_case_insensitive() {
        let ignore = HashSet::from_iter(b"N".to_vec());
        let t = Transformer {
            ignored_chars: ignore,
            case_sensitive: false,
            ..Default::default()
        };
        let mut s = b"ACGnt".to_vec();

        t.transform(&mut s);
        let expected = vec![b'A', b'C', b'G', IGNORE, b'T'];

        assert_eq!(s, expected)
    }

    #[test]
    fn test_hamming_distance() {
        let a = vec![b'A', IGNORE, b't', b'C', b'-'];
        let b = vec![b'A', b'T', b'T', b'C', b'G'];

        let actual = hamming_distance(&a, &b);
        let expected = 2;

        assert_eq!(actual, expected)
    }
}
