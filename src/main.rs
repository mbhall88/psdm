use std::ffi::{OsStr, OsString};
use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use anyhow::{anyhow, Context, Result};
use noodles_fasta as fasta;
use structopt::StructOpt;

/// A utility function that allows the CLI to error if a path doesn't exist
fn path_exists<S: AsRef<OsStr> + ?Sized>(s: &S) -> Result<PathBuf, OsString> {
    let path = PathBuf::from(s);
    if path.exists() {
        Ok(path)
    } else {
        Err(OsString::from(format!("{:?} does not exist", path)))
    }
}

#[derive(StructOpt, Debug)]
struct Opt {
    /// Alignment file(s) to compute the pairwise distance for.
    ///
    /// Providing two files will compute the distances for all sequences in one file against all
    /// sequences from the other file - i.e., not between sequences in the same file.
    /// The alignment file(s) can be compressed.
    #[structopt(required = true, min_values = 1, max_values = 2, parse(try_from_os_str = path_exists))]
    alignments: Vec<PathBuf>,

    /// Output file name [default: stdout]
    #[structopt(short, long, parse(from_os_str))]
    output: Option<PathBuf>,

    /// Sort the alignment(s) sequences by ID.
    #[structopt(short, long)]
    sort: bool,
}

fn main() -> Result<()> {
    let opt = Opt::from_args();

    let mut ostream: Box<dyn Write> = match opt.output {
        None => Box::new(stdout()),
        Some(p) => {
            let file = File::create(p).context("Failed to create output file")?;
            Box::new(BufWriter::new(file))
        }
    };

    let mut reader1 = niffler::from_path(&opt.alignments[0])
        .map(|(r, _)| BufReader::new(r))
        .map(fasta::Reader::new)
        .context("Could not open first alignment file")?;

    let (names1, seqs1) =
        load_alignment(&mut reader1, 0, opt.sort).context("Failed to load first alignment file")?;

    for n in &names1 {
        writeln!(&mut ostream, "{}", n)?;
    }

    for s in &seqs1 {
        writeln!(&mut ostream, "{}", s.len())?;
    }

    let (names2, seqs2) = match opt.alignments.get(1) {
        Some(p) => {
            let mut reader2 = niffler::from_path(&p)
                .map(|(r, _)| BufReader::new(r))
                .map(fasta::Reader::new)
                .context("Could not open second alignment file")?;
            let (n, s) = load_alignment(&mut reader2, seqs1[0].len(), opt.sort)
                .context("Failed to load second alignment file")?;
            (Some(n), Some(s))
        }
        None => (None, None),
    };

    if let (Some(ns), Some(ss)) = (names2, seqs2) {
        for n in ns {
            writeln!(&mut ostream, "{}", n)?;
        }

        for s in ss {
            writeln!(&mut ostream, "{}", s.len())?;
        }
    }

    Ok(())
}

fn load_alignment<R: BufRead>(
    reader: &mut fasta::Reader<R>,
    starting_seqlen: usize,
    sort: bool,
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

    if sort {
        let mut indices = names.argsort();
        names.sort();
        seqs.sort_by_indices(&mut indices);
    }

    Ok((names, seqs))
}

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

#[cfg(test)]
mod tests {
    use std::ffi::OsStr;

    use super::*;

    #[test]
    fn check_path_exists_it_doesnt() {
        let result = path_exists(OsStr::new("fake.path"));
        assert!(result.is_err())
    }

    #[test]
    fn check_path_it_does() {
        let actual = path_exists(OsStr::new("Cargo.toml")).unwrap();
        let expected = PathBuf::from("Cargo.toml");
        assert_eq!(actual, expected)
    }

    #[test]
    fn alignments_all_have_same_length() {
        let data = b">s1\nACGT\n>s0\nCCCC\n";
        let mut reader = fasta::Reader::new(&data[..]);

        let actual = load_alignment(&mut reader, 0, false).unwrap();
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

        let actual = load_alignment(&mut reader, 0, false).unwrap_err();
        assert!(actual.to_string().contains("[id: s1]"))
    }

    #[test]
    fn alignments_do_not_have_same_length_as_starting_seqlen() {
        let data = b">s0\nACGT\n>s1\nCCCC\n";
        let mut reader = fasta::Reader::new(&data[..]);

        let actual = load_alignment(&mut reader, 1, false).unwrap_err();
        assert!(actual.to_string().contains("[id: s0]"))
    }

    #[test]
    fn alignments_sorted_by_id() {
        let data = b">s10\nACGT\n>s51\nCCCC\n>s0\nGGCC\n";
        let mut reader = fasta::Reader::new(&data[..]);

        let actual = load_alignment(&mut reader, 0, true).unwrap();
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
