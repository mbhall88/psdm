use itertools::iproduct;
use rayon::prelude::*;
use std::ffi::{OsStr, OsString};
use std::fs::File;
use std::io::{stdout, BufReader, BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use noodles_fasta as fasta;
use psdm::{hamming_distance, Transformer};
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

    /// Number of threads to use. Setting to 0 will use all available
    #[structopt(short, long, default_value = "1")]
    threads: usize,

    #[structopt(flatten)]
    transformer: Transformer,
}

fn main() -> Result<()> {
    let opts = Opt::from_args();

    // set the global default number of threads for rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(opts.threads)
        .build_global()?;

    let mut ostream: Box<dyn Write> = match opts.output {
        None => Box::new(stdout()),
        Some(p) => {
            let file = File::create(p).context("Failed to create output file")?;
            Box::new(BufWriter::new(file))
        }
    };

    let mut reader1 = niffler::from_path(&opts.alignments[0])
        .map(|(r, _)| BufReader::new(r))
        .map(fasta::Reader::new)
        .context("Could not open first alignment file")?;

    let (_names1, seqs1) = opts
        .transformer
        .load_alignment(&mut reader1, 0)
        .context("Failed to load first alignment file")?;

    let (_names2, seqs2) = match opts.alignments.get(1) {
        Some(p) => {
            let mut reader2 = niffler::from_path(&p)
                .map(|(r, _)| BufReader::new(r))
                .map(fasta::Reader::new)
                .context("Could not open second alignment file")?;
            let (n, s) = opts
                .transformer
                .load_alignment(&mut reader2, seqs1[0].len())
                .context("Failed to load second alignment file")?;
            (Some(n), Some(s))
        }
        None => (None, None),
    };

    let n_seqs1 = seqs1.len();
    let n_seqs2 = match seqs2 {
        None => n_seqs1,
        Some(ref s) => s.len(),
    };

    let pairwise_indices: Vec<(usize, usize)> = iproduct!(0..n_seqs1, 0..n_seqs2).collect();

    let dists = pairwise_indices
        .into_par_iter()
        .map(|(i, j)| {
            match &seqs2 {
                None if i == j => 0, // distance between a sequence and itself
                None => hamming_distance(&seqs1[i], &seqs1[j]),
                Some(ref s) => hamming_distance(&seqs1[i], &s[j]),
            }
        })
        .collect::<Vec<u64>>();

    writeln!(&mut ostream, "{}", format!("{:?}", dists))?;

    Ok(())
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
}
