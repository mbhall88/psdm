use itertools::{iproduct, Itertools};
use ndarray::Array;
use rayon::prelude::*;
use std::ffi::{OsStr, OsString};
use std::fs::File;
use std::io::{stdout, BufReader, BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use noodles_fasta as fasta;
use psdm::{hamming_distance, ToTable, Transformer};
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

fn parse_delim(s: &str) -> Result<char, String> {
    let strip = &['\'', '"', ' '][..];
    let stripped = s.replace(strip, "").replace("\\\\", "\\");
    let chars: Vec<char> = stripped.chars().collect();
    if chars.len() > 1 {
        if chars[0] == '\\' && chars[1] == 't' {
            Ok('\t')
        } else {
            Err(format!("Too many delimiter characters given {:?}", chars))
        }
    } else {
        Ok(chars[0])
    }
}

/// Compute a pairwise SNP distance matrix from one or two alignment(s)
#[derive(StructOpt, Debug)]
#[structopt()]
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

    /// Output as long-form ("melted") table
    ///
    /// By default the output is a N x N or N x M table where N is the number of sequences in the
    /// first alignment and M is the number of sequences in the (optional) second alignment.
    #[structopt(short, long = "long")]
    long_form: bool,

    /// Delimiting character for the output table
    #[structopt(short, long = "delim", default_value = ",", parse(try_from_str=parse_delim))]
    delimiter: char,

    /// Show a progress bar
    #[structopt(short = "P", long = "progress")]
    show_progress: bool,

    #[structopt(flatten)]
    transformer: Transformer,
}

fn main() -> Result<()> {
    let opts = Opt::from_args();
    // todo: add logging
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

    let (names1, seqs1) = opts
        .transformer
        .load_alignment(&mut reader1, 0)
        .context("Failed to load first alignment file")?;

    let (names2, seqs2) = match opts.alignments.get(1) {
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
    let n_seqs2: usize = match seqs2 {
        None => 0,
        Some(ref s) => s.len(),
    };

    // for intra-alignment distances, we don't need to compute the whole NxN matrix so we just
    // generate the lower-left triangle (and the diagonal for labelling reasons).
    let pairwise_indices: Vec<Vec<usize>> = match n_seqs2 {
        0 => (0..n_seqs1).combinations_with_replacement(2).collect(),
        i => iproduct!(0..n_seqs1, 0..i)
            .map(|t| vec![t.0, t.1])
            .collect(),
    };

    let pb = if opts.show_progress {
        let style = ProgressStyle::default_bar()
            .template("[{pos}/{len} ({percent}%) comparisons] {bar:40} [ETA {eta_precise}]");
        ProgressBar::new(pairwise_indices.len() as u64).with_style(style)
    } else {
        ProgressBar::hidden()
    };
    let dists: Vec<u64> = pairwise_indices
        .as_slice()
        .into_par_iter()
        .progress_with(pb)
        .map(|ix| {
            let i = ix[0];
            let j = ix[1];
            match &seqs2 {
                None if i == j => 0, // distance between a sequence and itself
                None => hamming_distance(&seqs1[i], &seqs1[j]),
                Some(ref s) => hamming_distance(&seqs1[i], &s[j]),
            }
        })
        .collect();

    let matrix =
        if n_seqs2 > 0 {
            Array::from_shape_vec((n_seqs1, n_seqs2), dists).context(
            "Failed to create matrix. This shouldn't happen, please raise an issue on GitHub",
        )?.t().to_owned()
        } else {
            let mut mtx = Array::zeros((n_seqs1, n_seqs1));
            for (ix, d) in pairwise_indices.iter().zip(dists) {
                let i = ix[0];
                let j = ix[1];
                mtx[[i, j]] = d;
                if i != j {
                    mtx[[j, i]] = d;
                }
            }
            mtx
        };

    let row_names: &Vec<String> = match &names2 {
        Some(n) => n,
        None => &names1,
    };
    let col_names: &Vec<String> = &names1;

    if opts.long_form {
        matrix
            .to_long(&mut ostream, opts.delimiter, col_names, row_names)
            .context("Failed to write output table")?;
    } else {
        matrix
            .to_csv(&mut ostream, opts.delimiter, col_names, row_names)
            .context("Failed to write output table")?;
    }

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
