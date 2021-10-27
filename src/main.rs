use anyhow::{Context, Result};
use noodles_fasta as fasta;
use std::ffi::{OsStr, OsString};
use std::fs::File;
use std::io::{stdout, BufReader, BufWriter, Write};
use std::path::PathBuf;
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
    #[structopt(required = true, min_values = 1, max_values = 2, parse(try_from_os_str = path_exists))]
    alignments: Vec<PathBuf>,

    /// Output file name [default: stdout]
    #[structopt(short, parse(from_os_str))]
    output: Option<PathBuf>,
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

    for p in opt.alignments {
        println!("{:?}", p);
        let mut reader = niffler::from_path(p)
            .map(|(r, _)| BufReader::new(r))
            .map(fasta::Reader::new)
            .context("Could not open alignment file")?;

        for result in reader.records() {
            let record = result?;
            writeln!(&mut ostream, "{}", record.sequence().len())?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::OsStr;

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
