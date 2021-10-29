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

    let mut reader1 = niffler::from_path(&opt.alignments[0])
        .map(|(r, _)| BufReader::new(r))
        .map(fasta::Reader::new)
        .context("Could not open first alignment file")?;

    let (names1, seqs1) = load_alignment(&mut reader1)?;

    for n in names1 {
        writeln!(&mut ostream, "{}", n)?;
    }

    for s in seqs1 {
        writeln!(&mut ostream, "{}", s.len())?;
    }

    

    Ok(())
}

fn load_alignment<R: BufRead>(
    reader: &mut fasta::Reader<R>,
) -> Result<(Vec<String>, Vec<Vec<u8>>), anyhow::Error> {
    let mut seqlen: usize = 0;
    let mut names: Vec<String> = vec![];
    let mut seqs: Vec<Vec<u8>> = vec![];

    for result in reader.records() {
        let record = result.context("Failed to parse record")?;
        names.push(record.name().to_owned());
        if seqlen > 0 && seqlen != record.sequence().len() {
            return Err(anyhow!(
                "Alignment sequences must all be the same length".to_string()
            ));
        } else if seqlen == 0 {
            seqlen = record.sequence().len();
        }
        seqs.push(record.sequence().to_vec());
    }
    Ok((names, seqs))
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
