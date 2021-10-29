use assert_cmd::Command;
use std::io::Write;

#[test]
fn input_file_does_not_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let err_msg = cmd.arg("nonexistent.fa").unwrap_err().to_string();

    assert!(err_msg.contains("does not exist"));

    Ok(())
}

#[test]
fn more_than_two_input_files_fails() -> Result<(), Box<dyn std::error::Error>> {
    let text = ">s0\nACGT\n>s1\nGCCC\n";
    let mut file = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
    file.write_all(text.as_bytes()).unwrap();
    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let err_msg = cmd
        .args(&[file.path(), file.path(), file.path()])
        .unwrap_err()
        .to_string();

    assert!(err_msg.contains("expecting any more values"));

    Ok(())
}

#[test]
fn trying_to_create_output_in_nonexistent_dir() -> Result<(), Box<dyn std::error::Error>> {
    let text = ">s0\nACGT\n>s1\nGCCC\n";
    let mut file = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
    file.write_all(text.as_bytes()).unwrap();
    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let p = "foo/bar/aln.fa";

    let err_msg = cmd
        .args(&["-o", p, file.path().to_str().unwrap()])
        .unwrap_err()
        .to_string();

    assert!(err_msg.contains("Failed to create output file"));

    Ok(())
}

#[test]
fn alignment_is_not_fasta() -> Result<(), Box<dyn std::error::Error>> {
    let text = "@s0\nACGT\n@s1\nGCCC\n";
    let mut file = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
    file.write_all(text.as_bytes()).unwrap();
    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let err_msg = cmd.args(&[file.path()]).unwrap_err().to_string();

    assert!(err_msg.contains("Failed to load first alignment"));

    Ok(())
}

#[test]
fn alignment_seqs_do_not_have_same_length() -> Result<(), Box<dyn std::error::Error>> {
    let text = ">s0\nAGCGT\n>s1\nGCCC\n";
    let mut file = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
    file.write_all(text.as_bytes()).unwrap();
    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let err_msg = cmd.args(&[file.path()]).unwrap_err().to_string();

    assert!(err_msg.contains("sequences must all be the same length"));

    Ok(())
}

#[test]
fn alignment_files_do_not_have_same_length() -> Result<(), Box<dyn std::error::Error>> {
    let mut file1 = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
    file1.write_all(b">s0\nAGCT\n>s1\nGCCC\n").unwrap();
    let mut file2 = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
    file2.write_all(b">s0\nAGNCT\n>s1\nGCCCC\n").unwrap();

    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let err_msg = cmd
        .args(&[file1.path(), file2.path()])
        .unwrap_err()
        .to_string();

    assert!(err_msg.contains("sequences must all be the same length"));

    Ok(())
}
