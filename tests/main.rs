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
        .args([file.path(), file.path(), file.path()])
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
        .args(["-o", p, file.path().to_str().unwrap()])
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
    let err_msg = cmd.args([file.path()]).unwrap_err().to_string();

    assert!(err_msg.contains("Failed to load first alignment"));

    Ok(())
}

#[test]
fn alignment_seqs_do_not_have_same_length() -> Result<(), Box<dyn std::error::Error>> {
    let text = ">s0\nAGCGT\n>s1\nGCCC\n";
    let mut file = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
    file.write_all(text.as_bytes()).unwrap();
    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let err_msg = cmd.args([file.path()]).unwrap_err().to_string();

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
        .args([file1.path(), file2.path()])
        .unwrap_err()
        .to_string();

    assert!(err_msg.contains("sequences must all be the same length"));

    Ok(())
}

#[test]
fn intra_alignment_with_defaults() -> Result<(), Box<dyn std::error::Error>> {
    let aln = "tests/cases/aln1.fa";

    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let output = cmd.args(["-c", aln]).unwrap().stdout;

    let expected = b",s1,s2,s0\ns1,0,3,3\ns2,3,0,5\ns0,3,5,0\n";
    assert_eq!(output, expected);

    Ok(())
}

#[test]
fn intra_alignment_with_ignore_case_ignore_no_characters() -> Result<(), Box<dyn std::error::Error>>
{
    let aln = "tests/cases/aln1.fa";

    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let output = cmd.args(["-e ''", aln]).unwrap().stdout;

    let expected = b",s1,s2,s0\ns1,0,4,1\ns2,4,0,5\ns0,1,5,0\n";
    assert_eq!(output, expected);

    Ok(())
}

#[test]
fn intra_alignment_with_use_case_ignore_characters() -> Result<(), Box<dyn std::error::Error>> {
    let aln = "tests/cases/aln1.fa";

    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let output = cmd.args(["-ce n-", aln]).unwrap().stdout;

    let expected = b",s1,s2,s0\ns1,0,3,3\ns2,3,0,5\ns0,3,5,0\n";
    assert_eq!(output, expected);

    Ok(())
}

#[test]
fn intra_alignment_with_ignore_case_ignore_characters() -> Result<(), Box<dyn std::error::Error>> {
    let aln = "tests/cases/aln1.fa";

    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let output = cmd.args(["-e aN", aln]).unwrap().stdout;

    let expected = b",s1,s2,s0\ns1,0,2,1\ns2,2,0,3\ns0,1,3,0\n";
    assert_eq!(output, expected);

    Ok(())
}

#[test]
fn intra_alignment_with_ignore_case_ignore_characters_sorted(
) -> Result<(), Box<dyn std::error::Error>> {
    let aln = "tests/cases/aln1.fa";

    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let output = cmd.args(["-s", "-e aN", aln]).unwrap().stdout;

    let expected = b",s0,s1,s2\ns0,0,1,3\ns1,1,0,2\ns2,3,2,0\n";
    assert_eq!(output, expected);

    Ok(())
}

#[test]
fn intra_alignment_with_defaults_long_form() -> Result<(), Box<dyn std::error::Error>> {
    let aln = "tests/cases/aln1.fa";

    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let output = cmd.args(["-cl", aln]).unwrap().stdout;

    let expected = b"s1,s1,0
s1,s2,3
s1,s0,3
s2,s1,3
s2,s2,0
s2,s0,5
s0,s1,3
s0,s2,5
s0,s0,0\n";
    assert_eq!(output, expected);

    Ok(())
}

#[test]
fn intra_alignment_with_defaults_long_form_sorted() -> Result<(), Box<dyn std::error::Error>> {
    let aln = "tests/cases/aln1.fa";

    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let output = cmd.args(["-cls", aln]).unwrap().stdout;

    let expected = b"s0,s0,0
s0,s1,3
s0,s2,5
s1,s0,3
s1,s1,0
s1,s2,3
s2,s0,5
s2,s1,3
s2,s2,0\n";
    assert_eq!(output, expected);

    Ok(())
}

#[test]
fn inter_alignment_with_tab_delim() -> Result<(), Box<dyn std::error::Error>> {
    let aln1 = "tests/cases/aln1.fa";
    let aln2 = "tests/cases/aln2.fa.gz";

    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let output = cmd.args(["-cd '\t'", aln1, aln2]).unwrap().stdout;

    let expected = b"\ts1\ts2\ts0\ns2\t6\t6\t5\ns5\t1\t4\t3\n";
    assert_eq!(output, expected);

    Ok(())
}
#[test]
fn inter_alignment_in_long_form_sorted() -> Result<(), Box<dyn std::error::Error>> {
    let aln1 = "tests/cases/aln1.fa";
    let aln2 = "tests/cases/aln2.fa.gz";

    let mut cmd = Command::cargo_bin("psdm").unwrap();
    let output = cmd.args(["-lsc", aln1, aln2]).unwrap().stdout;

    let expected = b"s0,s2,5\ns0,s5,3\ns1,s2,6\ns1,s5,1\ns2,s2,6\ns2,s5,4\n";
    assert_eq!(output, expected);

    Ok(())
}
