extern crate bio;
extern crate flate2;
extern crate indicatif;
extern crate pico_args;

use std::io::BufRead;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;

fn main() {
    let mut args = pico_args::Arguments::from_env();

    let input_path: std::path::PathBuf = args
        .value_from_str("--input")
        .expect("could not get value of `--input` flag");
    let list_path: std::path::PathBuf = args
        .value_from_str("--list")
        .expect("could not get value of `--list` flag");
    let output_path: std::path::PathBuf = args
        .value_from_str("--output")
        .expect("could not get value of `--output` flag");

    // load the list of sequences to extract
    let selection = std::fs::File::open(&list_path)
        .map(std::io::BufReader::new)
        .expect("could not open list file")
        .lines()
        .map(|s| s.unwrap().trim().to_string())
        .collect::<std::collections::HashSet<String>>();

    // create a progress bar
    let pb = ProgressBar::new(
        input_path
            .metadata()
            .expect("could not find input file")
            .len() as u64,
    );
    pb.set_style(ProgressStyle::default_bar()
        .template("{msg:>30} {percent:>3}% |{wide_bar:.cyan/blue}| {bytes}/{total_bytes} [{elapsed_precise}<{eta_precise}]"));

    // initialize the reader
    let mut reader: Box<dyn std::io::Read> = std::fs::File::open(&input_path)
        .map(|f| pb.wrap_read(f))
        .map(Box::new)
        .expect("could not open input file");
    if input_path.ends_with(".gz") {
        reader = Box::new(flate2::read::GzDecoder::new(reader));
    }

    // initialize the writer
    let mut writer =
        bio::io::fasta::Writer::to_file(output_path).expect("could not create output path");

    // write back the record if it's in the selection
    for res in bio::io::fasta::Reader::new(reader).records() {
        let record = res.unwrap();
        if selection.contains(record.id()) {
            writer
                .write_record(&record)
                .expect("failed to write record");
        }
    }
}
