extern crate bio;
extern crate flate2;
extern crate indicatif;
extern crate pico_args;
extern crate rand;
extern crate rand_chacha;

use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use rand::seq::SliceRandom;
use rand::SeedableRng;

fn main() {
    let mut args = pico_args::Arguments::from_env();

    let input_path: std::path::PathBuf = args
        .value_from_str("--input")
        .expect("could not get value of `--input` flag");
    let output_path: std::path::PathBuf = args
        .value_from_str("--output")
        .expect("could not get value of `--output` flag");
    let min_size: usize = args
        .opt_value_from_str("--min-size")
        .expect("could not get value of `--min-size` flag")
        .unwrap_or_default();
    let seed: u64 = args
        .opt_value_from_str("--seed")
        .expect("could not get value of `--seed` flag")
        .unwrap_or(42);
    let total: usize = args
        .opt_value_from_str("--total")
        .expect("could not get value of `--total`")
        .unwrap_or(5000);

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

    // collect sequence IDs if above minimum size
    let mut seq_ids = Vec::new();
    for res in bio::io::fasta::Reader::new(reader).records() {
        let record = res.unwrap();
        if record.seq().len() > min_size {
            seq_ids.push(record.id().to_string());
        }
    }

    // randomly select `--total` sequences
    let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(seed);
    let selected = seq_ids
        .choose_multiple(&mut rng, total)
        .cloned()
        .collect::<std::collections::HashSet<_>>();

    // recover the input sequences
    pb.reset();
    pb.reset_eta();
    let mut reader: Box<dyn std::io::Read> = std::fs::File::open(&input_path)
        .map(|f| pb.wrap_read(f))
        .map(Box::new)
        .expect("could not open input path");
    if input_path.ends_with(".gz") {
        reader = Box::new(flate2::read::GzDecoder::new(reader));
    }

    // write to output file
    let mut writer =
        bio::io::fasta::Writer::to_file(output_path).expect("could not create output path");
    for res in bio::io::fasta::Reader::new(reader).records() {
        let record = res.unwrap();
        if selected.contains(record.id()) {
            writer
                .write_record(&record)
                .expect("failed to write record");
        }
    }
}
