extern crate flate2;
extern crate bio;
extern crate pico_args;
extern crate rand;
extern crate rand_chacha;

use rand::SeedableRng;
use rand::seq::IteratorRandom;

fn main() {

    let mut args = pico_args::Arguments::from_env();

    let input_path: std::path::PathBuf = args.value_from_str("--input")
        .expect("could not get value of `--input` flag");
    let output_path: std::path::PathBuf = args.value_from_str("--output")
        .expect("could not get value of `--output` flag");
    let min_size: usize = args.opt_value_from_str("--min-size")
        .expect("could not get value of `--min-size` flag")
        .unwrap_or_default();
    let seed: u64 = args.opt_value_from_str("--seed")
        .expect("could not get value of `--seed` flag")
        .unwrap_or(42);
    let total: usize = args.opt_value_from_str("--total")
        .expect("could not get value of `--total`")
        .unwrap_or(5000);

    // initialize the reader
    let mut reader: Box<dyn std::io::Read> = std::fs::File::open(&input_path)
        .map(Box::new)
        .expect("could not open input path");
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
    let mut selected = Vec::with_capacity(total);
    seq_ids.into_iter().choose_multiple_fill(&mut rng, selected.as_mut());

    // hash the selected ids for faster lookup
    let selected_index = selected.iter().cloned().collect::<std::collections::HashSet<_>>();

    // recover the input sequences
    let mut reader: Box<dyn std::io::Read> = std::fs::File::open(&input_path)
        .map(Box::new)
        .expect("could not open input path");
    if input_path.ends_with(".gz") {
        reader = Box::new(flate2::read::GzDecoder::new(reader));
    }

    // write to output file
    let mut writer = bio::io::fasta::Writer::to_file(output_path)
        .expect("could not create output path");
    for res in bio::io::fasta::Reader::new(reader).records() {
        let record = res.unwrap();
        if selected_index.contains(record.id()) {
            writer.write_record(&record).expect("failed to write record");
        }
    }
}
