extern crate bio;
extern crate indicatif;
extern crate pico_args;

use std::io::Write;

use indicatif::ProgressBar;
use indicatif::ProgressStyle;

fn main() {
    let mut args = pico_args::Arguments::from_env();

    let input_path: std::path::PathBuf = args
        .value_from_str("--input")
        .expect("could not get value of `--input` flag");
    let output_path: std::path::PathBuf = args
        .value_from_str("--output")
        .expect("could not get value of `--output` flag");

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
    let reader = std::fs::File::open(&input_path)
        .map(|f| pb.wrap_read(f))
        .expect("could not open input file");

    // create a temporary work directory
    let workdir = std::path::PathBuf::from(".genomes");
    if !workdir.exists() {
        std::fs::create_dir(&workdir).unwrap();
    }

    // create one temporary file per sequence in input file
    let mut paths = Vec::new();
    for res in bio::io::fasta::Reader::new(reader).records() {
        let record = res.unwrap();

        let genome_file = workdir.join(format!("{}.fa", record.id()));
        let mut contig_writer = bio::io::fasta::Writer::to_file(&genome_file)
            .expect("could not create file");
        contig_writer.write_record(&record).unwrap();

        paths.push(genome_file);
    }

    // write genome list
    let mut list_file = std::fs::File::create(".genomes.txt")
        .expect("could not create temporary file");

    for path in paths.iter() {
        writeln!(list_file, "{}", path.display()).unwrap();
    }
    list_file.flush().unwrap();

    // run FastANI
    let proc = std::process::Command::new("/home/larralde/.local/bin/fastANI")
        .arg("--ql")
        .arg(".genomes.txt")
        .arg("--rl")
        .arg(".genomes.txt")
        .arg("-o")
        .arg(&output_path)
        .stdout(std::process::Stdio::inherit())
        .stderr(std::process::Stdio::inherit())
        .output()
        .expect("failed to run fastANI");

    std::fs::remove_file(".genomes.txt");
    std::fs::remove_dir_all(workdir);

    assert!(proc.status.success());
}
