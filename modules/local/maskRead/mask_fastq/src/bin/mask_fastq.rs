use std::io::{self, BufWriter, Write, IsTerminal};
use needletail::parse_fastx_stdin;
use flate2::{Compression, write::GzEncoder};
use clap::Parser;
use mask_fastq::mask_sequence;

/// Mask low-complexity regions in FASTQ reads using entropy calculation
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Window size for entropy calculation
    #[arg(short = 'w', long, default_value_t = 25)]
    window: usize,

    /// Entropy threshold (mask if entropy < threshold)
    #[arg(short = 'e', long, default_value_t = 0.55)]
    entropy: f64,

    /// K-mer size for entropy calculation (maximum k=8 for optimized u16 encoding)
    #[arg(short = 'k', long, default_value_t = 5)]
    kmer: usize,

    /// Gzip compression level (0-9, where 0=no compression, 1=fast/default, 9=max compression)
    #[arg(short = 'c', long, default_value_t = 1)]
    compression_level: u32,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Validate k-mer size (u16 encoding supports up to k=8)
    if args.kmer > 8 {
        eprintln!("Error: k-mer size k={} exceeds maximum supported value (k ≤ 8)", args.kmer);
        eprintln!("The optimized u16 encoding uses 2 bits per base, limiting k to 8 bases (16 bits).");
        eprintln!("For low-complexity masking, k=3 to k=7 is typically used.");
        std::process::exit(1);
    }

    if args.kmer < 1 {
        eprintln!("Error: k-mer size k={} is too small (k must be at least 1)", args.kmer);
        std::process::exit(1);
    }

    // Validate compression level
    if args.compression_level > 9 {
        eprintln!("Error: compression level {} is invalid (must be 0-9)", args.compression_level);
        std::process::exit(1);
    }

    // Check if stdin is a terminal (no piped input)
    if std::io::stdin().is_terminal() {
        eprintln!("Error: No input provided. This tool reads FASTQ data from stdin.");
        eprintln!();
        eprintln!("Usage:");
        eprintln!("  cat input.fastq | mask_fastq [OPTIONS] > output.fastq.gz");
        eprintln!("  zcat input.fastq.gz | mask_fastq [OPTIONS] > output.fastq.gz");
        eprintln!();
        eprintln!("Example:");
        eprintln!("  cat reads.fastq | mask_fastq -w 25 -e 0.55 -k 5 > masked.fastq.gz");
        eprintln!();
        eprintln!("For full help, use: mask_fastq --help");
        std::process::exit(1);
    }

    // Create gzip encoder for stdout
    let stdout = io::stdout();
    let gz_writer = GzEncoder::new(stdout, Compression::new(args.compression_level));
    let mut writer = BufWriter::new(gz_writer);

    // Parse FASTQ from stdin (handles both plain and gzipped input)
    let mut reader = parse_fastx_stdin()?;

    while let Some(record) = reader.next() {
        let rec = record?;

        // Get sequence and quality
        let sequence = rec.seq();
        let quality = rec.qual().unwrap_or(&[]);

        // Mask low-complexity regions
        let (masked_seq, masked_qual) = mask_sequence(
            sequence.as_ref(),
            quality,
            args.window,
            args.entropy,
            args.kmer
        );

        // Write masked record in FASTQ format
        writeln!(writer, "@{}", String::from_utf8_lossy(rec.id()))?;
        writeln!(writer, "{}", String::from_utf8_lossy(&masked_seq))?;
        writeln!(writer, "+")?;
        writeln!(writer, "{}", String::from_utf8_lossy(&masked_qual))?;
    }

    writer.flush()?;
    Ok(())
}
