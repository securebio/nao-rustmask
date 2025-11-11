use std::io::{self, BufWriter, Write, IsTerminal};
use std::fs::File;
use needletail::{parse_fastx_stdin, parse_fastx_file};
use flate2::{Compression, write::GzEncoder};
use clap::Parser;
use mask_fastq::mask_sequence_array;

/// Mask low-complexity regions in FASTQ reads using entropy calculation (array-based optimized version)
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input FASTQ file (plain or gzipped). If not specified, reads from stdin
    #[arg(short = 'i', long)]
    input: Option<String>,

    /// Output FASTQ file. If not specified, writes to stdout (gzipped)
    #[arg(short = 'o', long)]
    output: Option<String>,

    /// Window size for entropy calculation
    #[arg(short = 'w', long, default_value_t = 25)]
    window: usize,

    /// Entropy threshold (mask if entropy < threshold)
    #[arg(short = 'e', long, default_value_t = 0.55)]
    entropy: f64,

    /// K-mer size for entropy calculation (maximum k=8 for optimized u16 encoding)
    #[arg(short = 'k', long, default_value_t = 5)]
    kmer: usize,

    /// Gzip compression level (0-9, where 0=no compression, 1=fast, 9=max compression).
    /// If not specified: stdout is uncompressed, .gz files use level 1 (fast compression).
    #[arg(short = 'c', long)]
    compression_level: Option<u32>,
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

    // Validate compression level if specified
    if let Some(level) = args.compression_level {
        if level > 9 {
            eprintln!("Error: compression level {} is invalid (must be 0-9)", level);
            std::process::exit(1);
        }
    }

    // Check if stdin is a terminal and no input file specified
    if args.input.is_none() && std::io::stdin().is_terminal() {
        eprintln!("Error: No input provided. Use -i to specify input file or pipe data to stdin.");
        eprintln!();
        eprintln!("Usage:");
        eprintln!("  mask_fastq_array -i input.fastq[.gz] -o output.fastq [OPTIONS]");
        eprintln!("  cat input.fastq[.gz] | mask_fastq_array [OPTIONS] > output.fastq");
        eprintln!();
        eprintln!("Note: This is the array-based optimized version (O(1) entropy calculation)");
        eprintln!("Note: Input can be plain or gzipped FASTQ (auto-detected)");
        eprintln!();
        eprintln!("For full help, use: mask_fastq_array --help");
        std::process::exit(1);
    }

    // Create reader from file or stdin
    let mut reader = if let Some(input_path) = &args.input {
        parse_fastx_file(input_path)?
    } else {
        parse_fastx_stdin()?
    };

    // Create writer to file or stdout
    let writer: Box<dyn Write> = if let Some(output_path) = &args.output {
        let output_file = File::create(output_path)?;

        // Determine if we should compress based on extension and -c flag
        let should_compress = match args.compression_level {
            Some(0) => false,  // Explicit -c 0: no compression
            Some(_) => true,   // Explicit -c 1-9: compress
            None => output_path.ends_with(".gz"),  // No -c flag: auto-detect from extension
        };

        if should_compress {
            let level = args.compression_level.unwrap_or(1);  // Default to level 1 for .gz files
            Box::new(BufWriter::new(GzEncoder::new(output_file, Compression::new(level))))
        } else {
            Box::new(BufWriter::new(output_file))
        }
    } else {
        // stdout: compress only if -c flag is specified with value > 0
        let should_compress = match args.compression_level {
            Some(0) | None => false,  // No -c or -c 0: uncompressed (NEW default)
            Some(_) => true,          // -c 1-9: compressed
        };

        if should_compress {
            let level = args.compression_level.unwrap();
            let stdout = io::stdout();
            Box::new(BufWriter::new(GzEncoder::new(stdout, Compression::new(level))))
        } else {
            let stdout = io::stdout();
            Box::new(BufWriter::new(stdout))
        }
    };

    let mut writer = writer;

    while let Some(record) = reader.next() {
        let rec = record?;

        // Get sequence and quality
        let sequence = rec.seq();
        let quality = rec.qual().unwrap_or(&[]);

        // Mask low-complexity regions using array-based tracker
        let (masked_seq, masked_qual) = mask_sequence_array(
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
