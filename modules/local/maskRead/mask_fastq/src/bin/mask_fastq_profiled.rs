use std::io::{self, BufWriter, Write, IsTerminal};
use std::fs::File;
use std::time::Instant;
use needletail::{parse_fastx_stdin, parse_fastx_file};
use flate2::{Compression, write::GzEncoder};
use clap::Parser;
use mask_fastq::mask_sequence_array;

/// Profiled version of mask_fastq_array with timing instrumentation
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input FASTQ file (plain or gzipped). If not specified, reads from stdin
    #[arg(short = 'i', long)]
    input: Option<String>,

    /// Output FASTQ file. If not specified, writes to stdout
    #[arg(short = 'o', long)]
    output: Option<String>,

    /// Window size for entropy calculation
    #[arg(short = 'w', long, default_value_t = 25)]
    window: usize,

    /// Entropy threshold (mask if entropy < threshold)
    #[arg(short = 'e', long, default_value_t = 0.55)]
    entropy: f64,

    /// K-mer size for entropy calculation
    #[arg(short = 'k', long, default_value_t = 5)]
    kmer: usize,

    /// Gzip compression level
    #[arg(short = 'c', long)]
    compression_level: Option<u32>,

    /// Print detailed profiling information
    #[arg(long, default_value_t = true)]
    profile: bool,
}

struct ProfilingStats {
    io_read_time: u128,
    masking_time: u128,
    io_write_time: u128,
    total_reads: usize,
    total_bases: usize,
}

impl ProfilingStats {
    fn new() -> Self {
        Self {
            io_read_time: 0,
            masking_time: 0,
            io_write_time: 0,
            total_reads: 0,
            total_bases: 0,
        }
    }

    fn print(&self, total_time: u128) {
        let total_measured = self.io_read_time + self.masking_time + self.io_write_time;
        let other_time = total_time.saturating_sub(total_measured);

        eprintln!("\n========================================");
        eprintln!("Profiling Results");
        eprintln!("========================================");
        eprintln!("Total reads processed: {}", self.total_reads);
        eprintln!("Total bases processed: {}", self.total_bases);
        eprintln!();
        eprintln!("Time breakdown:");
        eprintln!("  I/O Reading:  {:8} ms  ({:5.1}%)",
            self.io_read_time,
            100.0 * self.io_read_time as f64 / total_time as f64);
        eprintln!("  Masking:      {:8} ms  ({:5.1}%)",
            self.masking_time,
            100.0 * self.masking_time as f64 / total_time as f64);
        eprintln!("  I/O Writing:  {:8} ms  ({:5.1}%)",
            self.io_write_time,
            100.0 * self.io_write_time as f64 / total_time as f64);
        eprintln!("  Other:        {:8} ms  ({:5.1}%)",
            other_time,
            100.0 * other_time as f64 / total_time as f64);
        eprintln!("  Total:        {:8} ms", total_time);
        eprintln!();
        eprintln!("Performance:");
        eprintln!("  Throughput:   {:.1} reads/sec",
            self.total_reads as f64 / (total_time as f64 / 1000.0));
        eprintln!("  Throughput:   {:.1} Mbases/sec",
            self.total_bases as f64 / (total_time as f64 / 1000.0) / 1_000_000.0);
        eprintln!("  Per-read avg: {:.3} ms/read",
            total_time as f64 / self.total_reads as f64);
        eprintln!();

        if self.masking_time > 0 {
            eprintln!("Masking breakdown:");
            eprintln!("  Masking only: {:.3} ms/read",
                self.masking_time as f64 / self.total_reads as f64);
            eprintln!("  Masking rate: {:.1} Mbases/sec",
                self.total_bases as f64 / (self.masking_time as f64 / 1000.0) / 1_000_000.0);
        }
        eprintln!("========================================\n");
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let profile = args.profile;

    // Validate args (same as original)
    if args.kmer > 8 || args.kmer < 1 {
        eprintln!("Error: k-mer size must be 1-8");
        std::process::exit(1);
    }

    if let Some(level) = args.compression_level {
        if level > 9 {
            eprintln!("Error: compression level must be 0-9");
            std::process::exit(1);
        }
    }

    if args.input.is_none() && std::io::stdin().is_terminal() {
        eprintln!("Error: No input provided");
        std::process::exit(1);
    }

    let total_start = Instant::now();
    let mut stats = ProfilingStats::new();

    // Create reader
    let mut reader = if let Some(input_path) = &args.input {
        parse_fastx_file(input_path)?
    } else {
        parse_fastx_stdin()?
    };

    // Create writer
    let writer: Box<dyn Write> = if let Some(output_path) = &args.output {
        let output_file = File::create(output_path)?;
        let should_compress = match args.compression_level {
            Some(0) => false,
            Some(_) => true,
            None => output_path.ends_with(".gz"),
        };

        if should_compress {
            let level = args.compression_level.unwrap_or(1);
            Box::new(BufWriter::new(GzEncoder::new(output_file, Compression::new(level))))
        } else {
            Box::new(BufWriter::new(output_file))
        }
    } else {
        let should_compress = match args.compression_level {
            Some(0) | None => false,
            Some(_) => true,
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

    // Process reads with timing
    while let Some(record) = {
        let start = Instant::now();
        let rec = reader.next();
        stats.io_read_time += start.elapsed().as_millis();
        rec
    } {
        let rec = record?;

        // Get sequence and quality
        let sequence = rec.seq();
        let quality = rec.qual().unwrap_or(&[]);
        stats.total_reads += 1;
        stats.total_bases += sequence.len();

        // Mask with timing
        let (masked_seq, masked_qual) = {
            let start = Instant::now();
            let result = mask_sequence_array(
                sequence.as_ref(),
                quality,
                args.window,
                args.entropy,
                args.kmer
            );
            stats.masking_time += start.elapsed().as_millis();
            result
        };

        // Write with timing
        let start = Instant::now();
        writeln!(writer, "@{}", String::from_utf8_lossy(rec.id()))?;
        writeln!(writer, "{}", String::from_utf8_lossy(&masked_seq))?;
        writeln!(writer, "+")?;
        writeln!(writer, "{}", String::from_utf8_lossy(&masked_qual))?;
        stats.io_write_time += start.elapsed().as_millis();
    }

    writer.flush()?;

    let total_time = total_start.elapsed().as_millis();

    if profile {
        stats.print(total_time);
    }

    Ok(())
}
