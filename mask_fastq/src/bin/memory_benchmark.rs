use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;
use clap::Parser;
use mask_fastq::mask_sequence_auto;
use rayon::prelude::*;

/// Benchmark memory usage and performance for different k-mer sizes
#[derive(Parser, Debug)]
#[command(author, version, about = "Memory and performance benchmarking tool")]
struct Args {
    /// K-mer size to test
    #[arg(short = 'k', long, default_value_t = 5)]
    kmer: usize,

    /// Window size
    #[arg(short = 'w', long, default_value_t = 25)]
    window: usize,

    /// Number of reads to generate
    #[arg(short = 'n', long, default_value_t = 10000)]
    num_reads: usize,

    /// Read length
    #[arg(short = 'l', long, default_value_t = 150)]
    read_length: usize,

    /// Chunk size for parallel processing
    #[arg(short = 's', long, default_value_t = 1000)]
    chunk_size: usize,

    /// Number of threads
    #[arg(short = 't', long)]
    threads: Option<usize>,

    /// Output file for synthetic FASTQ (optional, for inspection)
    #[arg(short = 'o', long)]
    output: Option<String>,
}

/// Generate a random DNA sequence
fn generate_random_sequence(length: usize, seed: u64) -> Vec<u8> {
    // Simple LCG random number generator
    let mut rng = seed;
    let bases = [b'A', b'C', b'G', b'T'];

    (0..length)
        .map(|_| {
            rng = rng.wrapping_mul(1103515245).wrapping_add(12345);
            bases[(rng >> 16) as usize % 4]
        })
        .collect()
}

/// Generate a low-complexity sequence (for testing masking)
fn generate_low_complexity_sequence(length: usize, pattern: &[u8]) -> Vec<u8> {
    pattern.iter().cycle().take(length).copied().collect()
}

/// Generate synthetic FASTQ data
fn generate_fastq_data(num_reads: usize, read_length: usize, mix_complexity: bool) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut data = Vec::with_capacity(num_reads);

    for i in 0..num_reads {
        let seq = if mix_complexity && i % 4 == 0 {
            // 25% low-complexity reads
            if i % 8 == 0 {
                generate_low_complexity_sequence(read_length, b"AAAA")
            } else {
                generate_low_complexity_sequence(read_length, b"GCGC")
            }
        } else {
            // 75% random reads
            generate_random_sequence(read_length, i as u64)
        };

        let qual = vec![b'I'; read_length]; // High quality scores
        data.push((seq, qual));
    }

    data
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Validate k-mer size
    if args.kmer > 8 {
        eprintln!("Error: k-mer size k={} exceeds maximum (k ≤ 8)", args.kmer);
        std::process::exit(1);
    }

    // Configure thread pool
    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap();
    }

    let num_threads = rayon::current_num_threads();

    println!("=== Memory Benchmark Configuration ===");
    println!("K-mer size (k): {}", args.kmer);
    println!("Window size: {}", args.window);
    println!("Number of reads: {}", args.num_reads);
    println!("Read length: {} bp", args.read_length);
    println!("Chunk size: {}", args.chunk_size);
    println!("Threads: {}", num_threads);
    println!();

    // Calculate expected memory usage
    let kmer_space = 1usize << (2 * args.kmer); // 4^k
    let tracker_memory_bytes = kmer_space * 2; // u16 counts
    let window_kmers = args.window.saturating_sub(args.kmer).saturating_add(1);
    let auxiliary_memory = (window_kmers + 2) * (2 + 8); // count_counts + entropy_table
    let per_tracker_memory = tracker_memory_bytes + auxiliary_memory;

    let per_read_memory = args.read_length * 2 + 20; // seq + qual + id overhead
    let chunk_memory = args.chunk_size * per_read_memory * 2; // input + output

    let total_expected_mb = (per_tracker_memory * num_threads + chunk_memory) as f64 / 1_048_576.0;

    println!("=== Expected Memory Usage ===");
    println!("Tracker memory per thread:");
    println!("  - K-mer counts array (4^{}): {} entries × 2 bytes = {:.2} KB",
             args.kmer, kmer_space, tracker_memory_bytes as f64 / 1024.0);
    println!("  - Auxiliary arrays: {} bytes", auxiliary_memory);
    println!("  - Total per tracker: {:.2} KB", per_tracker_memory as f64 / 1024.0);
    println!("Total tracker memory ({} threads): {:.2} MB",
             num_threads, (per_tracker_memory * num_threads) as f64 / 1_048_576.0);
    println!();
    println!("Chunk memory:");
    println!("  - Per read: ~{} bytes (seq + qual + id)", per_read_memory);
    println!("  - Chunk size: {} reads", args.chunk_size);
    println!("  - Input + output chunks: {:.2} MB", chunk_memory as f64 / 1_048_576.0);
    println!();
    println!("Estimated total working memory: {:.2} MB", total_expected_mb);
    println!();

    // Generate synthetic data
    println!("Generating synthetic FASTQ data...");
    let data = generate_fastq_data(args.num_reads, args.read_length, true);
    println!("Generated {} reads of {} bp each", data.len(), args.read_length);
    println!();

    // Optionally write to file for inspection
    if let Some(output_path) = &args.output {
        println!("Writing synthetic FASTQ to {}...", output_path);
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);

        for (i, (seq, qual)) in data.iter().enumerate() {
            writeln!(writer, "@read_{}", i)?;
            writeln!(writer, "{}", String::from_utf8_lossy(seq))?;
            writeln!(writer, "+")?;
            writeln!(writer, "{}", String::from_utf8_lossy(qual))?;
        }
        writer.flush()?;
        println!("Wrote {} reads to {}", data.len(), output_path);
        println!();
    }

    // Benchmark processing
    println!("=== Benchmarking ===");
    println!("Processing reads in chunks of {}...", args.chunk_size);

    let start = Instant::now();
    let mut total_masked = 0;
    let mut chunks_processed = 0;

    for chunk in data.chunks(args.chunk_size) {
        // Process chunk in parallel
        let results: Vec<(Vec<u8>, Vec<u8>)> = chunk
            .par_iter()
            .map(|(seq, qual)| {
                mask_sequence_auto(
                    seq,
                    qual,
                    args.window,
                    0.55,
                    args.kmer,
                )
            })
            .collect();

        // Count masked bases
        for (masked_seq, _) in results.iter() {
            total_masked += masked_seq.iter().filter(|&&b| b == b'N').count();
        }

        chunks_processed += 1;
    }

    let elapsed = start.elapsed();

    println!();
    println!("=== Results ===");
    println!("Chunks processed: {}", chunks_processed);
    println!("Total reads: {}", args.num_reads);
    println!("Total bases: {}", args.num_reads * args.read_length);
    println!("Masked bases: {}", total_masked);
    println!("Masked percentage: {:.2}%",
             (total_masked as f64 / (args.num_reads * args.read_length) as f64) * 100.0);
    println!();
    println!("Time elapsed: {:.2} seconds", elapsed.as_secs_f64());
    println!("Throughput: {:.2} reads/sec", args.num_reads as f64 / elapsed.as_secs_f64());
    println!("Throughput: {:.2} Mbp/sec",
             (args.num_reads * args.read_length) as f64 / elapsed.as_secs_f64() / 1_000_000.0);

    Ok(())
}
