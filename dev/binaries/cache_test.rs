use std::time::Instant;
use mask_fastq::{encode_kmer, mask_sequence_array, mask_sequence};
use clap::Parser;

/// Test cache effects for different k-mer sizes
#[derive(Parser, Debug)]
#[command(about = "Test cache effects for array vs HashMap approaches")]
struct Args {
    /// K-mer size to test
    #[arg(short = 'k', long, default_value_t = 7)]
    kmer: usize,

    /// Number of iterations
    #[arg(short = 'n', long, default_value_t = 1000)]
    iterations: usize,

    /// Read length
    #[arg(short = 'l', long, default_value_t = 150)]
    read_length: usize,
}

/// Generate a random DNA sequence
fn generate_random_sequence(length: usize, seed: u64) -> Vec<u8> {
    let mut rng = seed;
    let bases = [b'A', b'C', b'G', b'T'];

    (0..length)
        .map(|_| {
            rng = rng.wrapping_mul(1103515245).wrapping_add(12345);
            bases[(rng >> 16) as usize % 4]
        })
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    if args.kmer > 8 {
        eprintln!("Error: This build only supports k ≤ 8 (u16 encoding)");
        eprintln!("Testing array approach beyond k=8 requires u32 encoding");
        std::process::exit(1);
    }

    println!("=== Cache Performance Test ===");
    println!("K-mer size: {}", args.kmer);
    println!("Read length: {} bp", args.read_length);
    println!("Iterations: {}", args.iterations);
    println!();

    let array_size = 1usize << (2 * args.kmer); // 4^k
    let array_memory_kb = (array_size * 2) as f64 / 1024.0;

    println!("Array size: {} entries ({:.2} KB)", array_size, array_memory_kb);
    println!();

    // Estimate cache sizes (typical modern CPU)
    println!("Typical cache sizes:");
    println!("  L1: 32-64 KB per core");
    println!("  L2: 256-512 KB per core");
    println!("  L3: 8-32 MB shared");
    println!();

    if array_memory_kb < 32.0 {
        println!("✓ Array fits in L1 cache - expect excellent performance");
    } else if array_memory_kb < 256.0 {
        println!("⚠ Array exceeds L1, fits in L2 - expect good performance");
    } else if array_memory_kb < 8192.0 {
        println!("⚠ Array exceeds L2, may fit in L3 - expect degraded performance");
    } else {
        println!("✗ Array exceeds typical L3 - expect poor performance");
    }
    println!();

    // Generate test data
    let mut sequences = Vec::new();
    for i in 0..args.iterations {
        let seq = generate_random_sequence(args.read_length, i as u64);
        let qual = vec![b'I'; args.read_length];
        sequences.push((seq, qual));
    }

    // Benchmark array-based approach
    println!("Benchmarking array-based approach...");
    let start = Instant::now();
    let mut total_masked = 0;

    for (seq, qual) in &sequences {
        let (masked, _) = mask_sequence_array(seq, qual, 25, 0.55, args.kmer);
        total_masked += masked.iter().filter(|&&b| b == b'N').count();
    }

    let array_time = start.elapsed();
    let array_throughput = (args.iterations * args.read_length) as f64 / array_time.as_secs_f64() / 1_000_000.0;

    println!("Array time: {:.3} seconds", array_time.as_secs_f64());
    println!("Array throughput: {:.2} Mbp/s", array_throughput);
    println!("Masked bases: {}", total_masked);
    println!();

    // Benchmark HashMap-based approach
    println!("Benchmarking HashMap-based approach...");
    let start = Instant::now();
    let mut total_masked_hash = 0;

    for (seq, qual) in &sequences {
        let (masked, _) = mask_sequence(seq, qual, 25, 0.55, args.kmer);
        total_masked_hash += masked.iter().filter(|&&b| b == b'N').count();
    }

    let hash_time = start.elapsed();
    let hash_throughput = (args.iterations * args.read_length) as f64 / hash_time.as_secs_f64() / 1_000_000.0;

    println!("HashMap time: {:.3} seconds", hash_time.as_secs_f64());
    println!("HashMap throughput: {:.2} Mbp/s", hash_throughput);
    println!("Masked bases: {}", total_masked_hash);
    println!();

    // Compare
    let speedup = array_time.as_secs_f64() / hash_time.as_secs_f64();
    println!("=== Comparison ===");
    if speedup < 1.0 {
        println!("HashMap is {:.2}x FASTER than array", 1.0 / speedup);
    } else {
        println!("Array is {:.2}x faster than HashMap", speedup);
    }
    println!();

    if speedup < 1.1 && array_memory_kb > 64.0 {
        println!("⚠ Recommendation: Array provides minimal benefit (<10% faster) but uses {:.2}x more memory",
                 array_memory_kb / 1.0);
        println!("   Consider using HashMap approach for k={}", args.kmer);
    } else if speedup > 1.5 {
        println!("✓ Array approach provides significant speedup (>50%)");
    }

    Ok(())
}
