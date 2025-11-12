use std::collections::HashMap;
use std::time::Instant;
use clap::Parser;

/// Compare u16 vs u32 encoding performance
#[derive(Parser, Debug)]
#[command(about = "Benchmark u16 vs u32 k-mer encoding")]
struct Args {
    /// K-mer size to test (must be ≤ 8)
    #[arg(short = 'k', long, default_value_t = 7)]
    kmer: usize,

    /// Number of iterations
    #[arg(short = 'n', long, default_value_t = 10000)]
    iterations: usize,

    /// Read length
    #[arg(short = 'l', long, default_value_t = 150)]
    read_length: usize,
}

/// Encode a k-mer into u16 (current implementation)
fn encode_kmer_u16(bases: &[u8]) -> Option<u16> {
    if bases.len() > 8 {
        return None;
    }
    let mut encoded: u16 = 0;
    for &base in bases {
        let bits = match base {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => return None,
        };
        encoded = (encoded << 2) | bits;
    }
    Some(encoded)
}

/// Encode a k-mer into u32 (proposed for k>8)
fn encode_kmer_u32(bases: &[u8]) -> Option<u32> {
    if bases.len() > 15 {
        return None;
    }
    let mut encoded: u32 = 0;
    for &base in bases {
        let bits = match base {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => return None,
        };
        encoded = (encoded << 2) | bits;
    }
    Some(encoded)
}

/// Process sequence with u16 HashMap
fn process_u16(sequence: &[u8], k: usize, window: usize) -> f64 {
    let seq_len = sequence.len();
    let mut kmer_counts: HashMap<u16, usize> = HashMap::new();
    let mut entropy_sum = 0.0;

    for i in 0..seq_len {
        let window_start = if i + 1 >= window {
            i + 1 - window
        } else {
            0
        };
        let window_end = i + 1;

        if window_end - window_start < window {
            continue;
        }

        // Build k-mer counts for window
        kmer_counts.clear();
        for j in window_start..=window_end.saturating_sub(k) {
            if let Some(encoded) = encode_kmer_u16(&sequence[j..j + k]) {
                *kmer_counts.entry(encoded).or_insert(0) += 1;
            }
        }

        // Calculate entropy (simplified)
        let total_kmers = window - k + 1;
        let mut entropy = 0.0;
        for &count in kmer_counts.values() {
            if count > 0 {
                let p = count as f64 / total_kmers as f64;
                entropy -= p * p.log2();
            }
        }
        entropy_sum += entropy;
    }

    entropy_sum
}

/// Process sequence with u32 HashMap
fn process_u32(sequence: &[u8], k: usize, window: usize) -> f64 {
    let seq_len = sequence.len();
    let mut kmer_counts: HashMap<u32, usize> = HashMap::new();
    let mut entropy_sum = 0.0;

    for i in 0..seq_len {
        let window_start = if i + 1 >= window {
            i + 1 - window
        } else {
            0
        };
        let window_end = i + 1;

        if window_end - window_start < window {
            continue;
        }

        // Build k-mer counts for window
        kmer_counts.clear();
        for j in window_start..=window_end.saturating_sub(k) {
            if let Some(encoded) = encode_kmer_u32(&sequence[j..j + k]) {
                *kmer_counts.entry(encoded).or_insert(0) += 1;
            }
        }

        // Calculate entropy (simplified)
        let total_kmers = window - k + 1;
        let mut entropy = 0.0;
        for &count in kmer_counts.values() {
            if count > 0 {
                let p = count as f64 / total_kmers as f64;
                entropy -= p * p.log2();
            }
        }
        entropy_sum += entropy;
    }

    entropy_sum
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
        eprintln!("Error: This test only supports k ≤ 8 (to compare u16 vs u32)");
        std::process::exit(1);
    }

    println!("=== u16 vs u32 Encoding Benchmark ===");
    println!("K-mer size: {}", args.kmer);
    println!("Read length: {} bp", args.read_length);
    println!("Iterations: {}", args.iterations);
    println!();

    // Memory analysis
    println!("=== Memory Analysis ===");
    println!("HashMap entry size:");
    println!("  u16: key=2 bytes, value=8 bytes, overhead=~8 bytes = ~18 bytes/entry");
    println!("  u32: key=4 bytes, value=8 bytes, overhead=~8 bytes = ~20 bytes/entry");
    println!("  Difference: +2 bytes per entry (+11%)");
    println!();

    let window = 25;
    let max_kmers_in_window = window - args.kmer + 1;
    println!("For window={}, k={}: max {} k-mers in window", window, args.kmer, max_kmers_in_window);
    println!("  u16 HashMap: ~{} bytes per window", max_kmers_in_window * 18);
    println!("  u32 HashMap: ~{} bytes per window", max_kmers_in_window * 20);
    println!("  Overhead: +{} bytes", max_kmers_in_window * 2);
    println!();

    // Generate test data
    println!("Generating test sequences...");
    let sequences: Vec<_> = (0..args.iterations)
        .map(|i| generate_random_sequence(args.read_length, i as u64))
        .collect();
    println!();

    // Benchmark u16
    println!("Benchmarking u16 HashMap...");
    let start = Instant::now();
    let mut total_entropy_u16 = 0.0;

    for seq in &sequences {
        total_entropy_u16 += process_u16(seq, args.kmer, window);
    }

    let u16_time = start.elapsed();
    let u16_throughput = (args.iterations * args.read_length) as f64 / u16_time.as_secs_f64() / 1_000_000.0;

    println!("u16 time: {:.3} seconds", u16_time.as_secs_f64());
    println!("u16 throughput: {:.2} Mbp/s", u16_throughput);
    println!("u16 total entropy: {:.2}", total_entropy_u16);
    println!();

    // Benchmark u32
    println!("Benchmarking u32 HashMap...");
    let start = Instant::now();
    let mut total_entropy_u32 = 0.0;

    for seq in &sequences {
        total_entropy_u32 += process_u32(seq, args.kmer, window);
    }

    let u32_time = start.elapsed();
    let u32_throughput = (args.iterations * args.read_length) as f64 / u32_time.as_secs_f64() / 1_000_000.0;

    println!("u32 time: {:.3} seconds", u32_time.as_secs_f64());
    println!("u32 throughput: {:.2} Mbp/s", u32_throughput);
    println!("u32 total entropy: {:.2}", total_entropy_u32);
    println!();

    // Verify correctness
    let entropy_diff = (total_entropy_u16 - total_entropy_u32).abs();
    if entropy_diff < 0.01 {
        println!("✓ Results match (difference: {:.6})", entropy_diff);
    } else {
        println!("⚠ Results differ by {:.2}!", entropy_diff);
    }
    println!();

    // Compare
    let slowdown_pct = ((u32_time.as_secs_f64() / u16_time.as_secs_f64()) - 1.0) * 100.0;

    println!("=== Comparison ===");
    if slowdown_pct.abs() < 1.0 {
        println!("Performance difference: <1% (negligible)");
    } else if slowdown_pct > 0.0 {
        println!("u32 is {:.1}% SLOWER than u16", slowdown_pct);
    } else {
        println!("u32 is {:.1}% FASTER than u16", -slowdown_pct);
    }
    println!();

    if slowdown_pct.abs() < 3.0 {
        println!("✓ Recommendation: Use u32 everywhere - negligible performance impact");
        println!("  Benefits:");
        println!("  - Simpler code (one implementation)");
        println!("  - Supports k up to 15");
        println!("  - Memory overhead: only +2 bytes per k-mer in HashMap (~11%)");
    } else if slowdown_pct > 5.0 {
        println!("⚠ Recommendation: Consider keeping both u16 and u32");
        println!("  u16 provides {:.1}% performance advantage for k≤8", -slowdown_pct);
    } else {
        println!("→ Marginal difference - either approach is reasonable");
    }

    Ok(())
}
