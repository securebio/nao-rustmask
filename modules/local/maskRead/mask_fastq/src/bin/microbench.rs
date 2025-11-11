use std::time::Instant;
use mask_fastq::{encode_kmer, ArrayEntropyTracker};

/// Microbenchmark to identify bottlenecks in masking operations
fn main() {
    println!("========================================");
    println!("Microbenchmark: mask_fastq Components");
    println!("========================================\n");

    // Test sequence (10Kbp, typical ONT read length)
    let test_seq: Vec<u8> = (0..10000)
        .map(|i| match i % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        })
        .collect();

    let k = 5;
    let window = 25;
    let iterations = 100;

    // Benchmark 1: encode_kmer performance
    println!("Benchmark 1: encode_kmer()");
    println!("  Testing {} k-mer encodings", test_seq.len() - k + 1);

    let start = Instant::now();
    let mut count = 0;
    for _ in 0..iterations {
        for i in 0..=test_seq.len() - k {
            if let Some(_kmer) = encode_kmer(&test_seq[i..i + k]) {
                count += 1;
            }
        }
    }
    let elapsed = start.elapsed();
    let per_encode = elapsed.as_nanos() / count as u128;

    println!("  Total time: {:?}", elapsed);
    println!("  Encodings: {}", count);
    println!("  Time per encoding: {} ns", per_encode);
    println!("  Rate: {:.1} M encodings/sec", 1000.0 / per_encode as f64);
    println!();

    // Benchmark 2: ArrayEntropyTracker operations
    println!("Benchmark 2: ArrayEntropyTracker operations");

    let mut tracker = ArrayEntropyTracker::new(k, window);

    // Test add_kmer
    let test_kmer = encode_kmer(b"ACGTA").unwrap();
    let start = Instant::now();
    for _ in 0..1000000 {
        tracker.add_kmer(test_kmer);
        tracker.remove_kmer(test_kmer);  // Keep it balanced
    }
    let elapsed = start.elapsed();

    println!("  add_kmer + remove_kmer (1M operations):");
    println!("    Total time: {:?}", elapsed);
    println!("    Time per op: {} ns", elapsed.as_nanos() / 2000000);
    println!();

    // Test entropy calculation
    tracker = ArrayEntropyTracker::new(k, window);
    for i in 0..20 {
        if let Some(kmer) = encode_kmer(&test_seq[i..i + k]) {
            tracker.add_kmer(kmer);
        }
    }

    let start = Instant::now();
    for _ in 0..1000000 {
        let _e = tracker.entropy();
    }
    let elapsed = start.elapsed();

    println!("  entropy() (1M calls):");
    println!("    Total time: {:?}", elapsed);
    println!("    Time per call: {} ns", elapsed.as_nanos() / 1000000);
    println!();

    // Benchmark 3: Realistic sliding window scenario
    println!("Benchmark 3: Realistic sliding window");
    println!("  Simulating 10Kbp read with window={}, k={}", window, k);

    let mut tracker = ArrayEntropyTracker::new(k, window);
    let start = Instant::now();

    for _ in 0..iterations {
        tracker.clear();

        // Initialize first window
        for j in 0..window - k + 1 {
            if let Some(kmer) = encode_kmer(&test_seq[j..j + k]) {
                tracker.add_kmer(kmer);
            }
        }

        // Slide window through sequence
        for i in window..test_seq.len() {
            // Remove leftmost k-mer
            let exit_pos = i - window;
            if let Some(kmer) = encode_kmer(&test_seq[exit_pos..exit_pos + k]) {
                tracker.remove_kmer(kmer);
            }

            // Add rightmost k-mer
            let enter_pos = i - k + 1;
            if let Some(kmer) = encode_kmer(&test_seq[enter_pos..enter_pos + k]) {
                tracker.add_kmer(kmer);
            }

            // Calculate entropy
            let _e = tracker.entropy();
        }
    }

    let elapsed = start.elapsed();
    let per_read = elapsed.as_millis() as f64 / iterations as f64;
    let windows_per_read = test_seq.len() - window + 1;
    let per_window = elapsed.as_nanos() as f64 / (iterations * windows_per_read) as f64;

    println!("  Total time: {:?}", elapsed);
    println!("  Time per 10Kbp read: {:.3} ms", per_read);
    println!("  Time per window: {:.0} ns", per_window);
    println!("  Throughput: {:.1} reads/sec", 1000.0 / per_read);
    println!();

    // Benchmark 4: Component breakdown
    println!("Benchmark 4: Component time breakdown (single pass)");

    let mut encode_time = 0u128;
    let mut add_time = 0u128;
    let mut remove_time = 0u128;
    let mut entropy_time = 0u128;

    let mut tracker = ArrayEntropyTracker::new(k, window);

    // Initialize first window
    for j in 0..window - k + 1 {
        let start = Instant::now();
        let kmer_opt = encode_kmer(&test_seq[j..j + k]);
        encode_time += start.elapsed().as_nanos();

        if let Some(kmer) = kmer_opt {
            let start = Instant::now();
            tracker.add_kmer(kmer);
            add_time += start.elapsed().as_nanos();
        }
    }

    // Slide through sequence
    for i in window..test_seq.len() {
        // Remove
        let exit_pos = i - window;
        let start = Instant::now();
        let kmer_opt = encode_kmer(&test_seq[exit_pos..exit_pos + k]);
        encode_time += start.elapsed().as_nanos();

        if let Some(kmer) = kmer_opt {
            let start = Instant::now();
            tracker.remove_kmer(kmer);
            remove_time += start.elapsed().as_nanos();
        }

        // Add
        let enter_pos = i - k + 1;
        let start = Instant::now();
        let kmer_opt = encode_kmer(&test_seq[enter_pos..enter_pos + k]);
        encode_time += start.elapsed().as_nanos();

        if let Some(kmer) = kmer_opt {
            let start = Instant::now();
            tracker.add_kmer(kmer);
            add_time += start.elapsed().as_nanos();
        }

        // Entropy
        let start = Instant::now();
        let _e = tracker.entropy();
        entropy_time += start.elapsed().as_nanos();
    }

    let total = encode_time + add_time + remove_time + entropy_time;

    println!("  encode_kmer:  {:8} ns  ({:5.1}%)", encode_time, 100.0 * encode_time as f64 / total as f64);
    println!("  add_kmer:     {:8} ns  ({:5.1}%)", add_time, 100.0 * add_time as f64 / total as f64);
    println!("  remove_kmer:  {:8} ns  ({:5.1}%)", remove_time, 100.0 * remove_time as f64 / total as f64);
    println!("  entropy():    {:8} ns  ({:5.1}%)", entropy_time, 100.0 * entropy_time as f64 / total as f64);
    println!("  Total:        {:8} ns", total);
    println!();

    println!("========================================");
    println!("Summary");
    println!("========================================");
    println!("The component breakdown shows where optimization");
    println!("efforts should be focused:");
    println!();

    if encode_time as f64 / total as f64 > 0.3 {
        println!("⚠️  encode_kmer takes >{:.0}% of time", 100.0 * encode_time as f64 / total as f64);
        println!("   → SIMD optimization could help significantly");
    } else {
        println!("✓  encode_kmer is not the bottleneck (<30%)");
        println!("   → SIMD optimization would have limited impact");
    }
    println!();

    if add_time + remove_time > encode_time {
        println!("⚠️  Tracker operations are slower than encoding");
        println!("   → Focus on optimizing ArrayEntropyTracker");
    } else {
        println!("✓  Tracker operations are efficient");
    }
    println!();

    if entropy_time as f64 / total as f64 > 0.1 {
        println!("⚠️  entropy() takes >{:.0}% despite being O(1)", 100.0 * entropy_time as f64 / total as f64);
        println!("   → Check for floating-point overhead");
    } else {
        println!("✓  entropy() is very fast (O(1) optimization working!)");
    }
    println!("========================================");
}
