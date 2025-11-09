use std::io::{self, BufWriter, Write};
use std::collections::HashMap;
use needletail::parse_fastx_stdin;
use flate2::{Compression, write::GzEncoder};
use clap::Parser;

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

    /// Gzip compression level (0-9, where 0=no compression, 6=default, 9=max compression)
    #[arg(short = 'c', long, default_value_t = 6)]
    compression_level: u32,
}

/// Encode a k-mer into a u16 using 2 bits per base (A=00, C=01, G=10, T=11)
/// Returns None if the k-mer contains N or invalid bases
/// Maximum k-mer size: 8 bases (16 bits / 2 bits per base)
fn encode_kmer(bases: &[u8]) -> Option<u16> {
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
            _ => return None,  // N or invalid base - skip this k-mer
        };
        encoded = (encoded << 2) | bits;
    }
    Some(encoded)
}

/// Calculate Shannon entropy from k-mer frequencies
/// Returns normalized entropy in range [0, 1]
fn shannon_entropy(kmer_counts: &HashMap<u16, usize>, total_kmers: usize) -> f64 {
    if total_kmers == 0 {
        return 0.0;
    }

    let mut entropy = 0.0;
    for &count in kmer_counts.values() {
        if count > 0 {
            let p = count as f64 / total_kmers as f64;
            entropy -= p * p.log2();
        }
    }

    // Normalize entropy to [0, 1] by dividing by maximum possible entropy
    // Maximum entropy occurs when all k-mers in the window are unique
    // For n k-mers in window: max_entropy = log2(n)
    // This matches BBMask's normalization which considers the actual window size
    let max_entropy = (total_kmers as f64).log2();

    if max_entropy > 0.0 {
        entropy / max_entropy
    } else {
        entropy
    }
}

/// Extract all k-mers from a sequence window (strand-specific, no canonicalization)
/// Matches BBMask behavior: counts k-mers as they appear in the sequence
/// Uses u16 bit-packed encoding for efficient HashMap operations
fn get_kmers(sequence: &[u8], k: usize) -> HashMap<u16, usize> {
    let mut kmer_counts = HashMap::new();

    if sequence.len() < k {
        return kmer_counts;
    }

    for i in 0..=sequence.len() - k {
        let kmer = &sequence[i..i + k];
        // Encode k-mer as u16; skip if contains N or invalid bases
        if let Some(encoded) = encode_kmer(kmer) {
            *kmer_counts.entry(encoded).or_insert(0) += 1;
        }
    }

    kmer_counts
}

/// Add a k-mer to the counts (used for incremental sliding window)
/// Uses u16 bit-packed encoding for efficient HashMap operations
fn add_kmer(kmer_counts: &mut HashMap<u16, usize>, kmer: &[u8]) {
    if let Some(encoded) = encode_kmer(kmer) {
        *kmer_counts.entry(encoded).or_insert(0) += 1;
    }
}

/// Remove a k-mer from the counts (used for incremental sliding window)
/// Uses u16 bit-packed encoding for efficient HashMap operations
fn remove_kmer(kmer_counts: &mut HashMap<u16, usize>, kmer: &[u8]) {
    if let Some(encoded) = encode_kmer(kmer) {
        if let Some(count) = kmer_counts.get_mut(&encoded) {
            *count -= 1;
            if *count == 0 {
                kmer_counts.remove(&encoded);
            }
        }
    }
}

/// Mask low-complexity regions in a sequence based on entropy
/// Matches BBMask behavior: masks entire window ranges when low entropy is detected
fn mask_sequence(sequence: &[u8], quality: &[u8], window: usize, entropy_threshold: f64, k: usize) -> (Vec<u8>, Vec<u8>) {
    let seq_len = sequence.len();
    let mut masked_seq = sequence.to_vec();
    let mut masked_qual = quality.to_vec();

    if seq_len < window {
        // If sequence is shorter than window, calculate entropy for the whole sequence
        let kmer_counts = get_kmers(sequence, k);
        let total_kmers = if seq_len >= k { seq_len - k + 1 } else { 0 };
        let entropy = shannon_entropy(&kmer_counts, total_kmers);

        if entropy < entropy_threshold {
            // Mask entire sequence
            for i in 0..seq_len {
                masked_seq[i] = b'N';
                masked_qual[i] = b'#';
            }
        }
        return (masked_seq, masked_qual);
    }

    // BBMask-style sliding window: mask entire window range when low entropy detected
    // Slide window forward one position at a time, checking entropy at each position
    // Use incremental k-mer tracking with u16 bit-packed keys for optimal performance

    let mut kmer_counts: HashMap<u16, usize> = HashMap::new();
    let mut first_full_window = true;

    for i in 0..seq_len {
        // Window extends from [window_start, window_end)
        // Build window up to position i (inclusive)
        let window_start = if i + 1 >= window {
            i + 1 - window
        } else {
            0
        };
        let window_end = i + 1;

        // Only check entropy once window is full (has reached target size)
        if window_end - window_start < window {
            continue;
        }

        if first_full_window {
            // First full window: initialize k-mer counts from scratch
            kmer_counts.clear();
            for j in window_start..=window_end.saturating_sub(k) {
                add_kmer(&mut kmer_counts, &sequence[j..j + k]);
            }
            first_full_window = false;
        } else {
            // Subsequent windows slide forward by 1 base
            // Remove the leftmost k-mer that just exited the window
            let exiting_kmer_pos = window_start - 1;
            if exiting_kmer_pos + k <= seq_len {
                remove_kmer(&mut kmer_counts, &sequence[exiting_kmer_pos..exiting_kmer_pos + k]);
            }

            // Add the new rightmost k-mer that just entered the window
            let entering_kmer_pos = window_end - k;
            if entering_kmer_pos < seq_len && entering_kmer_pos + k <= seq_len {
                add_kmer(&mut kmer_counts, &sequence[entering_kmer_pos..entering_kmer_pos + k]);
            }
        }

        // Calculate entropy for this window
        let total_kmers = if window >= k { window - k + 1 } else { 0 };
        let entropy = shannon_entropy(&kmer_counts, total_kmers);

        // If entropy is below threshold, mask the entire window range
        // This matches BBMask's behavior of masking complete windows
        if entropy < entropy_threshold {
            for pos in window_start..window_end {
                masked_seq[pos] = b'N';
                masked_qual[pos] = b'#';
            }
        }
    }

    (masked_seq, masked_qual)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shannon_entropy_uniform() {
        let mut counts = HashMap::new();
        counts.insert(encode_kmer(b"AA").unwrap(), 1);
        counts.insert(encode_kmer(b"CC").unwrap(), 1);
        counts.insert(encode_kmer(b"GG").unwrap(), 1);
        counts.insert(encode_kmer(b"TT").unwrap(), 1);

        let entropy = shannon_entropy(&counts, 4);
        // 4 unique k-mers out of 4 total → raw entropy = log2(4) = 2.0
        // Normalized by log2(4) = 2.0 → result = 1.0 (perfect entropy!)
        assert!((entropy - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_shannon_entropy_low_complexity() {
        let mut counts = HashMap::new();
        counts.insert(encode_kmer(b"AA").unwrap(), 10);

        let entropy = shannon_entropy(&counts, 10);
        assert_eq!(entropy, 0.0); // All same kmer = no entropy (still 0 after normalization)
    }

    #[test]
    fn test_get_kmers() {
        let sequence = b"ACGTACGT";
        let kmers = get_kmers(sequence, 3);

        // Without canonical k-mers (strand-specific):
        // ACG appears at positions 0 and 4
        // CGT appears at positions 1 and 5
        // GTA appears at position 2
        // TAC appears at position 3
        assert_eq!(kmers.get(&encode_kmer(b"ACG").unwrap()).unwrap(), &2);
        assert_eq!(kmers.get(&encode_kmer(b"CGT").unwrap()).unwrap(), &2);
        assert_eq!(kmers.get(&encode_kmer(b"GTA").unwrap()).unwrap(), &1);
        assert_eq!(kmers.get(&encode_kmer(b"TAC").unwrap()).unwrap(), &1);
    }

    #[test]
    fn test_gcgcgc_is_low_complexity() {
        // GCGCGC should be masked: only 2 distinct k-mers (GCGCG and CGCGC) in 26 total
        // Entropy = log2(2) / log2(26) ≈ 0.21 < 0.55 threshold
        let sequence = b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"; // 30 bases
        let quality = b"JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ";

        let (masked_seq, masked_qual) = mask_sequence(sequence, quality, 25, 0.55, 5);

        // Should be fully masked
        assert_eq!(masked_seq, b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
        assert_eq!(masked_qual, b"##############################");
    }

    #[test]
    fn test_mask_low_complexity() {
        let sequence = b"AAAAAAAA";
        let quality = b"IIIIIIII";

        let (masked_seq, masked_qual) = mask_sequence(sequence, quality, 5, 0.5, 3);

        assert_eq!(masked_seq, b"NNNNNNNN");
        assert_eq!(masked_qual, b"########");
    }

    #[test]
    fn test_no_mask_high_complexity() {
        // Use a longer random sequence with very high k-mer diversity
        let sequence = b"GACTGCATCGTAGCTGATCGACTGCAGTCGATCGACTGCAT";
        let quality = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";

        // With normalized entropy and low threshold, high-complexity sequence should not be masked
        let (masked_seq, masked_qual) = mask_sequence(sequence, quality, 10, 0.1, 3);

        // High complexity sequence should not be masked
        assert_eq!(masked_seq, sequence);
        assert_eq!(masked_qual, quality);
    }
}
