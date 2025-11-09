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

    /// K-mer size for entropy calculation
    #[arg(short = 'k', long, default_value_t = 5)]
    kmer: usize,
}

/// Calculate Shannon entropy from k-mer frequencies
/// Returns normalized entropy in range [0, 1]
fn shannon_entropy(kmer_counts: &HashMap<Vec<u8>, usize>, total_kmers: usize) -> f64 {
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
fn get_kmers(sequence: &[u8], k: usize) -> HashMap<Vec<u8>, usize> {
    let mut kmer_counts = HashMap::new();

    if sequence.len() < k {
        return kmer_counts;
    }

    for i in 0..=sequence.len() - k {
        let kmer = &sequence[i..i + k];
        // Only count k-mers with valid bases (ACGTN)
        if kmer.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n')) {
            *kmer_counts.entry(kmer.to_vec()).or_insert(0) += 1;
        }
    }

    kmer_counts
}

/// Check if a k-mer contains only valid bases
fn is_valid_kmer(kmer: &[u8]) -> bool {
    kmer.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n'))
}

/// Add a k-mer to the counts (used for incremental sliding window)
fn add_kmer(kmer_counts: &mut HashMap<Vec<u8>, usize>, kmer: &[u8]) {
    if is_valid_kmer(kmer) {
        *kmer_counts.entry(kmer.to_vec()).or_insert(0) += 1;
    }
}

/// Remove a k-mer from the counts (used for incremental sliding window)
fn remove_kmer(kmer_counts: &mut HashMap<Vec<u8>, usize>, kmer: &[u8]) {
    if is_valid_kmer(kmer) {
        if let Some(count) = kmer_counts.get_mut(&kmer.to_vec()) {
            *count -= 1;
            if *count == 0 {
                kmer_counts.remove(&kmer.to_vec());
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
    // Use incremental k-mer tracking to avoid recalculating all k-mers for overlapping windows

    let mut kmer_counts: HashMap<Vec<u8>, usize> = HashMap::new();
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

    // Create gzip encoder for stdout
    let stdout = io::stdout();
    let gz_writer = GzEncoder::new(stdout, Compression::default());
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
        counts.insert(vec![b'A', b'A'], 1);
        counts.insert(vec![b'C', b'C'], 1);
        counts.insert(vec![b'G', b'G'], 1);
        counts.insert(vec![b'T', b'T'], 1);

        let entropy = shannon_entropy(&counts, 4);
        // 4 unique k-mers out of 4 total → raw entropy = log2(4) = 2.0
        // Normalized by log2(4) = 2.0 → result = 1.0 (perfect entropy!)
        assert!((entropy - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_shannon_entropy_low_complexity() {
        let mut counts = HashMap::new();
        counts.insert(vec![b'A', b'A'], 10);

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
        assert_eq!(kmers.get(&vec![b'A', b'C', b'G']).unwrap(), &2);
        assert_eq!(kmers.get(&vec![b'C', b'G', b'T']).unwrap(), &2);
        assert_eq!(kmers.get(&vec![b'G', b'T', b'A']).unwrap(), &1);
        assert_eq!(kmers.get(&vec![b'T', b'A', b'C']).unwrap(), &1);
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
