// Shared library for mask_fastq and mask_fastq_parallel
use std::collections::HashMap;

/// Encode a k-mer into a u16 using 2 bits per base (A=00, C=01, G=10, T=11)
/// Returns None if the k-mer contains N or invalid bases
/// Maximum k-mer size: 8 bases (16 bits / 2 bits per base)
pub fn encode_kmer(bases: &[u8]) -> Option<u16> {
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
pub fn shannon_entropy(kmer_counts: &HashMap<u16, usize>, total_kmers: usize) -> f64 {
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
pub fn get_kmers(sequence: &[u8], k: usize) -> HashMap<u16, usize> {
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
pub fn add_kmer(kmer_counts: &mut HashMap<u16, usize>, kmer: &[u8]) {
    if let Some(encoded) = encode_kmer(kmer) {
        *kmer_counts.entry(encoded).or_insert(0) += 1;
    }
}

/// Remove a k-mer from the counts (used for incremental sliding window)
/// Uses u16 bit-packed encoding for efficient HashMap operations
pub fn remove_kmer(kmer_counts: &mut HashMap<u16, usize>, kmer: &[u8]) {
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
pub fn mask_sequence(sequence: &[u8], quality: &[u8], window: usize, entropy_threshold: f64, k: usize) -> (Vec<u8>, Vec<u8>) {
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
        let sequence = b"GCGCGCGCGCGCGCGCGCGCGCGCGC";
        let quality = vec![b'I'; 26];
        let (masked_seq, _) = mask_sequence(sequence, &quality, 25, 0.55, 5);

        // Should be entirely masked
        let masked_count = masked_seq.iter().filter(|&&b| b == b'N').count();
        assert_eq!(masked_count, 26);
    }

    #[test]
    fn test_mask_low_complexity() {
        // Low complexity: many repeats
        let sequence = b"AAAAAAAAAA";
        let quality = vec![b'I'; 10];
        let (masked_seq, _) = mask_sequence(sequence, &quality, 5, 0.55, 3);

        // Should be entirely masked
        let masked_count = masked_seq.iter().filter(|&&b| b == b'N').count();
        assert_eq!(masked_count, 10);
    }

    #[test]
    fn test_no_mask_high_complexity() {
        // High complexity: random sequence
        let sequence = b"ACGTACGTAGCTAGCT";
        let quality = vec![b'I'; 16];
        let (masked_seq, _) = mask_sequence(sequence, &quality, 5, 0.55, 3);

        // Should not be masked (high entropy)
        let masked_count = masked_seq.iter().filter(|&&b| b == b'N').count();
        assert_eq!(masked_count, 0);
    }
}
