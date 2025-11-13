// Shared library for rustmasker
use std::collections::HashMap;

/// Encode a k-mer into a u32 using 2 bits per base (A=00, C=01, G=10, T=11)
/// Returns None if the k-mer contains N or invalid bases
/// Maximum k-mer size: 15 bases (30 bits / 2 bits per base)
pub fn encode_kmer(bases: &[u8]) -> Option<u32> {
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
            _ => return None,  // N or invalid base - skip this k-mer
        };
        encoded = (encoded << 2) | bits;
    }
    Some(encoded)
}

/// Calculate Shannon entropy from k-mer frequencies
/// Returns normalized entropy in range [0, 1]
pub fn shannon_entropy(kmer_counts: &HashMap<u32, usize>, total_kmers: usize) -> f64 {
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
/// Uses u32 bit-packed encoding for efficient HashMap operations
pub fn get_kmers(sequence: &[u8], k: usize) -> HashMap<u32, usize> {
    let mut kmer_counts = HashMap::new();

    if sequence.len() < k {
        return kmer_counts;
    }

    for i in 0..=sequence.len() - k {
        let kmer = &sequence[i..i + k];
        // Encode k-mer as u32; skip if contains N or invalid bases
        if let Some(encoded) = encode_kmer(kmer) {
            *kmer_counts.entry(encoded).or_insert(0) += 1;
        }
    }

    kmer_counts
}

/// Add a k-mer to the counts (used for incremental sliding window)
/// Uses u32 bit-packed encoding for efficient HashMap operations
pub fn add_kmer(kmer_counts: &mut HashMap<u32, usize>, kmer: &[u8]) {
    if let Some(encoded) = encode_kmer(kmer) {
        *kmer_counts.entry(encoded).or_insert(0) += 1;
    }
}

/// Remove a k-mer from the counts (used for incremental sliding window)
/// Uses u32 bit-packed encoding for efficient HashMap operations
pub fn remove_kmer(kmer_counts: &mut HashMap<u32, usize>, kmer: &[u8]) {
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
    // Use incremental k-mer tracking with u32 bit-packed keys for optimal performance

    let mut kmer_counts: HashMap<u32, usize> = HashMap::new();
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

// ============================================================================
// Array-Based Entropy Tracker (BBMask-inspired optimization)
// ============================================================================

/// Array-based entropy tracker for efficient O(1) entropy calculations
/// Based on BBTools EntropyTracker design:
/// - Uses fixed-size arrays for k-mer counts (k ≤ 7 recommended)
/// - Maintains count-of-counts histogram for O(1) entropy updates
/// - Precalculates entropy values to avoid log() in hot path
pub struct ArrayEntropyTracker {
    window_kmers: usize,
    counts: Vec<u16>,           // K-mer counts (size 4^k)
    count_counts: Vec<u16>,     // Histogram of count frequencies (size window_kmers+2)
    entropy_table: Vec<f64>,    // Precalculated p*log2(p) for each possible count
    entropy_mult: f64,          // Normalization factor: -1/log2(window_kmers)
    current_esum: f64,          // Running entropy sum
    unique: usize,              // Number of unique k-mers
}

impl ArrayEntropyTracker {
    /// Create a new array-based entropy tracker
    ///
    /// # Arguments
    /// * `k` - K-mer size (recommended: k ≤ 7 for reasonable memory usage)
    /// * `window` - Window size in bases
    ///
    /// # Memory usage
    /// - k=5: ~4 KB (1024 kmers × 2 bytes)
    /// - k=6: ~16 KB (4096 kmers × 2 bytes)
    /// - k=7: ~64 KB (16384 kmers × 2 bytes)
    /// - k=8: ~256 KB (65536 kmers × 2 bytes)
    /// - k > 8: Not recommended (use HashMap instead for k > 7)
    pub fn new(k: usize, window: usize) -> Self {
        assert!(k > 0 && k <= 15, "k must be in range 1-15");
        assert!(window > k, "window must be larger than k");

        let window_kmers = window - k + 1;
        let kmer_space = 1 << (2 * k); // 4^k

        // Precalculate entropy table: entropy[count] = (count/window_kmers) * log2(count/window_kmers)
        // This matches BBMask's approach
        let mut entropy_table = vec![0.0; window_kmers + 2];
        for count in 1..entropy_table.len() {
            let p = count as f64 / window_kmers as f64;
            entropy_table[count] = p * p.log2();
        }

        // Normalization factor to convert entropy to 0-1 scale
        let entropy_mult = -1.0 / (window_kmers as f64).log2();

        // Initialize count_counts with all kmers having count 0
        let mut count_counts = vec![0u16; window_kmers + 2];
        count_counts[0] = window_kmers as u16;

        Self {
            window_kmers,
            counts: vec![0; kmer_space],
            count_counts,
            entropy_table,
            entropy_mult,
            current_esum: 0.0,
            unique: 0,
        }
    }

    /// Add a k-mer to the tracker (incremental window update)
    /// Updates counts, count_counts histogram, and running entropy sum
    /// Time complexity: O(1)
    pub fn add_kmer(&mut self, kmer_code: u32) {
        let old_count = self.counts[kmer_code as usize];
        let new_count = old_count + 1;

        // Update unique count
        if old_count == 0 {
            self.unique += 1;
        }

        // Update count_counts histogram
        self.count_counts[old_count as usize] -= 1;
        self.count_counts[new_count as usize] += 1;

        // Update k-mer count
        self.counts[kmer_code as usize] = new_count;

        // Update running entropy sum using precalculated values
        self.current_esum += self.entropy_table[new_count as usize]
                           - self.entropy_table[old_count as usize];
    }

    /// Remove a k-mer from the tracker (incremental window update)
    /// Updates counts, count_counts histogram, and running entropy sum
    /// Time complexity: O(1)
    pub fn remove_kmer(&mut self, kmer_code: u32) {
        let old_count = self.counts[kmer_code as usize];
        if old_count == 0 {
            return; // Nothing to remove
        }

        let new_count = old_count - 1;

        // Update count_counts histogram
        self.count_counts[old_count as usize] -= 1;
        self.count_counts[new_count as usize] += 1;

        // Update k-mer count
        self.counts[kmer_code as usize] = new_count;

        // Update running entropy sum using precalculated values
        self.current_esum += self.entropy_table[new_count as usize]
                           - self.entropy_table[old_count as usize];

        // Update unique count
        if new_count == 0 {
            self.unique -= 1;
        }
    }

    /// Get current entropy (normalized to 0-1 scale)
    /// Time complexity: O(1) - just returns cached value!
    #[inline]
    pub fn entropy(&self) -> f64 {
        let e = self.current_esum * self.entropy_mult;
        // Avoid negative zero due to floating point errors (branchless)
        e.max(0.0)
    }

    /// Clear the tracker for a new sequence
    pub fn clear(&mut self) {
        // Reset counts array
        for count in &mut self.counts {
            *count = 0;
        }

        // Reset count_counts (all kmers start with count 0)
        for cc in &mut self.count_counts {
            *cc = 0;
        }
        self.count_counts[0] = self.window_kmers as u16;

        // Reset accumulators
        self.current_esum = 0.0;
        self.unique = 0;
    }

    /// Get number of unique k-mers in current window
    pub fn unique(&self) -> usize {
        self.unique
    }
}

/// Mask low-complexity regions using array-based entropy tracker
/// Optimized version of mask_sequence() that uses O(1) entropy calculations
/// Recommended for k ≤ 7 (larger k uses more memory but still works)
pub fn mask_sequence_array(
    sequence: &[u8],
    quality: &[u8],
    window: usize,
    entropy_threshold: f64,
    k: usize
) -> (Vec<u8>, Vec<u8>) {
    let seq_len = sequence.len();
    let mut masked_seq = sequence.to_vec();
    let mut masked_qual = quality.to_vec();

    if seq_len < window {
        // If sequence is shorter than window, calculate entropy for the whole sequence
        // Fall back to HashMap for short sequences (not worth the array overhead)
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

    // Use array-based tracker for sliding window
    let mut tracker = ArrayEntropyTracker::new(k, window);
    let mut first_full_window = true;

    for i in 0..seq_len {
        // Window extends from [window_start, window_end)
        let window_start = if i + 1 >= window {
            i + 1 - window
        } else {
            0
        };
        let window_end = i + 1;

        // Only check entropy once window is full
        if window_end - window_start < window {
            continue;
        }

        if first_full_window {
            // First full window: initialize k-mer counts
            tracker.clear();
            for j in window_start..=window_end.saturating_sub(k) {
                if let Some(kmer_code) = encode_kmer(&sequence[j..j + k]) {
                    tracker.add_kmer(kmer_code);
                }
            }
            first_full_window = false;
        } else {
            // Subsequent windows: slide forward by 1 base
            // Remove the leftmost k-mer that just exited
            let exiting_kmer_pos = window_start - 1;
            if exiting_kmer_pos + k <= seq_len {
                if let Some(kmer_code) = encode_kmer(&sequence[exiting_kmer_pos..exiting_kmer_pos + k]) {
                    tracker.remove_kmer(kmer_code);
                }
            }

            // Add the new rightmost k-mer that just entered
            let entering_kmer_pos = window_end - k;
            if entering_kmer_pos < seq_len && entering_kmer_pos + k <= seq_len {
                if let Some(kmer_code) = encode_kmer(&sequence[entering_kmer_pos..entering_kmer_pos + k]) {
                    tracker.add_kmer(kmer_code);
                }
            }
        }

        // Get entropy - O(1) operation!
        let entropy = tracker.entropy();

        // If entropy is below threshold, mask the entire window range
        if entropy < entropy_threshold {
            for pos in window_start..window_end {
                masked_seq[pos] = b'N';
                masked_qual[pos] = b'#';
            }
        }
    }

    (masked_seq, masked_qual)
}

/// Automatically choose between array-based and HashMap-based masking based on k
/// - Uses array-based for k <= 7 (memory: 4KB for k=5, 16KB for k=6, 64KB for k=7)
/// - Uses HashMap-based for k > 7 (to avoid excessive memory usage)
///
/// This provides the best performance for typical k values while gracefully
/// handling larger k values that would require too much memory for arrays.
pub fn mask_sequence_auto(
    sequence: &[u8],
    quality: &[u8],
    window: usize,
    entropy_threshold: f64,
    k: usize
) -> (Vec<u8>, Vec<u8>) {
    if k <= 7 {
        // Use optimized array-based implementation (1.7-3.2x faster)
        mask_sequence_array(sequence, quality, window, entropy_threshold, k)
    } else {
        // Fall back to HashMap for k > 7 to avoid excessive memory (256KB+ for k=8)
        mask_sequence(sequence, quality, window, entropy_threshold, k)
    }
}

// ============================================================================
// SDUST Algorithm Implementation
// ============================================================================

/// Encode a triplet (3 bases) into a 6-bit word (0-63)
/// A=0b00, C=0b01, G=0b10, T=0b11
/// Returns None if sequence contains N or invalid bases
pub fn encode_triplet(bases: &[u8]) -> Option<u8> {
    if bases.len() != 3 {
        return None;
    }

    let mut word: u8 = 0;
    for &base in bases {
        let bits = match base {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => return None,  // N or invalid base
        };
        word = (word << 2) | bits;
    }
    Some(word)
}

/// Represents a region to be masked
#[derive(Debug, Clone, Copy)]
struct MaskRegion {
    start: usize,  // Start position (inclusive)
    end: usize,    // End position (exclusive)
}

/// SDUST window scorer
struct SdustScorer {
    threshold: i32,
}

impl SdustScorer {
    pub fn new(threshold: i32) -> Self {
        Self { threshold }
    }

    /// Calculate DUST score for a window of triplets
    /// DUST score = sum of count*(count-1)/2 for each unique triplet
    /// Returns the DUST score
    fn score_window(&self, triplets: &[u8]) -> i32 {
        let mut counts = [0u16; 64];  // 64 possible triplet values

        for &triplet in triplets {
            counts[triplet as usize] += 1;
        }

        // Calculate DUST score: sum of count*(count-1)/2 for each triplet
        let mut score = 0i32;
        for count in counts.iter() {
            let c = *count as i32;
            if c > 1 {
                score += c * (c - 1) / 2;
            }
        }
        score
    }

    /// Check if a window should be masked
    /// Formula from sdust: score * 10 / window_length > threshold
    fn should_mask(&self, score: i32, window_length: usize) -> bool {
        if window_length == 0 {
            return false;
        }
        score * 10 > self.threshold * window_length as i32
    }
}

/// Find regions to mask using SDUST algorithm
fn find_dust_regions(
    triplets: &[u8],
    positions: &[usize],
    window_size: usize,
    threshold: i32,
) -> Vec<MaskRegion> {
    let scorer = SdustScorer::new(threshold);
    let mut regions = Vec::new();

    if triplets.len() < window_size {
        // Score entire sequence
        let score = scorer.score_window(triplets);
        if scorer.should_mask(score, triplets.len()) {
            if let (Some(&start), Some(&end)) = (positions.first(), positions.last()) {
                regions.push(MaskRegion { start, end: end + 3 });
            }
        }
        return regions;
    }

    // Slide window over triplets
    let mut current_region: Option<MaskRegion> = None;

    for i in 0..=triplets.len().saturating_sub(window_size) {
        let window = &triplets[i..i + window_size];
        let score = scorer.score_window(window);

        if scorer.should_mask(score, window.len()) {
            let start_pos = positions[i];
            let end_pos = positions[i + window_size - 1] + 3;

            match current_region.as_mut() {
                Some(region) if region.end >= start_pos => {
                    // Extend existing region
                    region.end = region.end.max(end_pos);
                }
                _ => {
                    // Start new region
                    if let Some(region) = current_region.take() {
                        regions.push(region);
                    }
                    current_region = Some(MaskRegion { start: start_pos, end: end_pos });
                }
            }
        }
    }

    if let Some(region) = current_region {
        regions.push(region);
    }

    regions
}

/// Apply masks to sequence and quality scores
fn apply_masks(seq: &mut [u8], qual: &mut [u8], regions: &[MaskRegion]) {
    for region in regions {
        for i in region.start..region.end.min(seq.len()) {
            seq[i] = b'N';
            qual[i] = b'#';
        }
    }
}

/// Mask low-complexity regions using the SDUST algorithm
///
/// # Arguments
/// * `sequence` - DNA sequence bytes (A/C/G/T/N)
/// * `quality` - Quality scores (same length as sequence)
/// * `window_size` - Window size W (default: 64)
/// * `threshold` - Score threshold T (default: 20)
///
/// # Returns
/// Tuple of (masked_sequence, masked_quality) where low-complexity regions
/// are replaced with 'N' (sequence) and '#' (quality)
pub fn mask_sequence_sdust(
    sequence: &[u8],
    quality: &[u8],
    window_size: usize,
    threshold: i32,
) -> (Vec<u8>, Vec<u8>) {
    let mut masked_seq = sequence.to_vec();
    let mut masked_qual = quality.to_vec();

    if sequence.len() < 3 {
        return (masked_seq, masked_qual);
    }

    // Convert sequence to triplets
    let mut triplets = Vec::new();
    let mut triplet_positions = Vec::new();  // Track where each triplet starts

    for i in 0..sequence.len().saturating_sub(2) {
        if let Some(triplet) = encode_triplet(&sequence[i..i+3]) {
            triplets.push(triplet);
            triplet_positions.push(i);
        } else {
            // N or invalid base - process accumulated triplets
            if !triplets.is_empty() {
                let regions = find_dust_regions(&triplets, &triplet_positions, window_size, threshold);
                apply_masks(&mut masked_seq, &mut masked_qual, &regions);
                triplets.clear();
                triplet_positions.clear();
            }
        }
    }

    // Process final batch
    if !triplets.is_empty() {
        let regions = find_dust_regions(&triplets, &triplet_positions, window_size, threshold);
        apply_masks(&mut masked_seq, &mut masked_qual, &regions);
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

    // Tests for ArrayEntropyTracker

    #[test]
    fn test_array_tracker_basic() {
        let mut tracker = ArrayEntropyTracker::new(3, 10);

        // Add some k-mers
        let kmer_aaa = encode_kmer(b"AAA").unwrap();
        let kmer_ccc = encode_kmer(b"CCC").unwrap();

        tracker.add_kmer(kmer_aaa);
        assert_eq!(tracker.unique(), 1);

        tracker.add_kmer(kmer_ccc);
        assert_eq!(tracker.unique(), 2);

        tracker.add_kmer(kmer_aaa); // Add AAA again
        assert_eq!(tracker.unique(), 2); // Still 2 unique

        // Remove a k-mer
        tracker.remove_kmer(kmer_aaa);
        assert_eq!(tracker.unique(), 2); // Still 2 unique (AAA count is 1)

        tracker.remove_kmer(kmer_aaa);
        assert_eq!(tracker.unique(), 1); // Now only CCC remains
    }

    #[test]
    fn test_array_tracker_entropy() {
        let mut tracker = ArrayEntropyTracker::new(2, 10);
        // window_kmers = 10 - 2 + 1 = 9

        // All same k-mer = low entropy
        let kmer_aa = encode_kmer(b"AA").unwrap();
        for _ in 0..9 {
            tracker.add_kmer(kmer_aa);
        }
        let entropy = tracker.entropy();
        assert!(entropy < 0.01, "All same k-mer should have ~0 entropy, got {}", entropy);

        // Add diverse k-mers (fill the window with unique k-mers) = higher entropy
        tracker.clear();
        let kmers = [
            encode_kmer(b"AA").unwrap(),
            encode_kmer(b"AC").unwrap(),
            encode_kmer(b"AG").unwrap(),
            encode_kmer(b"AT").unwrap(),
            encode_kmer(b"CA").unwrap(),
            encode_kmer(b"CC").unwrap(),
            encode_kmer(b"CG").unwrap(),
            encode_kmer(b"CT").unwrap(),
            encode_kmer(b"GA").unwrap(),
        ];

        // Add 9 unique k-mers (fills the window)
        for &kmer in &kmers {
            tracker.add_kmer(kmer);
        }

        let entropy = tracker.entropy();
        // All unique k-mers = maximum entropy = 1.0
        assert!(entropy > 0.99, "All unique k-mers should have entropy ~1.0, got {}", entropy);
    }

    #[test]
    fn test_mask_sequence_array_matches_hashmap() {
        // Test that array-based and HashMap-based implementations produce identical results

        let test_cases = vec![
            (b"AAAAAAAAAA".as_ref(), "homopolymer"),
            (b"GCGCGCGCGCGCGCGCGCGCGCGCGC".as_ref(), "dinucleotide repeat"),
            (b"ACGTACGTAGCTAGCT".as_ref(), "random sequence"),
            (b"ATATATATATATATATAT".as_ref(), "AT repeat"),
        ];

        for (sequence, description) in test_cases {
            let quality = vec![b'I'; sequence.len()];

            let (masked_hashmap, qual_hashmap) = mask_sequence(sequence, &quality, 25, 0.55, 5);
            let (masked_array, qual_array) = mask_sequence_array(sequence, &quality, 25, 0.55, 5);

            assert_eq!(
                masked_hashmap, masked_array,
                "Sequences don't match for {}: \nHashMap: {}\nArray: {}",
                description,
                String::from_utf8_lossy(&masked_hashmap),
                String::from_utf8_lossy(&masked_array)
            );

            assert_eq!(
                qual_hashmap, qual_array,
                "Quality scores don't match for {}",
                description
            );
        }
    }

    #[test]
    fn test_mask_sequence_array_low_complexity() {
        let sequence = b"AAAAAAAAAA";
        let quality = vec![b'I'; 10];
        let (masked_seq, _) = mask_sequence_array(sequence, &quality, 5, 0.55, 3);

        // Should be entirely masked
        let masked_count = masked_seq.iter().filter(|&&b| b == b'N').count();
        assert_eq!(masked_count, 10);
    }

    #[test]
    fn test_mask_sequence_array_high_complexity() {
        let sequence = b"ACGTACGTAGCTAGCT";
        let quality = vec![b'I'; 16];
        let (masked_seq, _) = mask_sequence_array(sequence, &quality, 5, 0.55, 3);

        // Should not be masked
        let masked_count = masked_seq.iter().filter(|&&b| b == b'N').count();
        assert_eq!(masked_count, 0);
    }

    #[test]
    fn test_mask_sequence_array_gcgc() {
        let sequence = b"GCGCGCGCGCGCGCGCGCGCGCGCGC";
        let quality = vec![b'I'; 26];
        let (masked_seq, _) = mask_sequence_array(sequence, &quality, 25, 0.55, 5);

        // Should be entirely masked
        let masked_count = masked_seq.iter().filter(|&&b| b == b'N').count();
        assert_eq!(masked_count, 26);
    }

    // Tests for SDUST algorithm

    #[test]
    fn test_encode_triplet() {
        // Test valid triplets (A=00, C=01, G=10, T=11)
        assert_eq!(encode_triplet(b"AAA"), Some(0b000000));  // 0
        assert_eq!(encode_triplet(b"AAC"), Some(0b000001));  // 1
        assert_eq!(encode_triplet(b"AAG"), Some(0b000010));  // 2
        assert_eq!(encode_triplet(b"AAT"), Some(0b000011));  // 3
        assert_eq!(encode_triplet(b"CCC"), Some(0b010101));  // 21
        assert_eq!(encode_triplet(b"GGG"), Some(0b101010));  // 42
        assert_eq!(encode_triplet(b"TTT"), Some(0b111111));  // 63
        assert_eq!(encode_triplet(b"ACG"), Some(0b000110));  // 6

        // Test lowercase
        assert_eq!(encode_triplet(b"aaa"), Some(0b000000));
        assert_eq!(encode_triplet(b"acg"), Some(0b000110));

        // Test invalid cases
        assert_eq!(encode_triplet(b"AAN"), None);
        assert_eq!(encode_triplet(b"NNN"), None);
        assert_eq!(encode_triplet(b"AA"), None);
        assert_eq!(encode_triplet(b"AAAA"), None);
    }

    #[test]
    fn test_sdust_homopolymer() {
        // Homopolymers should be masked with default parameters
        let seq = b"AAAAAAAAAAAAAAAA";  // 16 A's
        let qual = vec![b'I'; 16];
        let (masked, masked_qual) = mask_sequence_sdust(seq, &qual, 10, 20);

        // Homopolymer should be masked
        assert!(masked.iter().all(|&b| b == b'N'), "Expected all bases to be masked");
        assert!(masked_qual.iter().all(|&b| b == b'#'), "Expected all quality to be masked");
    }

    #[test]
    fn test_sdust_high_complexity() {
        // High complexity sequence should not be masked
        let seq = b"ACGTACGTACGTACGT";  // High complexity
        let qual = vec![b'I'; 16];
        let (masked, masked_qual) = mask_sequence_sdust(seq, &qual, 10, 20);

        // Should not be masked
        assert_eq!(masked, seq);
        assert_eq!(masked_qual, qual);
    }

    #[test]
    fn test_sdust_with_n() {
        // N bases should break the sequence into segments
        // Use longer homopolymer runs and lower threshold
        let seq = b"AAAAAAAAANNNNGGGGGGGGGG";
        let qual = vec![b'I'; 23];
        let (masked, _) = mask_sequence_sdust(seq, &qual, 8, 20);

        // A's should be masked (homopolymer)
        let a_masked = masked[0..9].iter().filter(|&&b| b == b'N').count();
        assert!(a_masked > 5, "Expected most A's to be masked, got {}", a_masked);

        // N's stay as N (positions 9-12)
        assert!(masked[9..13].iter().all(|&b| b == b'N'));

        // G's should be masked (homopolymer)
        let g_masked = masked[13..23].iter().filter(|&&b| b == b'N').count();
        assert!(g_masked > 5, "Expected most G's to be masked, got {}", g_masked);
    }

    #[test]
    fn test_sdust_short_sequence() {
        // Very short sequences (< 3 bases) should be returned unchanged
        let seq = b"AA";
        let qual = vec![b'I'; 2];
        let (masked, masked_qual) = mask_sequence_sdust(seq, &qual, 10, 20);

        assert_eq!(masked, seq);
        assert_eq!(masked_qual, qual);
    }

    #[test]
    fn test_sdust_dinucleotide_repeat() {
        // Dinucleotide repeat (GCGC...) is low complexity and should be masked
        // GCGC produces alternating triplets: GCG and CGC
        let seq = b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";  // 32 bases
        let qual = vec![b'I'; 32];
        // Use window_size that captures the repetition pattern
        let (masked, _) = mask_sequence_sdust(seq, &qual, 16, 20);

        // Should be masked (or mostly masked)
        let masked_count = masked.iter().filter(|&&b| b == b'N').count();
        assert!(masked_count > 16, "Expected most of the dinucleotide repeat to be masked, got {} out of {}", masked_count, 32);
    }

    #[test]
    fn test_sdust_threshold_sensitivity() {
        // Test that higher threshold means less masking
        let seq = b"ATATATATATATAT";  // AT repeat
        let qual = vec![b'I'; 14];

        // Low threshold - should mask more
        let (masked_low, _) = mask_sequence_sdust(seq, &qual, 10, 10);
        let masked_count_low = masked_low.iter().filter(|&&b| b == b'N').count();

        // High threshold - should mask less
        let (masked_high, _) = mask_sequence_sdust(seq, &qual, 10, 50);
        let masked_count_high = masked_high.iter().filter(|&&b| b == b'N').count();

        assert!(masked_count_low >= masked_count_high,
            "Lower threshold should mask at least as much as higher threshold. Low: {}, High: {}",
            masked_count_low, masked_count_high);
    }
}
