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
fn shannon_entropy(kmer_counts: &HashMap<Vec<u8>, usize>, total_kmers: usize, k: usize) -> f64 {
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
    // Maximum entropy occurs when all 4^k possible k-mers are equally likely
    // max_entropy = log2(4^k) = k * log2(4) = k * 2
    let max_entropy = (k as f64) * 2.0;

    if max_entropy > 0.0 {
        entropy / max_entropy
    } else {
        entropy
    }
}

/// Get reverse complement of a DNA sequence
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&base| match base {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            _ => base,
        })
        .collect()
}

/// Get canonical k-mer (lexicographically smaller of k-mer and its reverse complement)
fn canonical_kmer(kmer: &[u8]) -> Vec<u8> {
    let rc = reverse_complement(kmer);
    if kmer <= rc.as_slice() {
        kmer.to_vec()
    } else {
        rc
    }
}

/// Extract all k-mers from a sequence window (using canonical form)
fn get_kmers(sequence: &[u8], k: usize) -> HashMap<Vec<u8>, usize> {
    let mut kmer_counts = HashMap::new();

    if sequence.len() < k {
        return kmer_counts;
    }

    for i in 0..=sequence.len() - k {
        let kmer = &sequence[i..i + k];
        // Only count k-mers with valid bases (ACGTN)
        if kmer.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n')) {
            let canonical = canonical_kmer(kmer);
            *kmer_counts.entry(canonical).or_insert(0) += 1;
        }
    }

    kmer_counts
}

/// Mask low-complexity regions in a sequence based on entropy
fn mask_sequence(sequence: &[u8], quality: &[u8], window: usize, entropy_threshold: f64, k: usize) -> (Vec<u8>, Vec<u8>) {
    let seq_len = sequence.len();
    let mut masked_seq = sequence.to_vec();
    let mut masked_qual = quality.to_vec();

    if seq_len < window {
        // If sequence is shorter than window, calculate entropy for the whole sequence
        let kmer_counts = get_kmers(sequence, k);
        let total_kmers = if seq_len >= k { seq_len - k + 1 } else { 0 };
        let entropy = shannon_entropy(&kmer_counts, total_kmers, k);

        if entropy < entropy_threshold {
            // Mask entire sequence
            for i in 0..seq_len {
                masked_seq[i] = b'N';
                masked_qual[i] = b'#';
            }
        }
        return (masked_seq, masked_qual);
    }

    // For each position, calculate entropy of the surrounding window
    for i in 0..seq_len {
        // Determine window boundaries
        let window_start = if i >= window / 2 {
            i - window / 2
        } else {
            0
        };
        let window_end = (i + window / 2 + 1).min(seq_len);

        // Get the window sequence
        let window_seq = &sequence[window_start..window_end];

        // Calculate entropy for this window
        let kmer_counts = get_kmers(window_seq, k);
        let total_kmers = if window_seq.len() >= k { window_seq.len() - k + 1 } else { 0 };
        let entropy = shannon_entropy(&kmer_counts, total_kmers, k);

        // Mask if entropy is below threshold
        if entropy < entropy_threshold {
            masked_seq[i] = b'N';
            masked_qual[i] = b'#';
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

        let k = 2; // k-mer size
        let entropy = shannon_entropy(&counts, 4, k);
        // 4 equal kmers → raw entropy = log2(4) = 2.0
        // Normalized: 2.0 / (k * 2) = 2.0 / 4.0 = 0.5
        assert!((entropy - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_shannon_entropy_low_complexity() {
        let mut counts = HashMap::new();
        counts.insert(vec![b'A', b'A'], 10);

        let k = 2; // k-mer size
        let entropy = shannon_entropy(&counts, 10, k);
        assert_eq!(entropy, 0.0); // All same kmer = no entropy (still 0 after normalization)
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GCGCG"), b"CGCGC");
        assert_eq!(reverse_complement(b"CGCGC"), b"GCGCG");
    }

    #[test]
    fn test_canonical_kmer() {
        // GCGCG and CGCGC should both canonicalize to CGCGC
        assert_eq!(canonical_kmer(b"GCGCG"), b"CGCGC");
        assert_eq!(canonical_kmer(b"CGCGC"), b"CGCGC");

        // ACG and CGT are reverse complements
        assert_eq!(canonical_kmer(b"ACG"), b"ACG");
        assert_eq!(canonical_kmer(b"CGT"), b"ACG");
    }

    #[test]
    fn test_get_kmers() {
        let sequence = b"ACGTACGT";
        let kmers = get_kmers(sequence, 3);

        // With canonical k-mers:
        // ACG and CGT both map to ACG
        // GTA and TAC both map to GTA
        assert_eq!(kmers.get(&vec![b'A', b'C', b'G']).unwrap(), &4);
        assert_eq!(kmers.get(&vec![b'G', b'T', b'A']).unwrap(), &2);
    }

    #[test]
    fn test_gcgcgc_is_low_complexity() {
        // GCGCGC should be masked because canonical k-mers collapse to one
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
