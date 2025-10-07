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
    entropy
}

/// Extract all k-mers from a sequence window
fn get_kmers(sequence: &[u8], k: usize) -> HashMap<Vec<u8>, usize> {
    let mut kmer_counts = HashMap::new();

    if sequence.len() < k {
        return kmer_counts;
    }

    for i in 0..=sequence.len() - k {
        let kmer = sequence[i..i + k].to_vec();
        // Only count k-mers with valid bases (ACGTN)
        if kmer.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n')) {
            *kmer_counts.entry(kmer).or_insert(0) += 1;
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
        let entropy = shannon_entropy(&kmer_counts, total_kmers);

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

        let entropy = shannon_entropy(&counts, 4);
        assert!((entropy - 2.0).abs() < 0.001); // Perfect entropy for 4 equal kmers
    }

    #[test]
    fn test_shannon_entropy_low_complexity() {
        let mut counts = HashMap::new();
        counts.insert(vec![b'A', b'A'], 10);

        let entropy = shannon_entropy(&counts, 10);
        assert_eq!(entropy, 0.0); // All same kmer = no entropy
    }

    #[test]
    fn test_get_kmers() {
        let sequence = b"ACGTACGT";
        let kmers = get_kmers(sequence, 3);

        assert_eq!(kmers.get(&vec![b'A', b'C', b'G']).unwrap(), &2);
        assert_eq!(kmers.get(&vec![b'C', b'G', b'T']).unwrap(), &2);
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
        // Use a longer sequence to avoid edge effects
        let sequence = b"ACGTACGTACGTACGTACGTACGT";
        let quality = b"IIIIIIIIIIIIIIIIIIIIIIII";

        let (masked_seq, masked_qual) = mask_sequence(sequence, quality, 10, 0.5, 3);

        // High complexity sequence should not be masked
        assert_eq!(masked_seq, sequence);
        assert_eq!(masked_qual, quality);
    }
}
