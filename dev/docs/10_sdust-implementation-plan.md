# SDUST Algorithm Implementation Plan for Rustmasker

## Project Context

**Repository**: `securebio/nao-rustmask`
**Branch**: `main`
**Current Status**: CLI flags have been renamed in preparation for multi-algorithm support

### What Has Been Completed

1. **Flag Renaming**
   - Renamed `-e/--entropy` → `-t/--threshold` (more generic for both algorithms)
   - Renamed `-t/--threads` → `-j/--threads` (follows GNU Make convention)
   - Updated all documentation in README.md
   - Updated benchmark script: `scripts/benchmark_vs_bbmask.sh`
   - Commit: `c6e0898 - Rename CLI flags in preparation for multi-algorithm support`

2. **Current Binary Capabilities**
   - Entropy-based masking (BBMask-compatible)
   - Two implementation methods: HashMap and Array-based (optimized for k≤7)
   - Multi-threaded processing with Rayon
   - Streaming architecture for memory efficiency
   - Parallel compression support

### Current CLI Interface

```bash
-i, --input              # Input FASTQ file
-o, --output             # Output FASTQ file
-w, --window             # Window size (default: 80)
-t, --threshold          # Threshold (default: 0.70)
-k, --kmer               # K-mer size (default: 5)
-m, --method             # Entropy method: auto/array/hashmap (default: auto)
-c, --compression-level  # Gzip level (default: auto)
-j, --threads            # Number of threads (default: auto)
-s, --chunk-size         # Reads per chunk (default: 1000)
```

## Objective

Add support for the symmetric DUST algorithm (sdust) as an alternative masking method, allowing users to choose between entropy-based masking (current) and sdust masking (new). Output when running in sdust mode should be identical to that from Heng Li's sdust program (https://github.com/lh3/sdust).

## Target CLI Interface

```bash
# Common flags
-i, --input              # Input FASTQ file
-o, --output             # Output FASTQ file

# Algorithm selection (NEW)
-a, --algorithm          # Choose: entropy (default) or sdust

# Masking parameters (algorithm-dependent defaults)
-w, --window             # Window size (entropy: 80, sdust: 64)
-t, --threshold          # Masking threshold (entropy: 0.70, sdust: 20)
-k, --kmer               # K-mer size for entropy only (default: 5)
-m, --method             # Entropy computation method only (default: auto)

# Performance/output options
-j, --threads            # Number of threads (default: auto)
-s, --chunk-size         # Reads per chunk (default: 1000)
-c, --compression-level  # Gzip level (default: auto)
```

### Usage Examples

```bash
# Entropy masking (default, current behavior)
rustmasker -i input.fastq.gz -o output.fastq.gz

# Entropy with custom parameters
rustmasker -i input.fastq.gz -o output.fastq.gz \
  --algorithm entropy -w 80 -t 0.70 -k 5 -m array

# SDUST masking
rustmasker -i input.fastq.gz -o output.fastq.gz \
  --algorithm sdust -w 64 -t 20

# Short form
rustmasker -i input.fastq.gz -o output.fastq.gz -a sdust -w 64 -t 20
```

## SDUST Algorithm Background

### What is SDUST?

sdust is a symmetric DUST algorithm implementation by Heng Li that identifies low-complexity regions in DNA sequences. It gives output nearly identical to NCBI's dustmasker but runs ~4x faster.

**Reference**: https://github.com/lh3/sdust

### Key Characteristics

1. **Fixed triplet-based**: Always uses 3-mers (triplets), not configurable like entropy k-mers
2. **Window-based scoring**: Slides a window and scores based on triplet repetition
3. **Integer threshold**: Uses integer scoring (default: 20), not float like entropy
4. **Perfect intervals**: Identifies contiguous regions that consistently exceed threshold

### Algorithm Parameters

- **W (window size)**: Default 64 bases
- **T (threshold)**: Default 20 (integer)
- **Word length**: Fixed at 3 (triplets)

### Core Algorithm Logic

From the C implementation analysis:

```c
// Triplet encoding: A=0, C=1, G=2, T=3 (2 bits per base)
// Window scoring: if (repetition_count * 10 > T * window_length) → mask

// Example scoring logic:
int should_mask = (r * 10 > T * l);  // r=repetitions, l=length, T=threshold
```

**Key steps**:
1. Convert bases to triplet words (3 consecutive bases → 6-bit code)
2. Slide window over triplet sequence
3. Score each window based on triplet repetition frequency
4. Identify "perfect intervals" where score consistently exceeds threshold
5. Merge overlapping intervals
6. Mask identified regions

### Differences from Entropy Method

| Aspect | Entropy | SDUST |
|--------|---------|-------|
| K-mer size | Configurable (1-15) | Fixed at 3 (triplets) |
| Threshold | Float (0.0-1.0) | Integer (typically 20) |
| Metric | Shannon entropy | Repetition scoring |
| Complexity | O(1) with array, O(k) with HashMap | O(window_size) per position |
| Output | BBMask-compatible | dustmasker-compatible |

## Implementation Plan

### Phase 1: Core Data Structures (lib.rs)

**Location**: `/home/user/nao-rustmask/rustmasker/src/lib.rs`

#### 1.1 Triplet Encoding

```rust
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
```

#### 1.2 SDUST Scoring Structure

```rust
/// Represents a region to be masked
#[derive(Debug, Clone, Copy)]
struct MaskRegion {
    start: usize,  // Start position (inclusive)
    end: usize,    // End position (exclusive)
}

/// SDUST window scorer
struct SdustScorer {
    window_size: usize,
    threshold: i32,
}

impl SdustScorer {
    pub fn new(window_size: usize, threshold: i32) -> Self {
        Self { window_size, threshold }
    }

    /// Calculate score for a window of triplets
    /// Returns (repetition_count, window_length)
    fn score_window(&self, triplets: &[u8]) -> (usize, usize) {
        let mut counts = [0u16; 64];  // 64 possible triplet values

        for &triplet in triplets {
            counts[triplet as usize] += 1;
        }

        // Find the most repeated triplet
        let max_count = counts.iter().max().copied().unwrap_or(0);
        (max_count as usize, triplets.len())
    }

    /// Check if a window should be masked
    fn should_mask(&self, r: usize, l: usize) -> bool {
        r * 10 > self.threshold as usize * l
    }
}
```

#### 1.3 Main SDUST Masking Function

```rust
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

/// Find regions to mask using SDUST algorithm
fn find_dust_regions(
    triplets: &[u8],
    positions: &[usize],
    window_size: usize,
    threshold: i32,
) -> Vec<MaskRegion> {
    let scorer = SdustScorer::new(window_size, threshold);
    let mut regions = Vec::new();

    if triplets.len() < window_size {
        // Score entire sequence
        let (r, l) = scorer.score_window(triplets);
        if scorer.should_mask(r, l) {
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
        let (r, l) = scorer.score_window(window);

        if scorer.should_mask(r, l) {
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
```

### Phase 2: CLI Integration (rustmasker.rs)

**Location**: `/home/user/nao-rustmask/rustmasker/src/bin/rustmasker.rs`

#### 2.1 Add Algorithm Enum

```rust
/// Algorithm for masking
#[derive(ValueEnum, Clone, Debug)]
enum Algorithm {
    /// Shannon entropy-based masking (BBMask-compatible)
    Entropy,
    /// Symmetric DUST algorithm (sdust/dustmasker-compatible)
    Sdust,
}
```

#### 2.2 Update Args Structure

```rust
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    // ... existing input/output args ...

    /// Masking algorithm to use
    #[arg(short = 'a', long, value_enum, default_value = "entropy")]
    algorithm: Algorithm,

    /// Window size (default: entropy=80, sdust=64)
    #[arg(short = 'w', long)]
    window: Option<usize>,

    /// Masking threshold (default: entropy=0.70, sdust=20)
    #[arg(short = 't', long)]
    threshold: Option<String>,  // String to handle both float and int

    // ... rest of args ...
}
```

#### 2.3 Add Validation Logic

```rust
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Set algorithm-specific defaults
    let window = args.window.unwrap_or(match args.algorithm {
        Algorithm::Entropy => 80,
        Algorithm::Sdust => 64,
    });

    // Parse threshold based on algorithm
    let (entropy_threshold, sdust_threshold) = match args.algorithm {
        Algorithm::Entropy => {
            let t = args.threshold
                .unwrap_or_else(|| "0.70".to_string())
                .parse::<f64>()?;
            if t < 0.0 || t > 1.0 {
                eprintln!("Error: Entropy threshold must be between 0.0 and 1.0");
                std::process::exit(1);
            }
            (t, 0)
        }
        Algorithm::Sdust => {
            let t = args.threshold
                .unwrap_or_else(|| "20".to_string())
                .parse::<i32>()?;
            if t < 0 {
                eprintln!("Error: SDUST threshold must be positive");
                std::process::exit(1);
            }
            (0.0, t)
        }
    };

    // Warn if algorithm-specific flags are used with wrong algorithm
    if matches!(args.algorithm, Algorithm::Sdust) {
        if args.kmer != 5 {  // 5 is default
            eprintln!("Warning: -k/--kmer is ignored with sdust algorithm (always uses triplets)");
        }
        if !matches!(args.method, Method::Auto) {
            eprintln!("Warning: -m/--method is ignored with sdust algorithm");
        }
    }

    // ... rest of main ...
}
```

#### 2.4 Update Processing Logic

```rust
fn process_and_write_chunk(
    chunk: &mut Vec<FastqRecord>,
    writer: &mut Box<dyn Write>,
    args: &Args,
    // Pass parsed parameters
    window: usize,
    entropy_threshold: f64,
    sdust_threshold: i32,
) -> Result<(), Box<dyn std::error::Error>> {
    let results: Vec<(Vec<u8>, Vec<u8>)> = chunk
        .par_iter()
        .map(|record| {
            match args.algorithm {
                Algorithm::Entropy => {
                    match args.method {
                        Method::Auto => mask_sequence_auto(
                            &record.seq,
                            &record.qual,
                            window,
                            entropy_threshold,
                            args.kmer,
                        ),
                        Method::Array => mask_sequence_array(
                            &record.seq,
                            &record.qual,
                            window,
                            entropy_threshold,
                            args.kmer,
                        ),
                        Method::Hashmap => mask_sequence(
                            &record.seq,
                            &record.qual,
                            window,
                            entropy_threshold,
                            args.kmer,
                        ),
                    }
                }
                Algorithm::Sdust => {
                    mask_sequence_sdust(
                        &record.seq,
                        &record.qual,
                        window,
                        sdust_threshold,
                    )
                }
            }
        })
        .collect();

    // ... write results ...
    Ok(())
}
```

### Phase 3: Testing

#### 3.1 Unit Tests (lib.rs)

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_triplet_encoding() {
        assert_eq!(encode_triplet(b"AAA"), Some(0b000000));
        assert_eq!(encode_triplet(b"AAC"), Some(0b000001));
        assert_eq!(encode_triplet(b"AAG"), Some(0b000010));
        assert_eq!(encode_triplet(b"AAT"), Some(0b000011));
        assert_eq!(encode_triplet(b"CCC"), Some(0b010101));
        assert_eq!(encode_triplet(b"GGG"), Some(0b101010));
        assert_eq!(encode_triplet(b"TTT"), Some(0b111111));
        assert_eq!(encode_triplet(b"AAN"), None);
        assert_eq!(encode_triplet(b"AA"), None);
    }

    #[test]
    fn test_sdust_homopolymer() {
        let seq = b"AAAAAAAAAAAAAAAA";  // 16 A's
        let qual = vec![b'I'; 16];
        let (masked, _) = mask_sequence_sdust(seq, &qual, 10, 20);

        // Homopolymer should be masked
        assert!(masked.iter().all(|&b| b == b'N'));
    }

    #[test]
    fn test_sdust_high_complexity() {
        let seq = b"ACGTACGTACGTACGT";  // High complexity
        let qual = vec![b'I'; 16];
        let (masked, _) = mask_sequence_sdust(seq, &qual, 10, 20);

        // Should not be masked
        assert_eq!(masked, seq);
    }

    #[test]
    fn test_sdust_with_n() {
        let seq = b"AAAANNNNGGGG";
        let qual = vec![b'I'; 12];
        let (masked, _) = mask_sequence_sdust(seq, &qual, 4, 20);

        // A's should be masked, N's stay as is, G's should be masked
        // (implementation dependent on exact algorithm)
    }
}
```

#### 3.2 Integration Testing

Create test data and compare with sdust reference implementation:

```bash
# Generate test sequences
cd scripts
./generate_test_data.py -n 1000 -l 150 -o test.fastq

# Run sdust reference
sdust test.fasta > sdust_output.bed

# Run rustmasker with sdust
rustmasker -i test.fastq -o rustmasker_sdust.fastq -a sdust -w 64 -t 20

# Compare masked regions
# (Write comparison script)
```

### Phase 4: Documentation Updates

#### 4.1 README.md Updates

Add sections:
- Update feature list to mention "entropy or sdust masking"
- Add "Choosing an Algorithm" section
- Document sdust parameters
- Add comparison table: entropy vs sdust
- Update usage examples

#### 4.2 Help Text Updates

Update the terminal help message in rustmasker.rs to explain the two algorithms.

### Phase 5: Benchmark Script Updates

Update `scripts/benchmark_vs_bbmask.sh` to optionally support sdust benchmarking:

```bash
# Add sdust comparison option
-a, --algorithm    Algorithm: entropy or sdust (default: entropy)

# When --algorithm sdust:
# - Compare against sdust binary instead of BBMask
# - Use appropriate parameters
```

## Implementation Notes

### Multi-threading Support

SDUST masking will automatically benefit from the existing Rayon-based parallelization:
- Each read is processed independently
- No shared state between reads
- Read-level parallelism is sufficient (don't need window-level parallelism)

### Performance Considerations

1. **Triplet encoding is fast**: Simple bitwise operations
2. **Window scoring**: O(window_size) per position, but window is small (64)
3. **Memory overhead**: Minimal - just stores triplet array temporarily
4. **Expected performance**: Should be comparable to entropy method

### Edge Cases to Handle

1. **Sequences shorter than window**: Score entire sequence
2. **N bases**: Break into segments, process each segment independently
3. **Very low threshold**: May mask everything
4. **Very high threshold**: May mask nothing

### Compatibility Goals

- **Input/Output**: Same FASTQ format as entropy method
- **Compression**: Same gzip support
- **Threading**: Same multi-core support
- **Quality masking**: Replace with '#' just like entropy

## Testing Strategy

### Correctness Testing

1. **Unit tests**: Test triplet encoding, scoring, region finding
2. **Known sequences**: Test on sequences with known complexity
3. **Reference comparison**: Compare with sdust binary output
4. **Round-trip**: Ensure masked regions are consistent

### Performance Testing

1. **Benchmark against sdust**: Compare speed on same input
2. **Thread scaling**: Verify multi-threading works correctly
3. **Memory profiling**: Ensure no memory leaks or excessive usage

### Regression Testing

1. **Entropy still works**: Ensure existing functionality unchanged
2. **All tests pass**: Run existing test suite
3. **Documentation accuracy**: Verify all examples work

## Files to Modify

1. **rustmasker/src/lib.rs** - Add SDUST implementation (~300-400 lines)
2. **rustmasker/src/bin/rustmasker.rs** - Add CLI integration (~100 lines of changes)
3. **README.md** - Add documentation (~100 lines)
4. **Cargo.toml** - No changes needed (no new dependencies)

## Files to Create

1. **rustmasker/src/sdust.rs** (optional) - Could separate SDUST code into own module
2. **scripts/benchmark_vs_sdust.sh** (optional) - Benchmark script for SDUST

## Estimated Complexity

- **Core algorithm**: 4-6 hours
- **CLI integration**: 2-3 hours
- **Testing**: 3-4 hours
- **Documentation**: 1-2 hours
- **Total**: 10-15 hours

## Success Criteria

1. ✅ SDUST algorithm produces reasonable masking output
2. ✅ CLI supports `-a/--algorithm` flag
3. ✅ Entropy masking still works (no regression)
4. ✅ Multi-threading works for both algorithms
5. ✅ Documentation is complete and accurate
6. ✅ Tests pass
7. ✅ Code compiles without warnings

## References

- **SDUST source**: https://github.com/lh3/sdust
- **NCBI dustmasker**: https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/dustmasker/
- **Paper describing symmetric DUST algorithm**: https://pubmed.ncbi.nlm.nih.gov/16796549/ (paywalled)

## Next Steps for Implementation

1. Start with triplet encoding function (simple, testable)
2. Implement window scoring logic
3. Implement region finding algorithm
4. Write unit tests for core functions
5. Integrate with CLI
6. Test with real data
7. Compare with sdust reference
8. Update documentation
9. Commit and push

## Open Questions / Decisions Needed

1. **Exact scoring formula**: Need to verify the exact SDUST scoring matches reference
2. **Perfect interval definition**: Need to understand how sdust defines "perfect intervals"
3. **Threshold mapping**: Confirm default threshold of 20 is appropriate
4. **Output format**: Should we log masking statistics like BBMask does?
5. **Benchmark targets**: What speedup should we target vs reference sdust?

## Known Issues / Risks

1. **Algorithm fidelity**: Getting exact match with sdust may be tricky
2. **Test data**: Need good test sequences with known complexity levels (consider if scripts/generate_test_data.py is sufficient for test data generation)
