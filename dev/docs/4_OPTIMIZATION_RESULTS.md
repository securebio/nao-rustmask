# Array-Based Entropy Tracker Optimization Results

## Summary

Implemented BBMask-inspired array-based entropy tracker that achieves **1.7-3.2x speedup** over the original HashMap implementation while maintaining identical results and streaming memory architecture.

## Key Optimizations Implemented

### 1. **Array-Based K-mer Tracking**
- Replaced `HashMap<u16, usize>` with fixed-size `Vec<u16>` arrays
- Direct array indexing: O(1) lookups vs HashMap hash computation
- Memory usage: ~4KB for k=5, ~64KB for k=7 (reasonable overhead)

### 2. **Precalculated Entropy Table**
- Pre-compute `p*log2(p)` values during initialization
- Eliminates expensive `log2()` calls from hot path
- Storage: ~400 bytes for typical window size of 80

### 3. **Count-of-Counts Histogram**
- Track how many k-mers have each count value
- Enables O(1) entropy updates instead of O(unique_kmers) iteration
- Key insight from BBMask's EntropyTracker design

### 4. **Incremental Entropy Updates**
```rust
// Old: O(unique_kmers) - iterate HashMap on every window slide
let entropy = shannon_entropy(&kmer_counts, total_kmers);

// New: O(1) - just return cached value!
let entropy = tracker.entropy();
```

## Performance Results

### Benchmark Setup
- **Hardware**: Standard CPU (details in environment)
- **Test Data**:
  - Small: 1K reads × 150bp (Illumina-like)
  - Medium: 10K reads × 1Kbp (mixed complexity)
  - Long: 1K reads × 10Kbp (ONT-like)
- **Parameters**: k=5, window=25, entropy=0.55 (defaults)

### Results

| Dataset | HashMap (baseline) | Array-Based | Speedup |
|---------|-------------------|-------------|---------|
| Small (1K × 150bp) | 0.041s | 0.024s | **1.71x** |
| Medium (10K × 1Kbp) | 2.026s | 0.707s | **2.87x** |
| Long ONT (1K × 10Kbp) | 2.074s | 0.648s | **3.20x** |

### Key Findings

1. **Speedup increases with read length**: Longer reads have more window slides, each benefiting from O(1) entropy calculation
2. **Identical output**: Verified with `diff` - array and HashMap implementations produce byte-for-byte identical results
3. **Memory stays bounded**: Array overhead is fixed per read (~4-64KB depending on k), not dependent on read length
4. **Best for typical k values**: k=5-7 provides good balance of memory usage and speedup

## Implementation Details

### ArrayEntropyTracker Structure
```rust
pub struct ArrayEntropyTracker {
    counts: Vec<u16>,           // K-mer counts (size 4^k)
    count_counts: Vec<u16>,     // Histogram of count frequencies
    entropy_table: Vec<f64>,    // Precalculated p*log2(p) values
    entropy_mult: f64,          // Normalization factor
    current_esum: f64,          // Running entropy sum
    unique: usize,              // Number of unique k-mers
}
```

### Memory Usage by K-mer Size
- k=5: ~4 KB (1,024 kmers × 2 bytes)
- k=6: ~16 KB (4,096 kmers × 2 bytes)
- k=7: ~64 KB (16,384 kmers × 2 bytes)
- k=8: ~256 KB (65,536 kmers × 2 bytes)

## Comparison with BBMask

### What We Adopted from BBMask
✅ Array-based k-mer counting
✅ Count-of-counts histogram
✅ Precalculated entropy tables
✅ Incremental O(1) entropy updates

### What We Kept from Original Design
✅ **Streaming architecture** (BBMask loads entire file into memory)
✅ **Bounded memory usage** (BBMask uses ~1 byte per base for whole file)
✅ **Gzip streaming support**
✅ **Chunk-based parallelism**

### Result
**Best of both worlds**: BBMask's computational efficiency with streaming memory architecture

## Files Modified

1. **modules/local/maskRead/mask_fastq/src/lib.rs**
   - Added `ArrayEntropyTracker` struct (lines 189-333)
   - Added `mask_sequence_array()` function (lines 335-425)
   - Added comprehensive tests (lines 507-625)

2. **modules/local/maskRead/mask_fastq/src/bin/mask_fastq_array.rs**
   - New binary using array-based implementation
   - Drop-in replacement for `mask_fastq`

3. **benchmark/test_data/**
   - Generated synthetic test datasets
   - Small, medium, and long read datasets with varying complexity

## Testing

All tests pass:
```
test result: ok. 12 passed; 0 failed
```

Critical test: `test_mask_sequence_array_matches_hashmap`
- Verifies array and HashMap implementations produce identical results
- Tests multiple sequence types: homopolymers, repeats, random sequences

## Recommendations

### For Production Use
1. **Default to array-based implementation** for k ≤ 7 (covers typical use cases)
2. **Keep HashMap as fallback** for k > 7 or when memory is extremely constrained
3. **Consider hybrid approach**: Auto-select based on k-mer size

### Future Optimizations
1. **SIMD for k-mer encoding**: Process multiple bases simultaneously
2. **Batch entropy calculations**: Only compute when needed for masking decisions
3. **Profile-guided optimization**: Different speed modes like BBMask (FAST/MEDIUM/SLOW)

## Conclusion

The array-based entropy tracker successfully demonstrates that BBMask's core optimization techniques can be adapted to a streaming architecture. The **3.2x speedup** on long reads with **zero change in output** validates the approach while maintaining the memory-bounded design critical for processing large sequencing files.

## References

- BBMask source: https://github.com/BioInfoTools/BBMap/blob/master/current/jgi/BBMask.java
- EntropyTracker: https://github.com/BioInfoTools/BBMap/blob/master/current/structures/EntropyTracker.java
- BBMask Guide: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmask-guide/
