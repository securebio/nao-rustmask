# Bit-Packing K-mer Optimization Analysis

## Context

After implementing the incremental sliding window optimization, `mask_fastq` achieved:
- **4.5x speedup** over initial implementation (5.060s → 1.120s)
- **2.1x slower** than BBMask (1.120s vs 0.537s)
- **96% memory reduction** vs BBMask (3.76MB vs 341MB)

This document analyzes the remaining performance gap and options for further optimization.

## Performance Gap Analysis

### Microbenchmark Results

Testing HashMap operations with different k-mer representations:

| K-mer Key Type | Time (ms) | Memory/Entry | Speedup |
|----------------|-----------|--------------|---------|
| `Vec<u8>` (current) | 4948 | ~48 bytes | 1.0x |
| `u16` bit-packed | 2083 | ~24 bytes | **2.38x** |

**Why bit-packing is faster:**
- Integer keys are stack-allocated (no heap allocations)
- Integer hashing is much faster than byte vector hashing
- Better cache locality
- Smaller HashMap memory footprint

### Current Implementation Bottlenecks

Estimated runtime breakdown:
- **60-70%**: HashMap operations (insert/remove/lookup with `Vec<u8>` keys)
- **10-15%**: Shannon entropy calculation
- **10-15%**: Window range masking
- **10-15%**: I/O and compression

**Projected impact of u16 bit-packing:**
- If HashMap ops are 60% of runtime: **1.56x overall speedup** → 0.72s (1.34x slower than BBMask)
- If HashMap ops are 70% of runtime: **1.72x overall speedup** → 0.65s (1.21x slower than BBMask)

## Bit-Packing Capacity

**Encoding scheme:** 2 bits per base (A=00, C=01, G=10, T=11)

| Integer Type | Bits Available | Max k-mer Size | Expected Speedup |
|--------------|----------------|----------------|------------------|
| u16 | 16 | k ≤ 8 | 2.4x (measured) |
| u32 | 32 | k ≤ 16 | 2.0x (estimated) |
| u64 | 64 | k ≤ 32 | 1.8x (estimated) |

**Current parameter:** k=5 (BBMask default) → requires 10 bits → fits in u16

## Implementation Options

### Option 1: u16 Only (Maximum Performance) ✅ **CHOSEN**

**Implementation:**
```rust
// Validate k parameter
if k > 8 {
    eprintln!("Error: k={} exceeds maximum supported value (k ≤ 8)", k);
    std::process::exit(1);
}

// Use HashMap<u16, usize> throughout
fn encode_kmer(bases: &[u8]) -> u16 {
    let mut encoded: u16 = 0;
    for &base in bases {
        let bits = match base {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => return None,  // Skip k-mers with N or invalid bases
        };
        encoded = (encoded << 2) | bits;
    }
    Some(encoded)
}
```

**Pros:**
- Maximum speedup: ~2.4x for HashMap operations
- Simplest implementation (one code path)
- Smallest memory footprint per HashMap entry
- Covers all common use cases (k=3 to k=7)

**Cons:**
- Limits k to ≤8
- Not an issue: k>8 is rare for low-complexity masking

**Use cases covered:**
- k=3: Very sensitive (detects "ATATAT")
- k=5: Standard (BBMask default)
- k=7: Less sensitive (longer repeats only)

### Option 2: u32 for All (Good Compromise)

**Implementation:**
```rust
// Validate k parameter
if k > 16 {
    eprintln!("Error: k={} exceeds maximum supported value (k ≤ 16)", k);
    std::process::exit(1);
}

// Use HashMap<u32, usize> for all k values
```

**Pros:**
- Handles k up to 16
- Still ~2.0x faster than Vec<u8>
- Simple single implementation
- More headroom for unusual use cases

**Cons:**
- Slightly slower than u16 for small k values
- Uses more memory per entry than u16

### Option 3: Dynamic Type Selection (Maximum Flexibility)

**Implementation:**
```rust
match k {
    1..=8 => mask_with_u16(sequence, k),    // 2.4x speedup
    9..=16 => mask_with_u32(sequence, k),   // 2.0x speedup
    17..=32 => mask_with_u64(sequence, k),  // 1.8x speedup
    _ => mask_with_vec(sequence, k),        // fallback
}
```

**Pros:**
- Optimal performance for any k value
- Supports very large k values
- Maximum flexibility

**Cons:**
- Code duplication across 4 implementations
- Complex to maintain
- Adds ~300-400 lines of code
- Overkill for this use case

## Recommendation

**Choose Option 1 (u16 only)** because:

1. **Covers all practical use cases:** k=3 to k=7 is standard for low-complexity masking
2. **Maximum performance:** 2.4x speedup for HashMap operations
3. **Simplest code:** One clear implementation path
4. **Smallest memory footprint:** Important for memory-constrained environments
5. **BBMask compatibility:** Matches BBMask's default k=5 perfectly

**Expected result after optimization:**
- Current: 1.120s (2.1x slower than BBMask)
- With u16: ~0.65-0.72s (1.2-1.3x slower than BBMask)
- Memory: Still 3.76MB (96% less than BBMask)

## Remaining Performance Gap

After u16 optimization, the remaining 1.2-1.3x gap vs BBMask is likely due to:

1. **Memory model** (20-30% of gap)
   - BBMask: Loads entire file into memory (341MB)
   - mask_fastq: Streams data (3.76MB)
   - In-memory processing has better cache locality

2. **Java JIT optimization** (10-20% of gap)
   - BBMask's hot loops get optimized at runtime
   - Aggressive inlining and vectorization

3. **I/O and compression** (10-20% of gap)
   - Output compression overhead
   - Different buffering strategies

**This is an acceptable tradeoff:** 96% memory savings with only 1.2-1.3x runtime overhead is excellent for the target use case (memory-constrained environments with large ONT reads).

## Implementation Notes

### Handling N Bases

K-mers containing N (ambiguous) bases should be skipped (not counted in entropy calculation):

```rust
fn encode_kmer(bases: &[u8]) -> Option<u16> {
    let mut encoded: u16 = 0;
    for &base in bases {
        let bits = match base {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => return None,  // N or invalid base
        };
        encoded = (encoded << 2) | bits;
    }
    Some(encoded)
}
```

### Backwards Compatibility

To maintain compatibility with existing tests:
- Keep the same command-line interface
- Keep the same default parameters (k=5, window=25, entropy=0.55)
- Output format unchanged
- All existing tests should pass without modification

## Testing Plan

1. **Unit tests:** Verify encoding/decoding functions
2. **Integration tests:** All existing cargo tests must pass
3. **Correctness:** Output must match BBMask exactly on all test datasets
4. **Performance:** Benchmark on ultralong_ont.fastq (expect ~0.65-0.72s)
5. **Edge cases:** Test with k=3, k=5, k=7, k=8 (and reject k=9 with error)
