# Memory Usage Analysis: nao-rustmask Parallel Implementation

## Executive Summary

Memory usage in the parallel array-based k-mer masking program scales according to:

```
Total Memory = (Tracker Memory × Active Threads) + (Chunk Memory × 2)
```

Where:
- **Tracker Memory** = `4^k × 2 bytes` (exponential in k)
- **Chunk Memory** = `chunk_size × (2 × read_length + overhead)` (linear in chunk size and read length)

For typical usage (k=7, 150bp reads, chunk_size=1000, 4 threads): **~0.74 MB total**

For k=15 with u32 encoding: **8 GB for 4 threads** (prohibitive) → HashMap approach required

---

## 1. Memory Components

### 1.1 ArrayEntropyTracker (Per-Thread)

Each thread maintains its own entropy tracker with these arrays:

#### Primary Memory: K-mer Counts Array
- **Size**: `4^k` entries × 2 bytes (u16 counts)
- **Purpose**: Stores count for each possible k-mer
- **Scaling**: Exponential base-4 in k

```
k=5:  1,024 entries = 2 KB
k=6:  4,096 entries = 8 KB
k=7:  16,384 entries = 32 KB
k=8:  65,536 entries = 128 KB
k=9:  262,144 entries = 512 KB
k=10: 1,048,576 entries = 2 MB
k=15: 1,073,741,824 entries = 2 GB (!)
```

#### Auxiliary Memory
- **count_counts**: `(window - k + 1 + 2)` × 2 bytes ≈ 40-200 bytes
- **entropy_table**: `(window - k + 1 + 2)` × 8 bytes ≈ 160-800 bytes
- **Total auxiliary**: <1 KB (negligible compared to counts array)

### 1.2 Chunk Memory (Shared Across Threads)

The parallel implementation processes reads in chunks:

```rust
struct FastqRecord {
    id: Vec<u8>,    // ~10-50 bytes
    seq: Vec<u8>,   // read_length bytes
    qual: Vec<u8>,  // read_length bytes
}
```

- **Per read**: ~`2 × read_length + 20` bytes
- **Input chunk**: `chunk_size` × per-read size
- **Output chunk**: `chunk_size` × per-read size (masked sequences)
- **Total chunk memory**: `2 × chunk_size × (2 × read_length + 20)` bytes

For 150bp reads:
```
Per read: ~320 bytes
chunk_size=1000: 0.61 MB total
chunk_size=5000: 3.05 MB total
```

---

## 2. Benchmark Results

### 2.1 Impact of K-mer Size (Exponential Scaling)

**Configuration**: 10,000 reads, 150bp, chunk_size=1000, 4 threads

| k | Tracker Memory (4 threads) | Chunk Memory | Total | Throughput |
|---|---------------------------|--------------|-------|------------|
| 5 | 0.01 MB | 0.61 MB | 0.62 MB | 48.28 Mbp/s |
| 7 | 0.13 MB | 0.61 MB | 0.74 MB | 39.33 Mbp/s |
| 8 | 0.50 MB | 0.61 MB | 1.11 MB | 15.90 Mbp/s |

**Key findings:**
- Memory increases by **4× for each increase in k**
- Performance **decreases** with larger k (cache effects)
- k=8 shows significant slowdown (128 KB array exceeds L1 cache)

### 2.2 Impact of Thread Count (Linear Scaling)

**Configuration**: k=7, 10,000 reads, 150bp, chunk_size=1000

| Threads | Tracker Memory | Chunk Memory | Total | Throughput |
|---------|---------------|--------------|-------|------------|
| 4 | 0.13 MB | 0.61 MB | 0.74 MB | 39.33 Mbp/s |
| 8 | 0.25 MB | 0.61 MB | 0.86 MB | 68.39 Mbp/s |

**Key finding:** Tracker memory scales **linearly** with thread count (each thread needs its own tracker)

### 2.3 Impact of Chunk Size (Linear Scaling)

**Configuration**: k=7, 10,000 reads, 150bp, 4 threads

| Chunk Size | Tracker Memory | Chunk Memory | Total |
|------------|---------------|--------------|-------|
| 1,000 | 0.13 MB | 0.61 MB | 0.74 MB |
| 5,000 | 0.13 MB | 3.05 MB | 3.18 MB |

**Key finding:** Chunk memory dominates for large chunk sizes

### 2.4 Impact of Read Length (Linear Scaling)

**Configuration**: k=7, 10,000 reads, chunk_size=1000, 4 threads

| Read Length | Tracker Memory | Chunk Memory | Total |
|-------------|---------------|--------------|-------|
| 150bp | 0.13 MB | 0.61 MB | 0.74 MB |
| 300bp | 0.13 MB | 1.18 MB | 1.31 MB |

**Key finding:** Read length only affects chunk memory, not tracker memory

---

## 3. Memory Scaling to k=15 (BBMask Compatibility)

### 3.1 The Problem with Array-Based Approach

To support k=15, we need u32 encoding (30 bits for 15 bases × 2 bits/base).

**Array size for k=15**: `4^15 = 1,073,741,824` entries × 2 bytes = **2 GB per tracker**

| k | Array Size | 1 thread | 4 threads | 8 threads |
|---|-----------|----------|-----------|-----------|
| 8 | 65,536 | 128 KB | 512 KB | 1 MB |
| 9 | 262,144 | 512 KB | 2 MB | 4 MB |
| 10 | 1,048,576 | 2 MB | 8 MB | 16 MB |
| 11 | 4,194,304 | 8 MB | 32 MB | 64 MB |
| 12 | 16,777,216 | 32 MB | 128 MB | 256 MB |
| 13 | 67,108,864 | 128 MB | 512 MB | 1 GB |
| 14 | 268,435,456 | 512 MB | 2 GB | 4 GB |
| 15 | 1,073,741,824 | **2 GB** | **8 GB** | **16 GB** |

### 3.2 Why This is Impractical

1. **Massive Memory Waste**: For a 25bp window with k=15:
   - Window contains only **11 k-mers** (25 - 15 + 1)
   - Array allocates space for **1 billion entries**
   - **Utilization: 0.000001%** (11 out of 1 billion used)

2. **Cache Performance**:
   - L1 cache: typically 32-64 KB per core
   - L2 cache: typically 256-512 KB per core
   - L3 cache: typically 8-32 MB shared
   - A 2 GB array would never fit in cache → constant main memory access → severe performance degradation

3. **System Requirements**:
   - 4 threads: 8 GB RAM just for trackers
   - 8 threads: 16 GB RAM just for trackers
   - Plus chunk memory, OS, other processes
   - Would require 16-32 GB total system RAM

### 3.3 The Solution: HashMap-Based Approach

The current implementation already uses `HashMap` for k > 7:

```rust
pub fn mask_sequence_auto(...) -> (Vec<u8>, Vec<u8>) {
    if k <= 7 {
        mask_sequence_array(...)  // Array-based (4^k memory)
    } else {
        mask_sequence(...)        // HashMap-based (O(actual_kmers) memory)
    }
}
```

For k=15 with HashMap:
- **Memory per tracker**: Stores only actual k-mers seen in window
- **Typical window (25bp)**: ~11 k-mers × 8 bytes (u32 key + u32 count) = **88 bytes**
- **With overhead**: <1 KB per tracker
- **4 threads**: <4 KB (vs 8 GB with array!)
- **Performance**: Good with modern hash tables (slightly slower than array for small k, but acceptable)

---

## 4. Recommendations

### 4.1 Current Implementation (k ≤ 8)

The current threshold (array for k ≤ 7, HashMap for k > 7) is **well-chosen**:
- k=7: 32 KB per tracker (reasonable)
- k=8: 128 KB per tracker (marginal, performance already degraded)
- Array approach becomes impractical beyond k=8

### 4.2 Extending to k=15 (BBMask Compatibility)

To support k ≤ 15:

1. **Change encoding**: Modify `encode_kmer()` to use u32 instead of u16
   ```rust
   pub fn encode_kmer(bases: &[u8]) -> Option<u32> {
       if bases.len() > 15 {  // was 8
           return None;
       }
       let mut encoded: u32 = 0;  // was u16
       // ... rest same
   }
   ```

2. **Keep HashMap for k > 7**: Do NOT use array-based approach for k > 8
   - The auto-selection logic is already correct
   - HashMap is the only viable approach for large k

3. **Update validation**: Change maximum k from 8 to 15 in CLIs

4. **Update HashMap types**: Change `HashMap<u16, usize>` to `HashMap<u32, usize>` throughout

### 4.3 Memory-Constrained Systems

For systems with limited memory, consider:
- **Reduce chunk_size**: Default 1000 is reasonable, but can go as low as 100
- **Limit threads**: Memory scales linearly with threads
- **Avoid k=8**: Use k=7 for better cache performance

### 4.4 Optimal Settings

For typical use cases:

| Scenario | k | Window | Chunk Size | Threads | Memory |
|----------|---|--------|------------|---------|--------|
| Fast, low memory | 5 | 25 | 1000 | 4 | ~0.6 MB |
| Balanced | 7 | 25 | 1000 | 4 | ~0.7 MB |
| BBMask compatible | 15 | 31 | 1000 | 4 | ~0.6 MB (HashMap) |
| High throughput | 7 | 25 | 5000 | 8 | ~3.3 MB |

---

## 5. Conclusion

The parallel array-based implementation has predictable memory scaling:

1. **For k ≤ 7**: Array approach is excellent (low memory, good cache performance)
2. **For k = 8**: Array approach is marginal (higher memory, cache misses)
3. **For k ≥ 9**: Array approach is impractical (prohibitive memory usage)

To support k ≤ 15 for BBMask compatibility:
- Change encoding to u32 (simple)
- Use HashMap for k > 7 (already implemented)
- Memory usage remains reasonable (< 1 MB for typical cases)

The current implementation's auto-selection strategy is sound and should be maintained when extending to k=15.
