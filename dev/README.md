# Development Archive

This directory contains historical development artifacts from the creation of `mask_fastq`. These materials document the development process, optimization decisions, and experimental implementations. **They are preserved for reference but are not actively maintained.**

## Contents

### `binaries/` - Development Binary Tools

Alternative implementations and profiling tools used during development:

- **`mask_fastq_array.rs`** - Forces array-based entropy calculation (vs. auto-selection). Useful for benchmarking the array optimization specifically.
- **`mask_fastq_profiled.rs`** - Instrumented version that tracks time spent in I/O, masking, and writing. Used to identify bottlenecks.
- **`microbench.rs`** - Component-level microbenchmarks for k-mer encoding, entropy tracker operations, and sliding window performance.
- **`memory_benchmark.rs`** - Memory usage benchmarking tool for analyzing memory scaling with different k values, chunk sizes, and thread counts.
- **`cache_test.rs`** - Cache performance testing tool for comparing array vs HashMap approaches at different k-mer sizes.
- **`encoding_benchmark.rs`** - Performance comparison of u16 vs u32 k-mer encoding to inform decision on supporting k≤15.

These binaries are **not maintained** and may not compile with the current library code without updates.

### `docs/` - Development Documentation

Technical analyses and decision logs created during development:

1. **`1_HANDOFF_CONTEXT.md`** - Project history, initial implementation notes, and handoff context
2. **`2_BITPACKING_ANALYSIS.md`** - Analysis of u16 bit-packing optimization for k-mer encoding
3. **`3_PARALLEL_IMPLEMENTATION.md`** - Design decisions for multi-threaded `mask_fastq_parallel`
4. **`4_OPTIMIZATION_RESULTS.md`** - Performance improvements from various optimizations
5. **`5_PROFILING_ANALYSIS.md`** - Bottleneck identification through profiling
6. **`6_SIMD_PORTABILITY.md`** - Investigation of SIMD optimizations and portability concerns
7. **`7_PARALLELIZATION_ANALYSIS.md`** - Parallel processing strategies and tradeoffs
8. **`MEMORY_ANALYSIS.md`** - Memory scaling analysis for extending k-mer support from k≤8 to k≤15, including benchmarking results and u32 migration recommendations

These documents capture important context about **why** certain design decisions were made, including:
- Why we use strand-specific (not canonical) k-mers
- Why entropy is normalized to [0,1]
- Why array-based tracking is limited to k≤7
- Why SIMD was not implemented
- Memory usage patterns and optimization strategies
- Why u32 encoding was chosen over maintaining both u16 and u32

### `benchmarks/` - Development Benchmark Scripts

14 shell scripts used during development for testing and optimization:

**Core benchmarking:**
- `run_benchmark.sh` - Single file comparison (most comprehensive)
- `run_all_benchmarks.sh` - Full benchmark suite (no memory limits)
- `run_safe_benchmarks.sh` - Safe version with memory limits
- `baseline_benchmark.sh` - Initial baseline measurements

**Specialized testing:**
- `bench_parallel.sh` - Parallel version benchmarking
- `profile_performance.sh` - Performance profiling
- `test_compression_impact.sh` - Compression performance analysis
- `test_compression_levels.sh` - Compression level comparison

**Debugging/Analysis:**
- `check_memory.sh` - Memory availability checker
- `debug_differences.sh` - Output difference debugging
- `extract_diff_sequences.sh` - Extract sequences that differ
- `analyze_diff_patterns.sh` - Pattern analysis in differences

Also includes `README.md` with extensive benchmarking documentation.

## Why Preserve This?

These materials are valuable for:

1. **Understanding design decisions** - Why certain approaches were chosen over alternatives
2. **Future optimization** - Context for potential future improvements
3. **Debugging** - Historical context if issues arise
4. **Education** - Learning resource for similar projects

## Actively Maintained Code

For production use, see:

- **`mask_fastq/src/`** - Production library and binaries
- **`scripts/`** - Maintained utilities (test data generation, simple benchmarking)
- **Main `README.md`** - User-facing documentation

## Status

**Last updated**: 2025-11-12 (updated with k≤15 extension artifacts)
**Status**: Archived (not actively maintained)
**Purpose**: Historical reference and context preservation

---

*These materials were created by Claude Code during the initial implementation and optimization of mask_fastq to address issue #323 in the NAO mgs-workflow pipeline, and later extended to support k≤15 for BBMask compatibility.*
