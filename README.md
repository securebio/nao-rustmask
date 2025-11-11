# nao-rustmask

Fast, memory-efficient low-complexity masking for FASTQ files, written in Rust.

## Background

This project was created to solve memory issues with [BBTools' BBMask program](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/) when processing large files of sequencing reads. The original tool loads entire FASTQ files into Java heap memory, which caused out-of-memory errors on large fastq files in the NAO's ONT pipeline.

### Key Features

- **Identical output** to BBMask's entropy masking
- **Streaming architecture** controls memory usage on large files
- **Parallel processing**: Multi-core support with `mask_fastq_parallel`
- **Compatible I/O**: Reads and writes plain or gzipped FASTQ files

## What is Low-Complexity Masking?

Low-complexity regions in DNA sequences (homopolymers, tandem repeats, etc.) can cause false-positive alignments and complicate downstream analysis. These tools identify and mask such regions by:

1. Computing Shannon entropy for k-mers in sliding windows
2. Masking windows that fall below an entropy threshold
3. Replacing masked bases with 'N' and quality scores with '#'

Example:
```
Original:  ACGTACGT AAAAAAAAAA GCTAGCTA
Masked:    ACGTACGT NNNNNNNNNN GCTAGCTA
           ↑ high    ↑ low      ↑ high
             entropy   entropy    entropy
```

## Installation

### Prerequisites

- Rust toolchain (1.70+)
- Cargo (comes with Rust)

### Build from Source

```bash
cd mask_fastq
cargo build --release
```

The compiled binaries will be in `target/release/`:
- `mask_fastq` - Single-threaded version
- `mask_fastq_parallel` - Multi-threaded version

## Usage

### Basic Usage: `mask_fastq`

Single-threaded version for moderate-sized files:

```bash
# Read from file, write to file
mask_fastq -i input.fastq.gz -o output.fastq.gz

# Use stdin/stdout
cat input.fastq.gz | mask_fastq > output.fastq

# Uncompressed output
mask_fastq -i input.fastq.gz -o output.fastq

# Custom parameters
mask_fastq -i input.fastq.gz -o output.fastq.gz \
  --window 50 \
  --entropy 0.6 \
  --kmer 7
```

### Parallel Processing: `mask_fastq_parallel`

Multi-threaded version for large files:

```bash
# Use all available CPU cores
mask_fastq_parallel -i large.fastq.gz -o masked.fastq.gz

# Specify thread count
mask_fastq_parallel -i large.fastq.gz -o masked.fastq.gz --threads 8

# Adjust chunk size for memory/performance tradeoff
mask_fastq_parallel -i large.fastq.gz -o masked.fastq.gz \
  --chunk-size 5000 \
  --threads 8
```

### Command-Line Options

#### Common Options (both tools)

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--input` | `-i` | stdin | Input FASTQ file (plain or gzipped) |
| `--output` | `-o` | stdout | Output FASTQ file |
| `--window` | `-w` | 25 | Window size for entropy calculation |
| `--entropy` | `-e` | 0.55 | Entropy threshold (mask if < threshold) |
| `--kmer` | `-k` | 5 | K-mer size (1-8) |
| `--compression-level` | `-c` | auto | Gzip compression level (0-9) |

#### Parallel-Only Options

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--threads` | `-t` | auto | Number of threads to use |
| `--chunk-size` | `-s` | 1000 | Reads per chunk (affects memory usage) |

### Compression Behavior

The tools automatically handle compression based on context:

- **Input**: Auto-detects plain or gzipped FASTQ
- **Output to stdout**: Uncompressed by default (use `-c` to compress)
- **Output to .gz file**: Compressed at level 1 by default
- **Output to other file**: Uncompressed (use `-c` to compress)

Examples:
```bash
# Compressed output (fast, level 1)
mask_fastq -i in.fastq.gz -o out.fastq.gz

# Compressed output (maximum, level 9)
mask_fastq -i in.fastq.gz -o out.fastq.gz -c 9

# Uncompressed output
mask_fastq -i in.fastq.gz -o out.fastq

# Compressed stdout
cat in.fastq | mask_fastq -c 6 > out.fastq.gz
```

## Algorithm Details

### Entropy Calculation

The tools use Shannon entropy normalized to [0, 1]:

```
H = -Σ(p_i × log2(p_i)) / log2(n)
```

Where:
- `p_i` = frequency of k-mer i in the window
- `n` = total k-mers in the window

### Masking Strategy

1. Slide a window across each sequence
2. Calculate entropy for k-mers in the window
3. If entropy < threshold, mask all bases in the window
4. Advance window by 1 base and repeat

This creates contiguous masked regions, matching BBMask behavior.

### Optimization Techniques

**Array-Based Entropy Tracker (k ≤ 7)**:
- Pre-allocates count arrays (4^k elements)
- Maintains count-of-counts histogram
- Pre-calculates entropy table
- O(1) entropy updates (vs O(k) for HashMap)
- Incremental sliding window with add/remove operations

**U16 K-mer Encoding**:
- 2 bits per base (A=00, C=01, G=10, T=11)
- Fits k ≤ 8 in 16 bits
- Fast bitwise operations
- Efficient HashMap/array indexing

## Benchmarking

The `benchmark/` directory contains benchmarking tools created by Claude Code during development.

```bash
cd benchmark

# Check system memory
./check_memory.sh

# Run safe benchmarks (recommended)
./run_safe_benchmarks.sh

# Generate custom test data
./generate_test_data.py -n 10000 -l 150 -o illumina.fastq

# Benchmark a specific file
./run_benchmark.sh illumina.fastq
```

See [benchmark/README.md](benchmark/README.md) for detailed documentation.

## Technical Documentation

The `benchmark/` directory also contains detailed technical analyses, also created by Claude code during development:

- `1_HANDOFF_CONTEXT.md` - Project history and development notes
- `2_BITPACKING_ANALYSIS.md` - U16 k-mer encoding optimization
- `3_PARALLEL_IMPLEMENTATION.md` - Multi-threading design
- `4_OPTIMIZATION_RESULTS.md` - Performance improvements
- `5_PROFILING_ANALYSIS.md` - Bottleneck identification
- `6_SIMD_PORTABILITY.md` - SIMD considerations
- `7_PARALLELIZATION_ANALYSIS.md` - Parallel processing strategies

## Correctness

Output should be identical to BBMask's `bbmask.sh` with equivalent parameters:

```bash
# BBMask command
bbmask.sh in=input.fastq out=bbmask_out.fastq \
  entropy=0.55 window=25 k=5

# Equivalent mask_fastq command
mask_fastq -i input.fastq -o mask_fastq_out.fastq \
  -e 0.55 -w 25 -k 5

# Verify identical output
diff bbmask_out.fastq mask_fastq_out.fastq
# (no output = files are identical)
```

## Limitations

- **K-mer size**: Maximum k=8 (u16 encoding limitation)
- Currently only works on Fastq files, not Fasta

## Acknowledgments

- Inspired by BBTools' `bbmask.sh` by Brian Bushnell
- Array-based entropy tracker design based on BBTools' `EntropyTracker.java`
- Created to address issue [#323](https://github.com/securebio/nao-mgs-workflow/issues/323) in the [NAO's mgs-workflow pipeline](https://github.com/securebio/nao-mgs-workflow/issues/323)
