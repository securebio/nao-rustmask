# nao-rustmask

Fast, memory-efficient low-complexity masking for FASTQ files, written in Rust.

## Background

This project was created to solve memory issues with [BBTools' BBMask program](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/) when processing large files of sequencing reads. The original tool loads entire FASTQ files into Java heap memory, which caused out-of-memory errors on large fastq files in the NAO's ONT pipeline.

### Key Features

- **Identical output** to BBMask's entropy masking
- **Streaming architecture** controls memory usage on large files
- **Parallel processing**: Multi-core support for fast processing
- **Compatible I/O**: Reads and writes plain or gzipped FASTQ files
- **Flexible method selection**: Choose between array-based (fast) or HashMap (memory-efficient) algorithms

## What is Low-Complexity Masking?

Low-complexity regions in DNA sequences (homopolymers, tandem repeats, etc.) can cause false-positive alignments and complicate downstream analysis. This tool identifies and masks such regions by:

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

The compiled binary will be in `mask_fastq/target/release/mask_fastq`

## Usage

### Basic Usage

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

# Force specific method (array is default for k≤7)
mask_fastq -i input.fastq.gz -o output.fastq.gz --kmer 9 --method array

# Use all available CPU cores (default)
mask_fastq -i large.fastq.gz -o masked.fastq.gz

# Specify thread count
mask_fastq -i large.fastq.gz -o masked.fastq.gz --threads 8

# Adjust chunk size for memory/performance tradeoff
mask_fastq -i large.fastq.gz -o masked.fastq.gz \
  --chunk-size 5000 \
  --threads 8
```

### Command-Line Options

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--input` | `-i` | stdin | Input FASTQ file (plain or gzipped) |
| `--output` | `-o` | stdout | Output FASTQ file |
| `--window` | `-w` | 80 | Window size for entropy calculation |
| `--entropy` | `-e` | 0.70 | Entropy threshold (mask if < threshold) |
| `--kmer` | `-k` | 5 | K-mer size (1-15) |
| `--method` | `-m` | auto | Method: `auto` (adaptive), `array` (fast), or `hashmap` (memory-efficient) |
| `--compression-level` | `-c` | auto | Gzip compression level (0-9) |
| `--threads` | `-t` | auto | Number of threads to use |
| `--chunk-size` | `-s` | 1000 | Reads per chunk (affects memory usage) |

### Compression Behavior

The tool automatically handles compression based on context:

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
- Used automatically for k≤7 (or can be forced with `--method array`)

**U32 K-mer Encoding**:
- 2 bits per base (A=00, C=01, G=10, T=11)
- Supports k ≤ 15 (30 bits) for full BBMask compatibility
- Fast bitwise operations
- Efficient HashMap/array indexing
- No performance penalty vs u16 (benchmarked 0-9% faster)

**Adaptive Method Selection** (`--method auto`, default):
- Uses array-based approach for k≤7 (optimal cache performance)
- Switches to HashMap for k>7 (memory-efficient for larger k)
- Can be overridden with `--method array` or `--method hashmap`

## Benchmarking

Compare `mask_fastq` performance against BBMask:

```bash
cd scripts

# Generate test data
./generate_test_data.py -n 10000 -l 150 -o illumina.fastq

# Benchmark against BBMask
./benchmark_vs_bbmask.sh illumina.fastq

# Custom parameters
./benchmark_vs_bbmask.sh illumina.fastq --window 50 --entropy 0.6
```

The benchmark script:
- Runs both tools with identical parameters
- Compares runtime and memory usage
- Verifies outputs are identical
- Requires BBMask installed (for comparison only)

## Development Archive

Historical development materials (analyses, alternative implementations, extensive benchmark scripts) are preserved in the `dev/` directory for reference but are **not actively maintained**. See [`dev/README.md`](dev/README.md) for details.

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

- Currently only works on Fastq files, not Fasta

## Acknowledgments

- Inspired by BBTools' `bbmask.sh` by Brian Bushnell
- Array-based entropy tracker design based on BBTools' `EntropyTracker.java`
- Created to address issue [#323](https://github.com/securebio/nao-mgs-workflow/issues/323) in the [NAO's mgs-workflow pipeline](https://github.com/securebio/nao-mgs-workflow/issues/323)
