# Benchmarking mask_fastq vs bbmask.sh

This directory contains scripts to benchmark the performance and memory usage of the new `mask_fastq` Rust utility compared to the original `bbmask.sh` from BBTools.

## Quick Start

### Run all benchmarks automatically:

```bash
cd benchmark
./run_all_benchmarks.sh
```

This will:
1. Generate 4 test datasets (small Illumina, medium Illumina, long ONT, ultra-long ONT)
2. Run benchmarks on each dataset
3. Compare speed and memory usage
4. Verify that outputs are identical

### Run a single benchmark:

```bash
# Generate test data
./generate_test_data.py --num-reads 1000 --read-length 10000 -o test.fastq

# Run benchmark
./run_benchmark.sh test.fastq
```

## Scripts

### 1. `generate_test_data.py`

Generates synthetic FASTQ data with configurable complexity.

**Usage:**
```bash
./generate_test_data.py [options] -o output.fastq

Options:
  -n, --num-reads N          Number of reads (default: 10000)
  -l, --read-length N        Average read length (default: 1000)
  -o, --output FILE          Output FASTQ file (required)
  --low-complexity FRAC      Fraction of low-complexity reads (default: 0.3)
  --medium-complexity FRAC   Fraction of medium-complexity reads (default: 0.3)
```

**Examples:**
```bash
# Illumina-like data (150bp)
./generate_test_data.py -n 10000 -l 150 -o illumina.fastq

# ONT-like data (50Kbp average)
./generate_test_data.py -n 1000 -l 50000 -o ont.fastq

# Ultra-long ONT reads (100Kbp+) - this is where memory issues occur!
./generate_test_data.py -n 100 -l 100000 -o ont_ultralong.fastq

# Custom complexity distribution (50% low, 20% medium, 30% high)
./generate_test_data.py -n 5000 -l 1000 --low-complexity 0.5 --medium-complexity 0.2 -o custom.fastq
```

### 2. `run_benchmark.sh`

Runs both tools on the same input and compares performance.

**Usage:**
```bash
./run_benchmark.sh <input.fastq> [options]

Options:
  -w, --window N    Window size (default: 25)
  -e, --entropy N   Entropy threshold (default: 0.55)
  -k, --kmer N      K-mer size (default: 5)
```

**Example:**
```bash
./run_benchmark.sh test.fastq
./run_benchmark.sh test.fastq --window 50 --entropy 0.6
```

### 3. `run_all_benchmarks.sh`

Runs comprehensive benchmarks across multiple scenarios.

**Usage:**
```bash
./run_all_benchmarks.sh
```

## Understanding the Results

### Output Example:

```
╔══════════════════════════════════════════════════════════════╗
║                      BENCHMARK RESULTS                       ║
╚══════════════════════════════════════════════════════════════╝

Metric                    bbmask.sh            mask_fastq
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Runtime                   2.456s               0.234s
Peak Memory (MB)          11895.00             4.52
Masked bases              150/480              150/480
Masked percentage         -                    31.250%
Outputs match             YES

Speedup: 10.50x faster
Memory reduction: 99.96%
```

### Key Metrics:

- **Runtime**: Total execution time (wall clock)
- **Peak Memory**: Maximum resident set size (RSS) in MB
- **Masked bases**: Number of bases converted to 'N'
- **Outputs match**: Whether both tools produce identical results
- **Speedup**: How many times faster mask_fastq is
- **Memory reduction**: Percentage reduction in memory usage

## Expected Results

### For short reads (Illumina-like, ~150bp):
- **Speed**: mask_fastq should be 2-5x faster
- **Memory**: mask_fastq should use 10-100x less memory
- **Output**: Should be identical

### For long reads (ONT, 10-50Kbp):
- **Speed**: mask_fastq should be 5-15x faster
- **Memory**: mask_fastq should use 100-1000x less memory
- **Output**: Should be identical

### For ultra-long reads (ONT, 100Kbp+):
- **Speed**: mask_fastq should be 10-20x faster
- **Memory**: **This is the critical test!**
  - bbmask.sh may run out of memory (OOM) or use many GB
  - mask_fastq should use only a few MB (streaming)
- **Output**: Should be identical

## Why mask_fastq is Faster and Uses Less Memory

### bbmask.sh (Java/BBTools):
- Loads entire FASTQ file into Java heap memory
- For a 100Kbp read, needs to allocate large arrays in memory
- Multiple reads = multiplicative memory usage
- JVM overhead adds additional memory

### mask_fastq (Rust):
- Streams FASTQ records one at a time
- Processes each record independently
- Memory usage is O(read_length), not O(file_size)
- No JVM overhead
- Efficient binary (793KB vs 300MB BBTools suite)

## Troubleshooting

### "mask_fastq binary not found"
```bash
cd ../modules/local/maskRead/mask_fastq
cargo build --release
```

### "bbmask.sh not found"
Install BBTools:
```bash
# On Linux/macOS with Homebrew
brew install bbtools

# Or download from: https://sourceforge.net/projects/bbmap/
```

### "/usr/bin/time not found"
Memory statistics will be limited. Install GNU time:
```bash
# Ubuntu/Debian
sudo apt-get install time

# macOS
brew install gnu-time
```

### "Outputs differ"
This shouldn't happen! If it does, investigate:
```bash
# Compare the outputs manually
cd results
diff bbmask_seqs.txt mask_fastq_seqs.txt | head -20
```

## Files Generated

```
benchmark/
├── test_data/              # Generated FASTQ files
│   ├── small_illumina.fastq
│   ├── medium_illumina.fastq
│   ├── long_ont.fastq
│   └── ultralong_ont.fastq
├── results/                # Benchmark outputs
│   ├── bbmask_output.fastq.gz
│   ├── mask_fastq_output.fastq.gz
│   ├── bbmask_log.txt
│   ├── bbmask_time.txt
│   └── mask_fastq_time.txt
└── *.sh, *.py             # Benchmark scripts
```

## Integration with Issue #323

These benchmarks directly test the fix for issue #323, which reported:
> "MASK_FASTQ_READS process memory usage is too high with large ONT reads"

The ultra-long read benchmark (100Kbp+) specifically targets this issue and demonstrates that mask_fastq solves the memory problem while maintaining identical output.
