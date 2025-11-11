#!/bin/bash
# Baseline performance benchmarks for mask_fastq
# Measures current HashMap-based implementation

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MASK_FASTQ="../modules/local/maskRead/mask_fastq/target/release/mask_fastq"
MASK_FASTQ_PARALLEL="../modules/local/maskRead/mask_fastq/target/release/mask_fastq_parallel"

cd "$SCRIPT_DIR"

# Check if binaries exist
if [ ! -f "$MASK_FASTQ" ]; then
    echo "Error: mask_fastq binary not found at $MASK_FASTQ"
    echo "Please build first: cd ../modules/local/maskRead/mask_fastq && cargo build --release"
    exit 1
fi

echo "================================================"
echo "Baseline Performance Benchmarks"
echo "HashMap-based implementation"
echo "================================================"
echo ""

# Function to run benchmark
run_benchmark() {
    local name=$1
    local input=$2
    local output=$3
    local binary=$4
    local extra_args=$5

    echo "Test: $name"
    echo "  Input: $input"

    # Get input file size
    local size=$(du -h "$input" | cut -f1)
    echo "  File size: $size"

    # Count reads
    local reads=$(($(wc -l < "$input") / 4))
    echo "  Reads: $reads"

    # Run with time measurement
    echo "  Running..."
    /usr/bin/time -f "  Time: %E (wall clock)\n  Memory: %M KB (max resident)" \
        $binary -i "$input" -o "$output" $extra_args 2>&1 | grep -E "(Time:|Memory:)"

    # Clean up output
    rm -f "$output"
    echo ""
}

# Create output directory
mkdir -p results

echo "=== Single-threaded (mask_fastq) ==="
echo ""
run_benchmark "Small Illumina (1K reads x 150bp)" \
    "test_data/small_illumina.fastq" \
    "results/small_out.fastq" \
    "$MASK_FASTQ"

run_benchmark "Medium mixed (10K reads x 1Kbp)" \
    "test_data/medium_mixed.fastq" \
    "results/medium_out.fastq" \
    "$MASK_FASTQ"

run_benchmark "Long ONT (1K reads x 10Kbp)" \
    "test_data/long_ont.fastq" \
    "results/long_out.fastq" \
    "$MASK_FASTQ"

echo "=== Multi-threaded (mask_fastq_parallel) ==="
echo ""

run_benchmark "Small Illumina (1K reads x 150bp, 4 threads)" \
    "test_data/small_illumina.fastq" \
    "results/small_out.fastq" \
    "$MASK_FASTQ_PARALLEL" \
    "-t 4"

run_benchmark "Medium mixed (10K reads x 1Kbp, 4 threads)" \
    "test_data/medium_mixed.fastq" \
    "results/medium_out.fastq" \
    "$MASK_FASTQ_PARALLEL" \
    "-t 4"

run_benchmark "Long ONT (1K reads x 10Kbp, 4 threads)" \
    "test_data/long_ont.fastq" \
    "results/long_out.fastq" \
    "$MASK_FASTQ_PARALLEL" \
    "-t 4"

echo "=== Chunk size comparison (parallel) ==="
echo ""

run_benchmark "Medium mixed (chunk=100)" \
    "test_data/medium_mixed.fastq" \
    "results/medium_out.fastq" \
    "$MASK_FASTQ_PARALLEL" \
    "-t 4 -s 100"

run_benchmark "Medium mixed (chunk=1000, default)" \
    "test_data/medium_mixed.fastq" \
    "results/medium_out.fastq" \
    "$MASK_FASTQ_PARALLEL" \
    "-t 4 -s 1000"

run_benchmark "Medium mixed (chunk=5000)" \
    "test_data/medium_mixed.fastq" \
    "results/medium_out.fastq" \
    "$MASK_FASTQ_PARALLEL" \
    "-t 4 -s 5000"

echo "================================================"
echo "Baseline benchmarks complete!"
echo "================================================"
