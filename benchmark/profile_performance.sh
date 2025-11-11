#!/bin/bash
# Profile mask_fastq implementations to identify bottlenecks
# Uses manual instrumentation and timing

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

MASK_FASTQ="../modules/local/maskRead/mask_fastq/target/release/mask_fastq"
MASK_FASTQ_ARRAY="../modules/local/maskRead/mask_fastq/target/release/mask_fastq_array"

echo "=========================================="
echo "Profiling mask_fastq Implementations"
echo "=========================================="
echo ""

# Check for profiling tools
echo "Checking for profiling tools..."
if command -v perf &> /dev/null; then
    echo "✓ perf available"
    USE_PERF=1
else
    echo "✗ perf not available"
    USE_PERF=0
fi

if command -v valgrind &> /dev/null; then
    echo "✓ valgrind available"
    USE_VALGRIND=1
else
    echo "✗ valgrind not available"
    USE_VALGRIND=0
fi

echo ""
echo "Using basic timing profiling..."
echo ""

# Function to profile with time breakdown
profile_with_time() {
    local name=$1
    local binary=$2
    local input=$3
    local output=$4

    echo "================================================"
    echo "Profiling: $name"
    echo "Input: $input"
    echo "================================================"

    # Get file stats
    local size=$(du -h "$input" | cut -f1)
    local reads=$(($(wc -l < "$input") / 4))
    echo "File size: $size"
    echo "Reads: $reads"
    echo ""

    # Run multiple times for stable measurements
    echo "Running 5 iterations for stable timing..."
    local total_time=0
    for i in {1..5}; do
        local start=$(date +%s.%N)
        $binary -i "$input" -o "$output" 2>/dev/null
        local end=$(date +%s.%N)
        local runtime=$(echo "$end - $start" | bc)
        echo "  Iteration $i: ${runtime}s"
        total_time=$(echo "$total_time + $runtime" | bc)
    done

    local avg_time=$(echo "scale=3; $total_time / 5" | bc)
    echo ""
    echo "Average time: ${avg_time}s"
    echo "Throughput: $(echo "scale=0; $reads / $avg_time" | bc) reads/sec"
    echo ""

    # Cleanup
    rm -f "$output"
}

# Profile both implementations on different datasets
echo "=========================================="
echo "HASHMAP IMPLEMENTATION"
echo "=========================================="
echo ""

profile_with_time "HashMap - Small dataset" \
    "$MASK_FASTQ" \
    "test_data/small_illumina.fastq" \
    "/tmp/profile_out.fastq"

profile_with_time "HashMap - Medium dataset" \
    "$MASK_FASTQ" \
    "test_data/medium_mixed.fastq" \
    "/tmp/profile_out.fastq"

profile_with_time "HashMap - Long ONT dataset" \
    "$MASK_FASTQ" \
    "test_data/long_ont.fastq" \
    "/tmp/profile_out.fastq"

echo ""
echo "=========================================="
echo "ARRAY-BASED IMPLEMENTATION"
echo "=========================================="
echo ""

profile_with_time "Array - Small dataset" \
    "$MASK_FASTQ_ARRAY" \
    "test_data/small_illumina.fastq" \
    "/tmp/profile_out.fastq"

profile_with_time "Array - Medium dataset" \
    "$MASK_FASTQ_ARRAY" \
    "test_data/medium_mixed.fastq" \
    "/tmp/profile_out.fastq"

profile_with_time "Array - Long ONT dataset" \
    "$MASK_FASTQ_ARRAY" \
    "test_data/long_ont.fastq" \
    "/tmp/profile_out.fastq"

echo ""
echo "=========================================="
echo "Profiling complete!"
echo "=========================================="
echo ""
echo "For more detailed profiling, consider:"
echo "  1. Install perf: sudo apt-get install linux-tools-generic"
echo "  2. Install cargo-flamegraph: cargo install flamegraph"
echo "  3. Run: cargo flamegraph --bin mask_fastq_array -- -i test.fastq -o /dev/null"
