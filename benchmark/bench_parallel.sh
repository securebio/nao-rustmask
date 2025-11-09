#!/bin/bash
set -euo pipefail

INPUT="synthetic_ont.fastq"
MASK_FASTQ="../modules/local/maskRead/mask_fastq/target/release/mask_fastq"
MASK_FASTQ_PARALLEL="../modules/local/maskRead/mask_fastq/target/release/mask_fastq_parallel"

if [[ ! -f "$INPUT" ]]; then
    echo "Error: $INPUT not found"
    exit 1
fi

echo "Benchmarking mask_fastq vs mask_fastq_parallel"
echo "Input: $INPUT (1000 reads, 5M bases)"
echo "=============================================="
echo ""

echo "Original (single-threaded):"
time (cat "$INPUT" | "$MASK_FASTQ" -w 25 -e 0.55 -k 5 -c 1 > /tmp/bench_original.fastq.gz) 2>&1 | grep "real"

echo ""
echo "Parallel - 1 thread:"
time (cat "$INPUT" | "$MASK_FASTQ_PARALLEL" -w 25 -e 0.55 -k 5 -c 1 -t 1 > /tmp/bench_parallel_1t.fastq.gz) 2>&1 | grep "real"

echo ""
echo "Parallel - 2 threads:"
time (cat "$INPUT" | "$MASK_FASTQ_PARALLEL" -w 25 -e 0.55 -k 5 -c 1 -t 2 > /tmp/bench_parallel_2t.fastq.gz) 2>&1 | grep "real"

echo ""
echo "Parallel - 4 threads:"
time (cat "$INPUT" | "$MASK_FASTQ_PARALLEL" -w 25 -e 0.55 -k 5 -c 1 -t 4 > /tmp/bench_parallel_4t.fastq.gz) 2>&1 | grep "real"

echo ""
echo "Parallel - 8 threads:"
time (cat "$INPUT" | "$MASK_FASTQ_PARALLEL" -w 25 -e 0.55 -k 5 -c 1 -t 8 > /tmp/bench_parallel_8t.fastq.gz) 2>&1 | grep "real"

echo ""
echo "Verifying outputs match..."
for file in /tmp/bench_parallel_*.fastq.gz; do
    diff <(zcat /tmp/bench_original.fastq.gz) <(zcat "$file") >/dev/null && echo "✓ $(basename $file)" || echo "✗ $(basename $file) DIFFERS!"
done
