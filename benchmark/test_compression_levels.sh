#!/bin/bash
set -euo pipefail

INPUT="synthetic_ont.fastq"
MASK_FASTQ="../modules/local/maskRead/mask_fastq/target/release/mask_fastq"

echo "Testing different compression levels on mask_fastq"
echo "=================================================="
echo "Input: $INPUT (1000 reads, 5M bases)"
echo ""

for level in 0 1 3 6 9; do
    echo "Compression level $level:"
    OUTPUT="/tmp/test_c${level}.fastq.gz"
    time (cat "$INPUT" | "$MASK_FASTQ" -w 25 -e 0.55 -k 5 -c "$level" > "$OUTPUT") 2>&1 | grep "real"
    SIZE=$(ls -lh "$OUTPUT" | awk '{print $5}')
    echo "  Output size: $SIZE"
    echo ""
done

echo "Summary:"
echo "--------"
ls -lh /tmp/test_c*.fastq.gz | awk '{print $5, $9}'
