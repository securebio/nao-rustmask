#!/bin/bash
# Test compression impact on performance

set -euo pipefail

INPUT="synthetic_ont.fastq"
MASK_FASTQ="../modules/local/maskRead/mask_fastq/target/release/mask_fastq"

if [[ ! -f "$INPUT" ]]; then
    echo "Error: $INPUT not found"
    exit 1
fi

echo "Testing compression impact on mask_fastq performance"
echo "===================================================="
echo ""

echo "Test 1: With gzip compression (current default)"
echo "------------------------------------------------"
time (cat "$INPUT" | "$MASK_FASTQ" -w 25 -e 0.55 -k 5 > /tmp/test_compressed.fastq.gz) 2>&1 | grep "real"

echo ""
echo "Test 2: Without compression (piped to /dev/null)"
echo "------------------------------------------------"
time (cat "$INPUT" | "$MASK_FASTQ" -w 25 -e 0.55 -k 5 | zcat > /dev/null) 2>&1 | grep "real"

echo ""
echo "Test 3: Uncompressed output (if we add the option)"
echo "------------------------------------------------"
echo "(Would need to modify mask_fastq to support --no-compress flag)"

echo ""
echo "File sizes:"
ls -lh /tmp/test_compressed.fastq.gz 2>/dev/null || echo "N/A"
