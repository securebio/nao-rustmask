#!/bin/bash
# Debug script to analyze differences between bbmask and mask_fastq outputs

set -euo pipefail

if [[ ! -f "results/bbmask_output.fastq.gz" ]] || [[ ! -f "results/mask_fastq_output.fastq.gz" ]]; then
    echo "Error: Run a benchmark first to generate output files"
    exit 1
fi

echo "Analyzing differences between bbmask.sh and mask_fastq outputs"
echo "================================================================"
echo ""

# Extract sequences with IDs
echo "Extracting sequences with read IDs..."
zcat results/bbmask_output.fastq.gz | paste - - - - | awk '{print $1 "\t" $2}' | sort > /tmp/bbmask_with_ids.txt
zcat results/mask_fastq_output.fastq.gz | paste - - - - | awk '{print $1 "\t" $2}' | sort > /tmp/mask_fastq_with_ids.txt

# Find differences
echo "Finding differences..."
diff /tmp/bbmask_with_ids.txt /tmp/mask_fastq_with_ids.txt > /tmp/diff_output.txt || true

if [[ ! -s /tmp/diff_output.txt ]]; then
    echo "✅ No differences found! Outputs are identical."
    exit 0
fi

# Count differences
NUM_DIFFS=$(grep "^<" /tmp/diff_output.txt | wc -l)
echo "❌ Found $NUM_DIFFS differing sequences"
echo ""

# Show first few differences
echo "First 5 differing sequences:"
echo "----------------------------"
echo ""

grep "^< @" /tmp/diff_output.txt | head -5 | while read line; do
    READ_ID=$(echo "$line" | awk '{print $2}')
    echo "Read: $READ_ID"

    # Get bbmask version
    BBMASK_SEQ=$(grep "^< $READ_ID" /tmp/diff_output.txt | awk '{print $3}')
    echo "  bbmask:    $BBMASK_SEQ"

    # Get mask_fastq version
    MASK_FASTQ_SEQ=$(grep "^> $READ_ID" /tmp/diff_output.txt | awk '{print $3}')
    echo "  mask_fastq: $MASK_FASTQ_SEQ"

    # Count N's
    BBMASK_NS=$(echo "$BBMASK_SEQ" | grep -o "N" | wc -l)
    MASK_FASTQ_NS=$(echo "$MASK_FASTQ_SEQ" | grep -o "N" | wc -l)

    echo "  N count: bbmask=$BBMASK_NS, mask_fastq=$MASK_FASTQ_NS"
    echo ""
done

echo "Analysis:"
echo "---------"

# Categorize differences
BBMASK_MORE_MASKED=0
MASK_FASTQ_MORE_MASKED=0

grep "^< @" /tmp/diff_output.txt | while read line; do
    READ_ID=$(echo "$line" | awk '{print $2}')
    BBMASK_SEQ=$(grep "^< $READ_ID" /tmp/diff_output.txt | awk '{print $3}')
    MASK_FASTQ_SEQ=$(grep "^> $READ_ID" /tmp/diff_output.txt | awk '{print $3}')

    BBMASK_NS=$(echo "$BBMASK_SEQ" | grep -o "N" | wc -l)
    MASK_FASTQ_NS=$(echo "$MASK_FASTQ_SEQ" | grep -o "N" | wc -l)

    if [[ $BBMASK_NS -gt $MASK_FASTQ_NS ]]; then
        echo "$READ_ID: bbmask masks MORE (${BBMASK_NS} vs ${MASK_FASTQ_NS})"
    else
        echo "$READ_ID: mask_fastq masks MORE (${MASK_FASTQ_NS} vs ${BBMASK_NS})"
    fi
done | tee /tmp/masking_analysis.txt

echo ""
BBMASK_MORE=$(grep "bbmask masks MORE" /tmp/masking_analysis.txt | wc -l || echo 0)
MASK_FASTQ_MORE=$(grep "mask_fastq masks MORE" /tmp/masking_analysis.txt | wc -l || echo 0)

echo "Summary:"
echo "  Sequences where bbmask masks MORE:      $BBMASK_MORE"
echo "  Sequences where mask_fastq masks MORE:  $MASK_FASTQ_MORE"
echo ""
echo "Full diff saved to: /tmp/diff_output.txt"
