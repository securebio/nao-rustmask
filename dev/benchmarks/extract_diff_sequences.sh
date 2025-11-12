#!/bin/bash

# Extract sequences where BBMask and mask_fastq outputs differ
# Usage: ./extract_diff_sequences.sh <original_fastq> [max_sequences]
#
# Outputs: diff_sequences.fastq with sequences that have different masking

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 <original_fastq> [max_sequences]"
    echo ""
    echo "Example: $0 test_data/small_illumina.fastq 20"
    echo ""
    echo "This will create diff_sequences.fastq with up to 20 sequences"
    echo "where BBMask and mask_fastq produce different results."
    exit 1
fi

INPUT_FASTQ="$1"
MAX_SEQS="${2:-10}"  # Default to 10 sequences if not specified
RESULTS_DIR="results"

# Check if benchmark results exist
if [ ! -f "$RESULTS_DIR/bbmask_seqs.txt" ] || [ ! -f "$RESULTS_DIR/mask_fastq_seqs.txt" ]; then
    echo "Error: Benchmark results not found in $RESULTS_DIR/"
    echo "Please run the benchmark first: ./run_benchmark.sh $INPUT_FASTQ"
    exit 1
fi

echo "Analyzing differences between BBMask and mask_fastq outputs..."

# Find sequences that differ
diff --changed-group-format='%<' --unchanged-group-format='' \
    "$RESULTS_DIR/bbmask_seqs.txt" "$RESULTS_DIR/mask_fastq_seqs.txt" \
    > "$RESULTS_DIR/bbmask_only_seqs.txt" || true

diff --changed-group-format='%<' --unchanged-group-format='' \
    "$RESULTS_DIR/mask_fastq_seqs.txt" "$RESULTS_DIR/bbmask_seqs.txt" \
    > "$RESULTS_DIR/mask_fastq_only_seqs.txt" || true

# Count differences
BBMASK_DIFF=$(wc -l < "$RESULTS_DIR/bbmask_only_seqs.txt")
MASK_FASTQ_DIFF=$(wc -l < "$RESULTS_DIR/mask_fastq_only_seqs.txt")

echo "Found $BBMASK_DIFF sequences with different BBMask output"
echo "Found $MASK_FASTQ_DIFF sequences with different mask_fastq output"

if [ "$BBMASK_DIFF" -eq 0 ] && [ "$MASK_FASTQ_DIFF" -eq 0 ]; then
    echo "No differences found! Outputs are identical."
    exit 0
fi

# Create a temporary file with unique differing sequences (take first occurrence from BBMask diff)
head -n "$MAX_SEQS" "$RESULTS_DIR/bbmask_only_seqs.txt" > "$RESULTS_DIR/diff_seqs_to_find.txt"

echo "Extracting up to $MAX_SEQS differing sequences from original FASTQ..."

# Extract matching FASTQ records
# This script reads the original FASTQ and outputs records where the sequence matches our diff list
zcat -f "$INPUT_FASTQ" | awk '
    BEGIN {
        # Read sequences to find into array
        while ((getline line < "'"$RESULTS_DIR"'/diff_seqs_to_find.txt") > 0) {
            target_seqs[line] = 1
            count++
        }
        close("'"$RESULTS_DIR"'/diff_seqs_to_find.txt")
        found = 0
        max_found = '"$MAX_SEQS"'
    }
    NR % 4 == 1 {
        # Save header
        header = $0
    }
    NR % 4 == 2 {
        # Check if this sequence is in our target list
        if ($0 in target_seqs && found < max_found) {
            # Print this record
            print header
            print $0
            getline
            print $0
            getline
            print $0
            found++
            if (found >= max_found) exit
        }
    }
' > diff_sequences.fastq

EXTRACTED=$(grep -c "^@" diff_sequences.fastq || echo 0)

echo ""
echo "✓ Created diff_sequences.fastq with $EXTRACTED sequences"
echo ""
echo "Now you can run the benchmark on just these sequences:"
echo "  ./run_benchmark.sh diff_sequences.fastq"
echo ""
echo "Or examine the sequences manually to understand the differences:"
echo "  less diff_sequences.fastq"
echo ""

# Create a summary file showing the differences
echo "Creating summary of differences..."
{
    echo "# Sequences where BBMask and mask_fastq disagree"
    echo "# Original file: $INPUT_FASTQ"
    echo "# Total differences: BBMask=$BBMASK_DIFF, mask_fastq=$MASK_FASTQ_DIFF"
    echo ""

    # Show first few examples
    echo "## Example BBMask output (first 5):"
    head -n 5 "$RESULTS_DIR/bbmask_only_seqs.txt" | while read seq; do
        echo "Length: ${#seq} bases"
        echo "$seq" | head -c 100
        if [ ${#seq} -gt 100 ]; then echo "..."; else echo ""; fi
        echo ""
    done

    echo "## Example mask_fastq output (first 5):"
    head -n 5 "$RESULTS_DIR/mask_fastq_only_seqs.txt" | while read seq; do
        echo "Length: ${#seq} bases"
        echo "$seq" | head -c 100
        if [ ${#seq} -gt 100 ]; then echo "..."; else echo ""; fi
        echo ""
    done
} > diff_summary.txt

echo "✓ Created diff_summary.txt with detailed comparison"
echo ""
