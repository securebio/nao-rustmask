#!/bin/bash

# Analyze patterns in the differences between BBMask and mask_fastq
# Usage: ./analyze_diff_patterns.sh
#
# Requires: diff_sequences.fastq (created by extract_diff_sequences.sh)

set -e

if [ ! -f "diff_sequences.fastq" ]; then
    echo "Error: diff_sequences.fastq not found"
    echo "Please run extract_diff_sequences.sh first"
    exit 1
fi

RESULTS_DIR="results"

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Analyzing masking differences between BBMask and mask_fastq"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Process each sequence in diff_sequences.fastq
awk 'NR % 4 == 1 {header = $0}
     NR % 4 == 2 {seq = $0; seqs[header] = seq}
     END {
         for (h in seqs) {
             print h
             print seqs[h]
         }
     }' diff_sequences.fastq > /tmp/original_seqs.txt

# Get BBMask and mask_fastq outputs for the diff sequences
# We need to run both tools on diff_sequences.fastq to get their outputs

echo "Running BBMask on diff sequences..."
if command -v bbmask.sh &> /dev/null; then
    bbmask.sh in=diff_sequences.fastq out=/tmp/bbmask_diff.fastq.gz \
        entropy=0.55 ke=5 window=25 maskrepeats=f masklowentropy=t 2>&1 | grep -v "^java" | head -5

    echo ""
    echo "Running mask_fastq on diff sequences..."
    cat diff_sequences.fastq | ../modules/local/maskRead/mask_fastq/target/release/mask_fastq \
        -w 25 -e 0.55 -k 5 > /tmp/mask_fastq_diff.fastq.gz

    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Detailed comparison of each sequence:"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""

    # Extract sequences from both outputs
    zcat /tmp/bbmask_diff.fastq.gz | awk 'NR % 4 == 1 {id = substr($0, 2)}
                                           NR % 4 == 2 {print id "\t" $0}' \
        > /tmp/bbmask_output.tsv

    zcat /tmp/mask_fastq_diff.fastq.gz | awk 'NR % 4 == 1 {id = substr($0, 2)}
                                               NR % 4 == 2 {print id "\t" $0}' \
        > /tmp/mask_fastq_output.tsv

    # Compare each sequence
    awk 'NR % 2 == 1 {id = substr($0, 2)}
         NR % 2 == 0 {original[id] = $0}' /tmp/original_seqs.txt > /tmp/originals.tmp

    awk '{bbmask[$1] = $2}' /tmp/bbmask_output.tsv > /tmp/bbmask.tmp
    awk '{mask_fastq[$1] = $2}' /tmp/mask_fastq_output.tsv > /tmp/mask_fastq.tmp

    # Show differences
    paste /tmp/bbmask_output.tsv /tmp/mask_fastq_output.tsv | head -20 | while IFS=$'\t' read id1 bbmask_seq id2 mask_fastq_seq; do
        echo "Sequence: $id1"

        # Count N's
        bbmask_n=$(echo "$bbmask_seq" | tr -cd 'N' | wc -c)
        mask_fastq_n=$(echo "$mask_fastq_seq" | tr -cd 'N' | wc -c)
        seq_len=${#bbmask_seq}

        bbmask_pct=$(awk "BEGIN {printf \"%.1f\", $bbmask_n * 100.0 / $seq_len}")
        mask_fastq_pct=$(awk "BEGIN {printf \"%.1f\", $mask_fastq_n * 100.0 / $seq_len}")

        echo "  Length: $seq_len bp"
        echo "  BBMask masked:     $bbmask_n / $seq_len ($bbmask_pct%)"
        echo "  mask_fastq masked: $mask_fastq_n / $seq_len ($mask_fastq_pct%)"

        # Show first and last 50 bases comparison
        echo "  First 50 bases:"
        echo "    BBMask:     ${bbmask_seq:0:50}"
        echo "    mask_fastq: ${mask_fastq_seq:0:50}"

        if [ $seq_len -gt 50 ]; then
            echo "  Last 50 bases:"
            echo "    BBMask:     ${bbmask_seq: -50}"
            echo "    mask_fastq: ${mask_fastq_seq: -50}"
        fi

        echo ""
    done

    # Clean up
    rm -f /tmp/bbmask_diff.fastq.gz /tmp/mask_fastq_diff.fastq.gz
    rm -f /tmp/bbmask_output.tsv /tmp/mask_fastq_output.tsv
    rm -f /tmp/originals.tmp /tmp/bbmask.tmp /tmp/mask_fastq.tmp
else
    echo "Warning: bbmask.sh not found. Showing original sequences only."
    echo ""
    head -40 /tmp/original_seqs.txt | while IFS= read -r line; do
        if [[ $line == @* ]]; then
            echo ""
            echo "Sequence: ${line:1}"
        else
            echo "  Length: ${#line} bp"
            echo "  First 80 bp: ${line:0:80}"
            if [ ${#line} -gt 80 ]; then
                echo "  Last 80 bp:  ${line: -80}"
            fi
        fi
    done
fi

rm -f /tmp/original_seqs.txt

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Analysis complete!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
