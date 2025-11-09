#!/bin/bash
# Test script to demonstrate help output behavior

echo "=== Test 1: -h flag (shows clap's auto-generated help) ==="
./target/release/mask_fastq -h
echo

echo "=== Test 2: --help flag (shows clap's auto-generated help) ==="
./target/release/mask_fastq_parallel --help
echo

echo "=== Test 3: Normal operation with piped input ==="
echo '@test
ACGTACGT
+
IIIIIIII' | ./target/release/mask_fastq -w 5 -e 0.5 -k 3 | zcat
echo

echo "=== Test 4: Run without arguments in interactive terminal ==="
echo "When you run './target/release/mask_fastq' directly in your terminal"
echo "(not piped, not in a script), you should see a helpful error message."
echo ""
echo "Try it yourself: ./target/release/mask_fastq"
