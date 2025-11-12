#!/bin/bash
# Run comprehensive benchmarks comparing bbmask.sh and mask_fastq
#
# This script:
# 1. Generates test datasets of various sizes
# 2. Runs benchmarks on each dataset
# 3. Summarizes results

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}════════════════════════════════════════════════════════════${NC}"
echo -e "${BLUE}   Comprehensive Benchmark: bbmask.sh vs mask_fastq         ${NC}"
echo -e "${BLUE}════════════════════════════════════════════════════════════${NC}"
echo ""

# Make scripts executable
chmod +x "$SCRIPT_DIR/generate_test_data.py"
chmod +x "$SCRIPT_DIR/run_benchmark.sh"

# Create test data directory
TEST_DATA_DIR="$SCRIPT_DIR/test_data"
mkdir -p "$TEST_DATA_DIR"

# Define test scenarios
echo -e "${GREEN}Generating test datasets...${NC}"
echo ""

# Scenario 1: Small dataset (Illumina-like, quick test)
echo -e "${YELLOW}[1/4] Small dataset: 1K reads, 150bp (Illumina-like)${NC}"
python3 "$SCRIPT_DIR/generate_test_data.py" \
    --num-reads 1000 \
    --read-length 150 \
    -o "$TEST_DATA_DIR/small_illumina.fastq"
echo ""

# Scenario 2: Medium dataset (Illumina-like)
echo -e "${YELLOW}[2/4] Medium dataset: 10K reads, 150bp${NC}"
python3 "$SCRIPT_DIR/generate_test_data.py" \
    --num-reads 10000 \
    --read-length 150 \
    -o "$TEST_DATA_DIR/medium_illumina.fastq"
echo ""

# Scenario 3: Long reads (ONT-like, this is where memory issues occur)
echo -e "${YELLOW}[3/4] Long reads: 1K reads, 10Kbp (ONT-like)${NC}"
python3 "$SCRIPT_DIR/generate_test_data.py" \
    --num-reads 1000 \
    --read-length 10000 \
    -o "$TEST_DATA_DIR/long_ont.fastq"
echo ""

# Scenario 4: Very long reads (ONT ultra-long, stress test)
echo -e "${YELLOW}[4/4] Very long reads: 100 reads, 100Kbp (ONT ultra-long)${NC}"
python3 "$SCRIPT_DIR/generate_test_data.py" \
    --num-reads 100 \
    --read-length 100000 \
    -o "$TEST_DATA_DIR/ultralong_ont.fastq"
echo ""

# Run benchmarks
echo -e "${BLUE}════════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}Running benchmarks...${NC}"
echo -e "${BLUE}════════════════════════════════════════════════════════════${NC}"
echo ""

# Scenario 1
echo -e "${YELLOW}Benchmark 1: Small Illumina dataset${NC}"
"$SCRIPT_DIR/run_benchmark.sh" "$TEST_DATA_DIR/small_illumina.fastq"
echo ""

# Scenario 2
echo -e "${YELLOW}Benchmark 2: Medium Illumina dataset${NC}"
"$SCRIPT_DIR/run_benchmark.sh" "$TEST_DATA_DIR/medium_illumina.fastq"
echo ""

# Scenario 3
echo -e "${YELLOW}Benchmark 3: Long ONT reads${NC}"
"$SCRIPT_DIR/run_benchmark.sh" "$TEST_DATA_DIR/long_ont.fastq"
echo ""

# Scenario 4
echo -e "${YELLOW}Benchmark 4: Ultra-long ONT reads (this is the stress test!)${NC}"
"$SCRIPT_DIR/run_benchmark.sh" "$TEST_DATA_DIR/ultralong_ont.fastq"
echo ""

echo -e "${GREEN}════════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}All benchmarks complete!${NC}"
echo -e "${GREEN}════════════════════════════════════════════════════════════${NC}"
echo ""
echo "Results saved to: $SCRIPT_DIR/results/"
echo ""
echo "Key files:"
echo "  - bbmask_log.txt      : bbmask.sh output"
echo "  - bbmask_time.txt     : bbmask.sh timing/memory"
echo "  - mask_fastq_time.txt : mask_fastq timing/memory"
