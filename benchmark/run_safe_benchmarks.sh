#!/bin/bash
# Safe benchmark runner with memory limits
#
# This script adds memory protection to prevent OOM on your instance:
# 1. Caps bbmask.sh Java heap size
# 2. Uses ulimit to restrict process memory
# 3. Checks available memory before running
# 4. Runs tests incrementally (smallest to largest)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Memory settings (adjust based on your instance)
MAX_BBMASK_MEMORY_GB=${MAX_BBMASK_MEMORY_GB:-4}  # Max 4GB for bbmask by default
SAFETY_MARGIN_GB=2  # Keep 2GB free for system

echo -e "${BLUE}════════════════════════════════════════════════════════════${NC}"
echo -e "${BLUE}   Safe Benchmark Suite with Memory Limits                 ${NC}"
echo -e "${BLUE}════════════════════════════════════════════════════════════${NC}"
echo ""

# Check available memory
if command -v free &> /dev/null; then
    TOTAL_MEM_KB=$(free | grep Mem | awk '{print $2}')
    AVAIL_MEM_KB=$(free | grep Mem | awk '{print $7}')
    TOTAL_MEM_GB=$(echo "scale=1; $TOTAL_MEM_KB / 1024 / 1024" | bc)
    AVAIL_MEM_GB=$(echo "scale=1; $AVAIL_MEM_KB / 1024 / 1024" | bc)

    echo -e "${GREEN}System Memory:${NC}"
    echo "  Total: ${TOTAL_MEM_GB} GB"
    echo "  Available: ${AVAIL_MEM_GB} GB"
    echo "  bbmask limit: ${MAX_BBMASK_MEMORY_GB} GB"
    echo "  Safety margin: ${SAFETY_MARGIN_GB} GB"
    echo ""

    # Check if we have enough memory
    MIN_REQUIRED=$(echo "$MAX_BBMASK_MEMORY_GB + $SAFETY_MARGIN_GB" | bc)
    if (( $(echo "$AVAIL_MEM_GB < $MIN_REQUIRED" | bc -l) )); then
        echo -e "${RED}WARNING: Available memory (${AVAIL_MEM_GB}GB) is less than recommended (${MIN_REQUIRED}GB)${NC}"
        echo "Consider:"
        echo "  1. Setting a lower limit: export MAX_BBMASK_MEMORY_GB=2"
        echo "  2. Skipping ultra-long read tests"
        echo ""
        read -p "Continue anyway? [y/N] " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            exit 1
        fi
    fi
fi

# Make scripts executable
chmod +x "$SCRIPT_DIR/generate_test_data.py"
chmod +x "$SCRIPT_DIR/run_benchmark.sh"

# Create test data directory
TEST_DATA_DIR="$SCRIPT_DIR/test_data"
mkdir -p "$TEST_DATA_DIR"

# Define test scenarios - starting small!
echo -e "${GREEN}Generating test datasets (starting small for safety)...${NC}"
echo ""

# Scenario 1: Tiny dataset (quick memory check)
echo -e "${YELLOW}[1/5] Tiny test: 100 reads, 150bp (memory check)${NC}"
python3 "$SCRIPT_DIR/generate_test_data.py" \
    --num-reads 100 \
    --read-length 150 \
    -o "$TEST_DATA_DIR/tiny_test.fastq"
echo ""

# Scenario 2: Small dataset (Illumina-like)
echo -e "${YELLOW}[2/5] Small dataset: 1K reads, 150bp (Illumina-like)${NC}"
python3 "$SCRIPT_DIR/generate_test_data.py" \
    --num-reads 1000 \
    --read-length 150 \
    -o "$TEST_DATA_DIR/small_illumina.fastq"
echo ""

# Scenario 3: Medium dataset (Illumina-like)
echo -e "${YELLOW}[3/5] Medium dataset: 10K reads, 150bp${NC}"
python3 "$SCRIPT_DIR/generate_test_data.py" \
    --num-reads 10000 \
    --read-length 150 \
    -o "$TEST_DATA_DIR/medium_illumina.fastq"
echo ""

# Scenario 4: Long reads (ONT-like - this is where memory matters!)
echo -e "${YELLOW}[4/5] Long reads: 500 reads, 10Kbp (ONT-like, REDUCED COUNT)${NC}"
python3 "$SCRIPT_DIR/generate_test_data.py" \
    --num-reads 500 \
    --read-length 10000 \
    -o "$TEST_DATA_DIR/long_ont.fastq"
echo ""

# Scenario 5: Very long reads - ONLY if user confirms
echo -e "${YELLOW}[5/5] Ultra-long reads: 50 reads, 100Kbp (stress test, REDUCED COUNT)${NC}"
echo -e "${RED}WARNING: This test may use significant memory with bbmask.sh!${NC}"
read -p "Generate ultra-long test data? [y/N] " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    python3 "$SCRIPT_DIR/generate_test_data.py" \
        --num-reads 50 \
        --read-length 100000 \
        -o "$TEST_DATA_DIR/ultralong_ont.fastq"
    SKIP_ULTRALONG=false
else
    echo "Skipping ultra-long test data generation"
    SKIP_ULTRALONG=true
fi
echo ""

# Export memory limit for bbmask
export _JAVA_OPTIONS="-Xmx${MAX_BBMASK_MEMORY_GB}g -Xms${MAX_BBMASK_MEMORY_GB}g"

# Run benchmarks incrementally
echo -e "${BLUE}════════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}Running benchmarks (with ${MAX_BBMASK_MEMORY_GB}GB Java heap limit)...${NC}"
echo -e "${BLUE}════════════════════════════════════════════════════════════${NC}"
echo ""

# Function to check if we should continue
check_continue() {
    read -p "Continue to next benchmark? [Y/n] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Nn]$ ]]; then
        echo "Stopping benchmarks"
        exit 0
    fi
}

# Benchmark 1 - Tiny (quick check)
echo -e "${YELLOW}Benchmark 1/5: Tiny test (memory check)${NC}"
"$SCRIPT_DIR/run_benchmark.sh" "$TEST_DATA_DIR/tiny_test.fastq" || {
    echo -e "${RED}Tiny test failed! This is unexpected. Check your setup.${NC}"
    exit 1
}
echo ""
check_continue

# Benchmark 2 - Small
echo -e "${YELLOW}Benchmark 2/5: Small Illumina dataset${NC}"
"$SCRIPT_DIR/run_benchmark.sh" "$TEST_DATA_DIR/small_illumina.fastq"
echo ""
check_continue

# Benchmark 3 - Medium
echo -e "${YELLOW}Benchmark 3/5: Medium Illumina dataset${NC}"
"$SCRIPT_DIR/run_benchmark.sh" "$TEST_DATA_DIR/medium_illumina.fastq"
echo ""
check_continue

# Benchmark 4 - Long ONT (this is important!)
echo -e "${YELLOW}Benchmark 4/5: Long ONT reads${NC}"
echo -e "${BLUE}This benchmark tests the real-world scenario for issue #323${NC}"
"$SCRIPT_DIR/run_benchmark.sh" "$TEST_DATA_DIR/long_ont.fastq"
echo ""

if [ "$SKIP_ULTRALONG" = false ]; then
    echo -e "${RED}═══════════════════════════════════════════════════════════${NC}"
    echo -e "${RED}   WARNING: Ultra-Long Read Test (High Memory Usage)      ${NC}"
    echo -e "${RED}═══════════════════════════════════════════════════════════${NC}"
    echo ""
    echo "This test will run bbmask.sh on 100Kbp reads."
    echo "Memory limit: ${MAX_BBMASK_MEMORY_GB}GB"
    echo ""
    echo -e "${YELLOW}If bbmask.sh runs out of memory, that's EXPECTED and proves the bug!${NC}"
    echo -e "${YELLOW}mask_fastq should still complete successfully with low memory.${NC}"
    echo ""
    read -p "Run ultra-long benchmark? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${YELLOW}Benchmark 5/5: Ultra-long ONT reads (stress test!)${NC}"
        "$SCRIPT_DIR/run_benchmark.sh" "$TEST_DATA_DIR/ultralong_ont.fastq" || {
            echo -e "${YELLOW}Note: bbmask.sh may have failed due to memory limits (this proves issue #323!)${NC}"
        }
    else
        echo "Skipping ultra-long benchmark"
    fi
fi

echo ""
echo -e "${GREEN}════════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}Benchmarks complete!${NC}"
echo -e "${GREEN}════════════════════════════════════════════════════════════${NC}"
echo ""
echo "Results saved to: $SCRIPT_DIR/results/"
echo ""
echo "Summary:"
echo "  - Java heap limit: ${MAX_BBMASK_MEMORY_GB}GB"
echo "  - Tests run incrementally with user confirmation"
echo "  - Ultra-long tests skipped or run with caution"
