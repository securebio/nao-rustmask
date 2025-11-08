#!/bin/bash
# Benchmark comparison between bbmask.sh and mask_fastq
#
# Usage:
#   ./run_benchmark.sh test.fastq
#   ./run_benchmark.sh test.fastq --window 25 --entropy 0.55 --kmer 5

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default parameters (matching bbmask defaults)
WINDOW=25
ENTROPY=0.55
KMER=5

# Parse command line arguments
INPUT_FASTQ=""
while [[ $# -gt 0 ]]; do
    case $1 in
        -w|--window)
            WINDOW="$2"
            shift 2
            ;;
        -e|--entropy)
            ENTROPY="$2"
            shift 2
            ;;
        -k|--kmer)
            KMER="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 <input.fastq> [options]"
            echo ""
            echo "Options:"
            echo "  -w, --window N    Window size for entropy calculation (default: 25)"
            echo "  -e, --entropy N   Entropy threshold (default: 0.55)"
            echo "  -k, --kmer N      K-mer size (default: 5)"
            echo "  -h, --help        Show this help"
            exit 0
            ;;
        *)
            if [[ -z "$INPUT_FASTQ" ]]; then
                INPUT_FASTQ="$1"
            else
                echo "Error: Unknown argument: $1" >&2
                exit 1
            fi
            shift
            ;;
    esac
done

if [[ -z "$INPUT_FASTQ" ]]; then
    echo "Error: No input FASTQ file specified" >&2
    echo "Usage: $0 <input.fastq>" >&2
    exit 1
fi

if [[ ! -f "$INPUT_FASTQ" ]]; then
    echo "Error: Input file not found: $INPUT_FASTQ" >&2
    exit 1
fi

# Find the mask_fastq binary
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
MASK_FASTQ="$PROJECT_DIR/modules/local/maskRead/mask_fastq/target/release/mask_fastq"

if [[ ! -f "$MASK_FASTQ" ]]; then
    echo -e "${RED}Error: mask_fastq binary not found at: $MASK_FASTQ${NC}" >&2
    echo "Please build it first:" >&2
    echo "  cd $PROJECT_DIR/modules/local/maskRead/mask_fastq" >&2
    echo "  cargo build --release" >&2
    exit 1
fi

# Check if bbmask.sh is available
if ! command -v bbmask.sh &> /dev/null; then
    echo -e "${RED}Error: bbmask.sh not found in PATH${NC}" >&2
    exit 1
fi

# Check if /usr/bin/time is available (for detailed memory stats)
TIME_CMD="/usr/bin/time"
if [[ ! -x "$TIME_CMD" ]]; then
    # Fall back to bash built-in time
    TIME_CMD="time"
    echo -e "${YELLOW}Warning: /usr/bin/time not found, memory stats will be limited${NC}"
fi

# Create output directory
OUTPUT_DIR="$SCRIPT_DIR/results"
mkdir -p "$OUTPUT_DIR"

# Get input file stats
INPUT_SIZE=$(du -h "$INPUT_FASTQ" | cut -f1)
NUM_READS=$(grep -c "^@" "$INPUT_FASTQ" || echo "unknown")

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║           Benchmark: bbmask.sh vs mask_fastq                 ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${GREEN}Input file:${NC} $INPUT_FASTQ"
echo -e "${GREEN}File size:${NC} $INPUT_SIZE"
echo -e "${GREEN}Num reads:${NC} $NUM_READS"
echo -e "${GREEN}Parameters:${NC} window=$WINDOW, entropy=$ENTROPY, k=$KMER"
echo ""

# Benchmark bbmask.sh
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}Running bbmask.sh...${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

BBMASK_OUT="$OUTPUT_DIR/bbmask_output.fastq.gz"
BBMASK_LOG="$OUTPUT_DIR/bbmask_log.txt"
BBMASK_TIME="$OUTPUT_DIR/bbmask_time.txt"

if [[ "$TIME_CMD" == "/usr/bin/time" ]]; then
    /usr/bin/time -v bbmask.sh \
        in="$INPUT_FASTQ" \
        out="$BBMASK_OUT" \
        entropy="$ENTROPY" \
        k="$KMER" \
        window="$WINDOW" \
        2>&1 | tee "$BBMASK_LOG" | tee "$BBMASK_TIME"
else
    { time bbmask.sh \
        in="$INPUT_FASTQ" \
        out="$BBMASK_OUT" \
        entropy="$ENTROPY" \
        k="$KMER" \
        window="$WINDOW" \
        > "$BBMASK_LOG" 2>&1 ; } 2> "$BBMASK_TIME"
fi

# Extract bbmask stats
BBMASK_RUNTIME=$(grep "Total Time:" "$BBMASK_LOG" | awk '{print $3}' || echo "unknown")
BBMASK_MASKED=$(grep "Total Bases Masked:" "$BBMASK_LOG" | awk '{print $4}' || echo "unknown")
if [[ "$TIME_CMD" == "/usr/bin/time" ]]; then
    BBMASK_MEMORY=$(grep "Maximum resident set size" "$BBMASK_TIME" | awk '{print $6}' || echo "unknown")
    BBMASK_MEMORY_MB=$(echo "scale=2; $BBMASK_MEMORY / 1024" | bc)
else
    BBMASK_MEMORY_MB="unavailable"
fi

echo -e "${GREEN}✓ bbmask.sh completed${NC}"
echo ""

# Benchmark mask_fastq
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}Running mask_fastq...${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

MASK_FASTQ_OUT="$OUTPUT_DIR/mask_fastq_output.fastq.gz"
MASK_FASTQ_TIME="$OUTPUT_DIR/mask_fastq_time.txt"

if [[ "$TIME_CMD" == "/usr/bin/time" ]]; then
    /usr/bin/time -v sh -c "cat '$INPUT_FASTQ' | '$MASK_FASTQ' -w $WINDOW -e $ENTROPY -k $KMER > '$MASK_FASTQ_OUT'" 2>&1 | tee "$MASK_FASTQ_TIME"
else
    { time sh -c "cat '$INPUT_FASTQ' | '$MASK_FASTQ' -w $WINDOW -e $ENTROPY -k $KMER > '$MASK_FASTQ_OUT'" ; } 2> "$MASK_FASTQ_TIME"
fi

# Extract mask_fastq stats
if [[ "$TIME_CMD" == "/usr/bin/time" ]]; then
    MASK_FASTQ_RUNTIME=$(grep "Elapsed (wall clock)" "$MASK_FASTQ_TIME" | awk '{print $8}')
    MASK_FASTQ_MEMORY=$(grep "Maximum resident set size" "$MASK_FASTQ_TIME" | awk '{print $6}')
    MASK_FASTQ_MEMORY_MB=$(echo "scale=2; $MASK_FASTQ_MEMORY / 1024" | bc)
else
    MASK_FASTQ_RUNTIME=$(grep "real" "$MASK_FASTQ_TIME" | awk '{print $2}')
    MASK_FASTQ_MEMORY_MB="unavailable"
fi

# Count masked bases in mask_fastq output
MASK_FASTQ_TOTAL_BASES=$(zcat "$MASK_FASTQ_OUT" | awk 'NR%4==2 {sum+=length($0)} END {print sum}')
MASK_FASTQ_MASKED_BASES=$(zcat "$MASK_FASTQ_OUT" | awk 'NR%4==2 {gsub(/[^N]/,"",$0); sum+=length($0)} END {print sum}')
MASK_FASTQ_MASKED_PCT=$(echo "scale=3; $MASK_FASTQ_MASKED_BASES * 100 / $MASK_FASTQ_TOTAL_BASES" | bc)

echo -e "${GREEN}✓ mask_fastq completed${NC}"
echo ""

# Compare outputs
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}Comparing outputs...${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

# Extract sequences only for comparison
zcat "$BBMASK_OUT" | awk 'NR%4==2' | sort > "$OUTPUT_DIR/bbmask_seqs.txt"
zcat "$MASK_FASTQ_OUT" | awk 'NR%4==2' | sort > "$OUTPUT_DIR/mask_fastq_seqs.txt"

if diff -q "$OUTPUT_DIR/bbmask_seqs.txt" "$OUTPUT_DIR/mask_fastq_seqs.txt" > /dev/null; then
    echo -e "${GREEN}✓ Outputs are IDENTICAL${NC}"
    OUTPUTS_MATCH="YES"
else
    echo -e "${RED}✗ Outputs DIFFER${NC}"
    OUTPUTS_MATCH="NO"
    echo "  Run: diff $OUTPUT_DIR/bbmask_seqs.txt $OUTPUT_DIR/mask_fastq_seqs.txt"
fi
echo ""

# Print results table
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║                      BENCHMARK RESULTS                       ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
printf "%-25s %-20s %-20s\n" "Metric" "bbmask.sh" "mask_fastq"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
printf "%-25s %-20s %-20s\n" "Runtime" "${BBMASK_RUNTIME}s" "${MASK_FASTQ_RUNTIME}"
printf "%-25s %-20s %-20s\n" "Peak Memory (MB)" "${BBMASK_MEMORY_MB}" "${MASK_FASTQ_MEMORY_MB}"
printf "%-25s %-20s %-20s\n" "Masked bases" "$BBMASK_MASKED" "${MASK_FASTQ_MASKED_BASES}/${MASK_FASTQ_TOTAL_BASES}"
printf "%-25s %-20s %-20s\n" "Masked percentage" "-" "${MASK_FASTQ_MASKED_PCT}%"
printf "%-25s %-20s\n" "Outputs match" "$OUTPUTS_MATCH"
echo ""

# Calculate speedup
if [[ "$BBMASK_RUNTIME" != "unknown" && "$MASK_FASTQ_RUNTIME" != *":"* ]]; then
    # Convert time to seconds if needed
    BBMASK_SEC=$(echo "$BBMASK_RUNTIME" | bc)
    MASK_FASTQ_SEC=$(echo "$MASK_FASTQ_RUNTIME" | sed 's/[^0-9.]//g' | head -c 10)

    if [[ -n "$MASK_FASTQ_SEC" && $(echo "$MASK_FASTQ_SEC > 0" | bc) -eq 1 ]]; then
        SPEEDUP=$(echo "scale=2; $BBMASK_SEC / $MASK_FASTQ_SEC" | bc)
        echo -e "${GREEN}Speedup: ${SPEEDUP}x faster${NC}"
    fi
fi

# Calculate memory reduction
if [[ "$BBMASK_MEMORY_MB" != "unavailable" && "$MASK_FASTQ_MEMORY_MB" != "unavailable" ]]; then
    MEMORY_REDUCTION=$(echo "scale=2; (1 - $MASK_FASTQ_MEMORY_MB / $BBMASK_MEMORY_MB) * 100" | bc)
    echo -e "${GREEN}Memory reduction: ${MEMORY_REDUCTION}%${NC}"
fi

echo ""
echo -e "${BLUE}Results saved to: $OUTPUT_DIR/${NC}"
