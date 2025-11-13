#!/bin/bash
# Benchmark comparison between bbmask.sh and rustmasker
#
# Usage:
#   ./benchmark_vs_bbmask.sh test.fastq
#   ./benchmark_vs_bbmask.sh test.fastq --window 25 --entropy 0.55 --kmer 5

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default parameters (matching bbmask defaults)
WINDOW=80
ENTROPY=0.70
KMER=5
THREADS=4  # Default thread count for parallel version

# Parse command line arguments
INPUT_FASTQ=""
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            echo "Usage: $0 <input.fastq> [options]"
            echo ""
            echo "Options:"
            echo "  -w, --window N    Window size for entropy calculation (default: $WINDOW)"
            echo "  -e, --entropy N   Entropy threshold (default: $ENTROPY)"
            echo "  -k, --kmer N      K-mer size (default: $KMER)"
            echo "  -t, --threads N   Number of threads for parallel version (default: $THREADS)"
            echo "  -h, --help        Show this help"
            exit 0
            ;;
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
        -t|--threads)
            THREADS="$2"
            shift 2
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

# Find the rustmasker binary
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
RUSTMASKER="$PROJECT_DIR/rustmasker/target/release/rustmasker"

if [[ ! -f "$RUSTMASKER" ]]; then
    echo -e "${RED}Error: rustmasker binary not found at: $RUSTMASKER${NC}" >&2
    echo "Please build it first:" >&2
    echo "  cd $PROJECT_DIR/rustmasker" >&2
    echo "  cargo build --release" >&2
    exit 1
fi

# Check if bbmask.sh is available
if ! command -v bbmask.sh &> /dev/null; then
    echo -e "${RED}Error: bbmask.sh not found in PATH${NC}" >&2
    exit 1
fi

# Check if GNU time is available for memory stats. On macOS, look for
# Homebrew's gtime.
OS="$(uname -s)"

if [ "$OS" = "Darwin" ] && command -v gtime >/dev/null 2>&1; then
    TIME_CMD="$(command -v gtime) -v"
    HAS_MEMORY=1
elif [ "$OS" != "Darwin" ] && command -v /usr/bin/time >/dev/null 2>&1; then
    TIME_CMD="/usr/bin/time -v"
    HAS_MEMORY=1
else
    TIME_CMD="time"
    HAS_MEMORY=0
    echo -e "${YELLOW}Note: Install GNU Time for memory statistics.${NC}"
    if [ "$OS" = "Darwin" ]; then
        echo "  macOS: brew install gnu-time   # then use 'gtime -v'"
    else
        echo "  Linux: sudo apt-get install time"
    fi
    echo ""
fi

# Check if we need to use gzcat (macOS) or zcat (Linux) for decompressing gzip files
if [ "$OS" = "Darwin" ]; then
    ZCAT_CMD="gzcat"
else
    ZCAT_CMD="zcat"
fi

# Create output directory
OUTPUT_DIR="$SCRIPT_DIR/results"
mkdir -p "$OUTPUT_DIR"

# Get input file stats
INPUT_SIZE=$(du -h "$INPUT_FASTQ" | cut -f1)
NUM_READS=$(grep -c "^@" "$INPUT_FASTQ" || echo "unknown")

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║           Benchmark: bbmask.sh vs rustmasker                 ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${GREEN}Input file:${NC} $INPUT_FASTQ"
echo -e "${GREEN}File size:${NC} $INPUT_SIZE"
echo -e "${GREEN}Num reads:${NC} $NUM_READS"
echo -e "${GREEN}Parameters:${NC} window=$WINDOW, entropy=$ENTROPY, ke=$KMER"
echo -e "${GREEN}Threads:${NC} $THREADS"
echo -e "${GREEN}BBMask mode:${NC} entropy-only (maskrepeats=f, masklowentropy=t)"
echo ""

# Benchmark bbmask.sh
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}Running bbmask.sh...${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

BBMASK_OUT="$OUTPUT_DIR/bbmask_output.fastq.gz"
BBMASK_LOG="$OUTPUT_DIR/bbmask_log.txt"
BBMASK_TIME="$OUTPUT_DIR/bbmask_time.txt"

if [[ "$HAS_MEMORY" -eq 1 ]]; then
    $TIME_CMD bbmask.sh \
        in="$INPUT_FASTQ" \
        out="$BBMASK_OUT" \
        entropy="$ENTROPY" \
        ke="$KMER" \
        window="$WINDOW" \
        maskrepeats=f \
        masklowentropy=t \
        2>&1 | tee "$BBMASK_LOG" | tee "$BBMASK_TIME"
else
    { time bbmask.sh \
        in="$INPUT_FASTQ" \
        out="$BBMASK_OUT" \
        entropy="$ENTROPY" \
        ke="$KMER" \
        window="$WINDOW" \
        maskrepeats=f \
        masklowentropy=t \
        > "$BBMASK_LOG" 2>&1 ; } 2> "$BBMASK_TIME"
fi

# Extract bbmask stats
BBMASK_RUNTIME=$(grep "Total Time:" "$BBMASK_LOG" | awk '{print $3}' || echo "unknown")
BBMASK_MASKED=$(grep "Total Bases Masked:" "$BBMASK_LOG" | awk '{print $4}' || echo "unknown")

# Calculate BBMask masked percentage
if [[ "$BBMASK_MASKED" != "unknown" && "$BBMASK_MASKED" == *"/"* ]]; then
    BBMASK_MASKED_BASES=$(echo "$BBMASK_MASKED" | cut -d'/' -f1)
    BBMASK_TOTAL_BASES=$(echo "$BBMASK_MASKED" | cut -d'/' -f2)
    BBMASK_MASKED_PCT=$(echo "scale=3; $BBMASK_MASKED_BASES * 100 / $BBMASK_TOTAL_BASES" | bc)
else
    BBMASK_MASKED_PCT="-"
fi

if [[ "$HAS_MEMORY" -eq 1 ]]; then
    BBMASK_MEMORY=$(grep "Maximum resident set size" "$BBMASK_TIME" | awk '{print $6}' || echo "unknown")
    BBMASK_MEMORY_MB=$(echo "scale=2; $BBMASK_MEMORY / 1024" | bc)
else
    BBMASK_MEMORY_MB="unavailable"
fi

echo -e "${GREEN}✓ bbmask.sh completed${NC}"
echo ""

# Benchmark rustmasker
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}Running rustmasker ($THREADS threads)...${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

RUSTMASKER_OUT="$OUTPUT_DIR/rustmasker_output.fastq.gz"
RUSTMASKER_TIME="$OUTPUT_DIR/rustmasker_time.txt"

if [[ "$HAS_MEMORY" -eq 1 ]]; then
    $TIME_CMD sh -c "cat '$INPUT_FASTQ' | '$RUSTMASKER' -w $WINDOW -e $ENTROPY -k $KMER -c 1 -t $THREADS > '$RUSTMASKER_OUT'" 2>&1 | tee "$RUSTMASKER_TIME"
else
    { time sh -c "cat '$INPUT_FASTQ' | '$RUSTMASKER' -w $WINDOW -e $ENTROPY -k $KMER -c 1 -t $THREADS > '$RUSTMASKER_OUT'" ; } 2> "$RUSTMASKER_TIME"
fi

# Extract rustmasker stats
if [[ "$HAS_MEMORY" -eq 1 ]]; then
    RUSTMASKER_RUNTIME_RAW=$(grep "Elapsed (wall clock)" "$RUSTMASKER_TIME" | awk '{print $8}')
    # Convert h:mm:ss or m:ss format to seconds
    if [[ "$RUSTMASKER_RUNTIME_RAW" == *":"*":"* ]]; then
        IFS=':' read -r hours minutes seconds <<< "$RUSTMASKER_RUNTIME_RAW"
        RUSTMASKER_RUNTIME=$(echo "$hours * 3600 + $minutes * 60 + $seconds" | bc)
    elif [[ "$RUSTMASKER_RUNTIME_RAW" == *":"* ]]; then
        IFS=':' read -r minutes seconds <<< "$RUSTMASKER_RUNTIME_RAW"
        RUSTMASKER_RUNTIME=$(echo "$minutes * 60 + $seconds" | bc)
    else
        RUSTMASKER_RUNTIME="$RUSTMASKER_RUNTIME_RAW"
    fi
    RUSTMASKER_MEMORY=$(grep "Maximum resident set size" "$RUSTMASKER_TIME" | awk '{print $6}')
    RUSTMASKER_MEMORY_MB=$(echo "scale=2; $RUSTMASKER_MEMORY / 1024" | bc)
else
    RUSTMASKER_RUNTIME=$(grep "real" "$RUSTMASKER_TIME" | awk '{print $2}')
    RUSTMASKER_MEMORY_MB="unavailable"
fi

# Count masked bases in rustmasker output
RUSTMASKER_TOTAL_BASES=$($ZCAT_CMD "$RUSTMASKER_OUT" | awk 'NR%4==2 {sum+=length($0)} END {print sum}')
RUSTMASKER_MASKED_BASES=$($ZCAT_CMD "$RUSTMASKER_OUT" | awk 'NR%4==2 {gsub(/[^N]/,"",$0); sum+=length($0)} END {print sum}')
RUSTMASKER_MASKED_PCT=$(echo "scale=3; $RUSTMASKER_MASKED_BASES * 100 / $RUSTMASKER_TOTAL_BASES" | bc)

echo -e "${GREEN}✓ rustmasker completed${NC}"
echo ""

# Compare outputs
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}Comparing outputs...${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

# Extract sequences only for comparison
$ZCAT_CMD "$BBMASK_OUT" | awk 'NR%4==2' | sort > "$OUTPUT_DIR/bbmask_seqs.txt"
$ZCAT_CMD "$RUSTMASKER_OUT" | awk 'NR%4==2' | sort > "$OUTPUT_DIR/rustmasker_seqs.txt"

# Compare outputs
OUTPUTS_MATCH="YES"
if diff -q "$OUTPUT_DIR/bbmask_seqs.txt" "$OUTPUT_DIR/rustmasker_seqs.txt" > /dev/null; then
    echo -e "${GREEN}✓ BBMask vs rustmasker: IDENTICAL${NC}"
else
    echo -e "${RED}✗ BBMask vs rustmasker: DIFFER${NC}"
    OUTPUTS_MATCH="NO"
fi
echo ""

# Print results table
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║                              BENCHMARK RESULTS                                   ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Format runtimes with proper decimal places
BBMASK_RUNTIME_FMT=$(printf "%.3f" "$BBMASK_RUNTIME" 2>/dev/null || echo "$BBMASK_RUNTIME")
RUSTMASKER_RUNTIME_FMT=$(printf "%.3f" "$RUSTMASKER_RUNTIME" 2>/dev/null || echo "$RUSTMASKER_RUNTIME")

printf "%-25s %-20s %-20s\n" "Metric" "bbmask.sh" "rustmasker"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
printf "%-25s %-20s %-20s\n" "Runtime" "${BBMASK_RUNTIME_FMT}s" "${RUSTMASKER_RUNTIME_FMT}s"
printf "%-25s %-20s %-20s\n" "Peak Memory (MB)" "${BBMASK_MEMORY_MB}" "${RUSTMASKER_MEMORY_MB}"
printf "%-25s %-20s %-20s\n" "Masked bases" "$BBMASK_MASKED" "${RUSTMASKER_MASKED_BASES}/${RUSTMASKER_TOTAL_BASES}"
printf "%-25s %-20s %-20s\n" "Masked percentage" "${BBMASK_MASKED_PCT}%" "${RUSTMASKER_MASKED_PCT}%"
printf "%-25s %-20s\n" "Outputs match" "$OUTPUTS_MATCH"
echo ""

# Calculate speedup
if [[ "$BBMASK_RUNTIME" != "unknown" && "$RUSTMASKER_RUNTIME" != "unknown" ]]; then
    BBMASK_SEC=$(echo "$BBMASK_RUNTIME" | bc)
    RUSTMASKER_SEC=$(echo "$RUSTMASKER_RUNTIME" | bc)

    if [[ $(echo "$BBMASK_SEC > 0" | bc) -eq 1 ]]; then
        SPEEDUP=$(echo "scale=2; $BBMASK_SEC / $RUSTMASKER_SEC" | bc)
        echo -e "${GREEN}Speedup (rustmasker vs BBMask): ${SPEEDUP}x${NC}"
    fi
fi

# Calculate memory reduction
if [[ "$BBMASK_MEMORY_MB" != "unavailable" && "$RUSTMASKER_MEMORY_MB" != "unavailable" ]]; then
    MEMORY_REDUCTION=$(echo "scale=2; (1 - $RUSTMASKER_MEMORY_MB / $BBMASK_MEMORY_MB) * 100" | bc)
    echo -e "${GREEN}Memory reduction: ${MEMORY_REDUCTION}%${NC}"
fi

echo ""
echo -e "${BLUE}Results saved to: $OUTPUT_DIR/${NC}"
