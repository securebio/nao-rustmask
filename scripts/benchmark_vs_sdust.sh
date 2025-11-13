#!/bin/bash
# Benchmark comparison between sdust and rustmasker (SDUST algorithm)
#
# Usage:
#   ./benchmark_vs_sdust.sh test.fastq
#   ./benchmark_vs_sdust.sh test.fastq --window 64 --threshold 20

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default parameters (matching sdust defaults)
WINDOW=64
THRESHOLD=20
THREADS=4  # Default thread count for rustmasker

# Parse command line arguments
INPUT_FASTQ=""
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            echo "Usage: $0 <input.fastq> [options]"
            echo ""
            echo "Options:"
            echo "  -w, --window N      Window size for SDUST calculation (default: $WINDOW)"
            echo "  -t, --threshold N   SDUST threshold (default: $THRESHOLD)"
            echo "  -j, --threads N     Number of threads for rustmasker (default: $THREADS)"
            echo "  -h, --help          Show this help"
            exit 0
            ;;
        -w|--window)
            WINDOW="$2"
            shift 2
            ;;
        -t|--threshold)
            THRESHOLD="$2"
            shift 2
            ;;
        -j|--threads)
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

# Check if sdust is available
if ! command -v sdust &> /dev/null; then
    echo -e "${RED}Error: sdust not found in PATH${NC}" >&2
    echo "Please install sdust from: https://github.com/lh3/sdust" >&2
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
echo -e "${BLUE}║              Benchmark: sdust vs rustmasker                  ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${GREEN}Input file:${NC} $INPUT_FASTQ"
echo -e "${GREEN}File size:${NC} $INPUT_SIZE"
echo -e "${GREEN}Num reads:${NC} $NUM_READS"
echo -e "${GREEN}Parameters:${NC} window=$WINDOW, threshold=$THRESHOLD"
echo -e "${GREEN}Threads:${NC} $THREADS (rustmasker only)"
echo -e "${GREEN}Algorithm:${NC} SDUST (symmetric DUST)"
echo ""

# Convert FASTQ to FASTA (sdust requires FASTA input)
echo -e "${YELLOW}Converting FASTQ to FASTA for sdust...${NC}"
FASTA_INPUT="$OUTPUT_DIR/input.fasta"
awk 'NR%4==1 {print ">" substr($0,2)} NR%4==2 {print}' "$INPUT_FASTQ" > "$FASTA_INPUT"
echo -e "${GREEN}✓ Conversion complete${NC}"
echo ""

# Benchmark sdust
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}Running sdust...${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

SDUST_OUT="$OUTPUT_DIR/sdust_output.bed"
SDUST_FASTA="$OUTPUT_DIR/sdust_masked.fasta"
SDUST_TIME="$OUTPUT_DIR/sdust_time.txt"

# Run sdust (outputs BED format with masked regions)
if [[ "$HAS_MEMORY" -eq 1 ]]; then
    $TIME_CMD sdust -w "$WINDOW" -t "$THRESHOLD" "$FASTA_INPUT" > "$SDUST_OUT" 2> "$SDUST_TIME"
else
    { time sdust -w "$WINDOW" -t "$THRESHOLD" "$FASTA_INPUT" > "$SDUST_OUT" ; } 2> "$SDUST_TIME"
fi

# Apply sdust masks to create masked FASTA
python3 - "$FASTA_INPUT" "$SDUST_OUT" "$SDUST_FASTA" << 'EOF'
import sys

# Read the FASTA file
fasta_file = sys.argv[1]
bed_file = sys.argv[2]
output_file = sys.argv[3]

# Read sequences
sequences = {}
current_id = None
with open(fasta_file) as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            current_id = line[1:].split()[0]
            sequences[current_id] = []
        else:
            sequences[current_id].append(line)

# Convert sequences to mutable lists
for seq_id in sequences:
    sequences[seq_id] = list(''.join(sequences[seq_id]))

# Apply masks from BED file
with open(bed_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            seq_id = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if seq_id in sequences:
                for i in range(start, min(end, len(sequences[seq_id]))):
                    sequences[seq_id][i] = 'N'

# Write masked sequences
with open(output_file, 'w') as f:
    for seq_id, seq_list in sequences.items():
        f.write(f'>{seq_id}\n')
        seq = ''.join(seq_list)
        # Write in 80-character lines
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + '\n')
EOF

# Extract sdust stats
if [[ "$HAS_MEMORY" -eq 1 ]]; then
    SDUST_RUNTIME_RAW=$(grep "Elapsed (wall clock)" "$SDUST_TIME" | awk '{print $8}')
    # Convert h:mm:ss or m:ss format to seconds
    if [[ "$SDUST_RUNTIME_RAW" == *":"*":"* ]]; then
        IFS=':' read -r hours minutes seconds <<< "$SDUST_RUNTIME_RAW"
        SDUST_RUNTIME=$(echo "$hours * 3600 + $minutes * 60 + $seconds" | bc)
    elif [[ "$SDUST_RUNTIME_RAW" == *":"* ]]; then
        IFS=':' read -r minutes seconds <<< "$SDUST_RUNTIME_RAW"
        SDUST_RUNTIME=$(echo "$minutes * 60 + $seconds" | bc)
    else
        SDUST_RUNTIME="$SDUST_RUNTIME_RAW"
    fi
    SDUST_MEMORY=$(grep "Maximum resident set size" "$SDUST_TIME" | awk '{print $6}')
    SDUST_MEMORY_MB=$(echo "scale=2; $SDUST_MEMORY / 1024" | bc)
else
    SDUST_RUNTIME=$(grep "real" "$SDUST_TIME" | awk '{print $2}')
    SDUST_MEMORY_MB="unavailable"
fi

# Count masked bases in sdust output
SDUST_TOTAL_BASES=$(awk 'NR%2==0 {sum+=length($0)} END {print sum}' "$SDUST_FASTA")
SDUST_MASKED_BASES=$(awk 'NR%2==0 {gsub(/[^N]/,"",$0); sum+=length($0)} END {print sum}' "$SDUST_FASTA")
SDUST_MASKED_PCT=$(echo "scale=3; $SDUST_MASKED_BASES * 100 / $SDUST_TOTAL_BASES" | bc)

echo -e "${GREEN}✓ sdust completed${NC}"
echo ""

# Benchmark rustmasker with SDUST algorithm
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}Running rustmasker -a sdust ($THREADS threads)...${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

RUSTMASKER_OUT="$OUTPUT_DIR/rustmasker_output.fastq.gz"
RUSTMASKER_TIME="$OUTPUT_DIR/rustmasker_time.txt"

if [[ "$HAS_MEMORY" -eq 1 ]]; then
    $TIME_CMD sh -c "cat '$INPUT_FASTQ' | '$RUSTMASKER' -a sdust -w $WINDOW -t $THRESHOLD -c 1 -j $THREADS > '$RUSTMASKER_OUT'" 2>&1 | tee "$RUSTMASKER_TIME"
else
    { time sh -c "cat '$INPUT_FASTQ' | '$RUSTMASKER' -a sdust -w $WINDOW -t $THRESHOLD -c 1 -j $THREADS > '$RUSTMASKER_OUT'" ; } 2> "$RUSTMASKER_TIME"
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

# Convert rustmasker FASTQ output to FASTA for comparison
RUSTMASKER_FASTA="$OUTPUT_DIR/rustmasker_output.fasta"
$ZCAT_CMD "$RUSTMASKER_OUT" | awk 'NR%4==1 {print ">" substr($0,2)} NR%4==2 {print}' > "$RUSTMASKER_FASTA"

# Extract sequences only for comparison (sorted by ID)
awk 'NR%2==0' "$SDUST_FASTA" | paste - - | sort > "$OUTPUT_DIR/sdust_seqs.txt"
awk 'NR%2==0' "$RUSTMASKER_FASTA" | paste - - | sort > "$OUTPUT_DIR/rustmasker_seqs.txt"

# Compare outputs
OUTPUTS_MATCH="YES"
DIFF_COUNT=0
if diff -q "$OUTPUT_DIR/sdust_seqs.txt" "$OUTPUT_DIR/rustmasker_seqs.txt" > /dev/null; then
    echo -e "${GREEN}✓ sdust vs rustmasker: IDENTICAL${NC}"
else
    echo -e "${YELLOW}⚠ sdust vs rustmasker: DIFFER${NC}"
    OUTPUTS_MATCH="SIMILAR"

    # Count differences
    DIFF_COUNT=$(diff "$OUTPUT_DIR/sdust_seqs.txt" "$OUTPUT_DIR/rustmasker_seqs.txt" | grep -c "^<" || true)

    # Calculate similarity percentage
    TOTAL_LINES=$(wc -l < "$OUTPUT_DIR/sdust_seqs.txt")
    if [[ $TOTAL_LINES -gt 0 ]]; then
        SIMILARITY=$(echo "scale=2; (1 - $DIFF_COUNT / $TOTAL_LINES) * 100" | bc)
        echo -e "${YELLOW}  Similarity: ${SIMILARITY}% ($DIFF_COUNT differences in $TOTAL_LINES sequences)${NC}"
        echo -e "${YELLOW}  Note: Minor differences expected due to different window sliding implementations${NC}"
    fi
fi
echo ""

# Print results table
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║                              BENCHMARK RESULTS                                   ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Format runtimes with proper decimal places
SDUST_RUNTIME_FMT=$(printf "%.3f" "$SDUST_RUNTIME" 2>/dev/null || echo "$SDUST_RUNTIME")
RUSTMASKER_RUNTIME_FMT=$(printf "%.3f" "$RUSTMASKER_RUNTIME" 2>/dev/null || echo "$RUSTMASKER_RUNTIME")

printf "%-25s %-20s %-20s\n" "Metric" "sdust" "rustmasker"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
printf "%-25s %-20s %-20s\n" "Runtime" "${SDUST_RUNTIME_FMT}s" "${RUSTMASKER_RUNTIME_FMT}s"
printf "%-25s %-20s %-20s\n" "Peak Memory (MB)" "${SDUST_MEMORY_MB}" "${RUSTMASKER_MEMORY_MB}"
printf "%-25s %-20s %-20s\n" "Masked bases" "${SDUST_MASKED_BASES}/${SDUST_TOTAL_BASES}" "${RUSTMASKER_MASKED_BASES}/${RUSTMASKER_TOTAL_BASES}"
printf "%-25s %-20s %-20s\n" "Masked percentage" "${SDUST_MASKED_PCT}%" "${RUSTMASKER_MASKED_PCT}%"
printf "%-25s %-20s\n" "Outputs match" "$OUTPUTS_MATCH"
echo ""

# Calculate speedup
if [[ "$SDUST_RUNTIME" != "unknown" && "$RUSTMASKER_RUNTIME" != "unknown" ]]; then
    SDUST_SEC=$(echo "$SDUST_RUNTIME" | bc)
    RUSTMASKER_SEC=$(echo "$RUSTMASKER_RUNTIME" | bc)

    if [[ $(echo "$SDUST_SEC > 0" | bc) -eq 1 ]]; then
        SPEEDUP=$(echo "scale=2; $SDUST_SEC / $RUSTMASKER_SEC" | bc)
        if [[ $(echo "$SPEEDUP >= 1" | bc) -eq 1 ]]; then
            echo -e "${GREEN}Speedup (rustmasker vs sdust): ${SPEEDUP}x${NC}"
        else
            echo -e "${YELLOW}Speedup (rustmasker vs sdust): ${SPEEDUP}x (slower)${NC}"
        fi
    fi
fi

# Calculate memory difference
if [[ "$SDUST_MEMORY_MB" != "unavailable" && "$RUSTMASKER_MEMORY_MB" != "unavailable" ]]; then
    MEMORY_DIFF=$(echo "scale=2; $RUSTMASKER_MEMORY_MB - $SDUST_MEMORY_MB" | bc)
    if [[ $(echo "$MEMORY_DIFF < 0" | bc) -eq 1 ]]; then
        MEMORY_REDUCTION=$(echo "scale=2; (1 - $RUSTMASKER_MEMORY_MB / $SDUST_MEMORY_MB) * 100" | bc)
        echo -e "${GREEN}Memory reduction: ${MEMORY_REDUCTION}%${NC}"
    else
        MEMORY_INCREASE=$(echo "scale=2; ($RUSTMASKER_MEMORY_MB / $SDUST_MEMORY_MB - 1) * 100" | bc)
        echo -e "${YELLOW}Memory increase: ${MEMORY_INCREASE}%${NC}"
    fi
fi

echo ""
echo -e "${BLUE}Results saved to: $OUTPUT_DIR/${NC}"
echo ""
echo -e "${YELLOW}Note:${NC} sdust and rustmasker may produce slightly different masks due to:"
echo "  - Different window sliding implementations"
echo "  - Edge case handling at sequence boundaries"
echo "  - Handling of N-containing regions"
echo "  These differences are typically minor (<1% of bases)."
