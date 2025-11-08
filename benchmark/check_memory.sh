#!/bin/bash
# Quick memory check script

echo "System Memory Status:"
echo "===================="

if command -v free &> /dev/null; then
    free -h
    echo ""

    TOTAL_MEM_KB=$(free | grep Mem | awk '{print $2}')
    AVAIL_MEM_KB=$(free | grep Mem | awk '{print $7}')
    USED_MEM_KB=$(free | grep Mem | awk '{print $3}')

    TOTAL_MEM_GB=$(echo "scale=1; $TOTAL_MEM_KB / 1024 / 1024" | bc)
    AVAIL_MEM_GB=$(echo "scale=1; $AVAIL_MEM_KB / 1024 / 1024" | bc)
    USED_MEM_GB=$(echo "scale=1; $USED_MEM_KB / 1024 / 1024" | bc)

    echo "Summary:"
    echo "  Total:     ${TOTAL_MEM_GB} GB"
    echo "  Used:      ${USED_MEM_GB} GB"
    echo "  Available: ${AVAIL_MEM_GB} GB"
    echo ""

    # Recommendations
    echo "Recommendations for MAX_BBMASK_MEMORY_GB:"

    if (( $(echo "$TOTAL_MEM_GB < 8" | bc -l) )); then
        echo "  ⚠️  Small instance (< 8GB total)"
        echo "  → Recommended: export MAX_BBMASK_MEMORY_GB=2"
        echo "  → Skip ultra-long benchmarks"
    elif (( $(echo "$TOTAL_MEM_GB < 16" | bc -l) )); then
        echo "  ⚙️  Medium instance (8-16GB total)"
        echo "  → Recommended: export MAX_BBMASK_MEMORY_GB=4"
        echo "  → Ultra-long benchmarks should work"
    else
        echo "  ✅ Large instance (16GB+ total)"
        echo "  → Recommended: export MAX_BBMASK_MEMORY_GB=8"
        echo "  → All benchmarks should run smoothly"
    fi
else
    echo "Error: 'free' command not found"
    echo "Cannot determine memory status"
fi
