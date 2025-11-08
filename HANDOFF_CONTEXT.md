# Handoff Context: mask_fastq Rust Utility Implementation

**Date:** 2025-10-07 (Updated: 2025-11-08)
**Branch:** claude/review-mike-issue323-011CUw8zmFc8fvpWUpyYjghu
**Issue:** #323 - MASK_FASTQ_READS process memory usage optimization

## Summary

Replaced bbmask.sh with a custom Rust utility (`mask_fastq`) to eliminate memory issues caused by bbmask.sh loading entire FASTQ files into Java heap memory. This was causing OOM errors with large ONT read files.

### Update 2025-11-08
- Binary built successfully (793KB)
- All 8 Rust unit tests passing (including canonical k-mer tests)
- Manual testing completed - masking works correctly for all low-complexity sequences
- **Issue discovered and RESOLVED:** GCGCGC alternating pattern was not masked initially
- **Root cause:** bbmask.sh uses canonical k-mers (lexicographically smaller of k-mer and reverse complement)
- **Fix applied:** Modified mask_fastq to use canonical k-mers in entropy calculation
- **Result:** ✅ All nf-tests passing (2/2)

## What Was Completed

### 1. Rust Utility Implementation
- **Location:** `modules/local/maskRead/mask_fastq/`
- **Binary:** `modules/local/maskRead/mask_fastq/target/release/mask_fastq`
- **Size:** 767KB

**Features:**
- Streaming FASTQ processing (stdin → gzipped stdout)
- Shannon entropy calculation for low-complexity masking
- Parameters: `-w` (window size), `-e` (entropy threshold), `-k` (kmer size)
- Masks low-entropy regions with 'N' bases and '#' quality scores
- All 5 unit tests passing

### 2. Nextflow Process Updated
- **File:** `modules/local/maskRead/main.nf:26`
- **Changes:**
  - Removed `label "BBTools"`
  - Replaced bbmask.sh command with mask_fastq utility
  - Command: `zcat -f !{reads} | ${mask_fastq} -w !{window_size} -e !{entropy} -k 5 > ${out}`
  - Empty file handling preserved (checked before calling mask_fastq)

### 3. Manual Testing Completed
Tested mask_fastq with `test-data/toy-data/test-random-low-complexity.fastq`:
- ✅ Low complexity sequences (all A's, all T's) correctly masked to N's
- ✅ High complexity sequences preserved unchanged
- ✅ Quality scores masked to '#' for low-complexity regions
- ✅ FASTQ format maintained
- ✅ Gzip output works correctly

## What Needs to Be Tested

### Required: nf-test Suite
Run the existing nf-test suite to verify integration:

```bash
nf-test test tests/modules/local/maskRead/main.nf.test
```

**Expected Results:**

#### Test 1: "Should run correctly on FASTQ data"
- Process should succeed
- Input and output should have same read IDs (16 reads)
- **First 5 sequences should contain 'N' (masked)**
- **Remaining 11 sequences should NOT be masked** (identical to input)

**✅ TEST ISSUE RESOLVED (2025-11-08):**

Initial issue: **sequence 3** (GCGCGC repeating pattern) was NOT being masked.

**Analysis:**
- Sequences 1, 2, 4, 5: All A's or T's → Entropy ≈ 0.0 → **MASKED** ✓
- Sequence 3: GCGCGC alternating → Two 5-mers (GCGCG, CGCGC) → Entropy ≈ 1.0 → **NOT MASKED** ✗

**Investigation:**
Compared with bbmask.sh, which DOES mask GCGCGC. Discovered that bbmask uses **canonical k-mers**.

**Solution:**
- GCGCG and CGCGC are reverse complements
- Canonical form: both map to CGCGC (lexicographically smaller)
- Result: Only 1 unique k-mer → Entropy = 0.0 → **MASKED** ✓

**Implementation:**
Added canonical k-mer support to mask_fastq:
- `reverse_complement()` function
- `canonical_kmer()` function
- Modified `get_kmers()` to use canonical forms
- Added comprehensive tests

**Verification:**
- ✅ cargo test: 8/8 tests pass
- ✅ Manual test: All 5 sequences masked correctly
- ✅ nf-test: 2/2 tests pass

#### Test 2: "Should run correctly on empty FASTQ data"
- Process should succeed
- Output should be empty (0 reads)

### Full Pipeline Testing
Consider running a minimal pipeline test with real data:
```bash
nextflow run main.nf -profile test
```

## Technical Details

### Algorithm Comparison
Both bbmask.sh and mask_fastq use Shannon entropy:
- **Window size:** 25 bases (default parameter)
- **Entropy threshold:** 0.55 (default parameter)
- **K-mer size:** 5 (hardcoded in both)
- **Calculation:** Same Shannon entropy formula: `-Σ(p_i * log2(p_i))` over k-mer frequencies

### Key Differences
| Aspect | bbmask.sh | mask_fastq |
|--------|-----------|------------|
| Memory | Loads entire file into Java heap | Streams records one-by-one |
| Language | Java (BBTools) | Rust |
| Binary size | ~300MB (BBTools suite) | 767KB |
| Dependency | BBTools installation | Single static binary |

### File Locations
```
modules/local/maskRead/
├── main.nf                           # Updated Nextflow process
└── mask_fastq/
    ├── Cargo.toml                    # Rust project config
    ├── src/
    │   └── main.rs                   # Implementation (5/5 tests passing)
    └── target/release/
        └── mask_fastq                # Compiled binary (767KB)
```

## Potential Issues to Watch

### 1. Empty File Handling
The mask_fastq utility returns an error on empty input. This is handled by the Nextflow process which checks for empty input first:
```bash
if [[ -z "$(zcat "!{reads}" | head)" ]]; then
    echo -n | gzip > ${out}
else
    # Call mask_fastq
fi
```

### 2. Binary Path Resolution
The process uses `!{projectDir}` to reference the binary:
```bash
mask_fastq=!{projectDir}/modules/local/maskRead/mask_fastq/target/release/mask_fastq
```

Verify this resolves correctly in your Nextflow environment.

### 3. Cargo/Rust Availability
If the binary needs to be rebuilt:
```bash
cd modules/local/maskRead/mask_fastq
cargo build --release
```

Note: The binary has been compiled and is available at the path above (793KB). The `target/` directory is gitignored, so you'll need to rebuild if checking out on a new machine.

### 4. nf-test Installation
If nf-test is not available, install it:
```bash
# See: https://code.askimed.com/nf-test/docs/getting-started/
curl -fsSL https://code.askimed.com/install/nf-test | bash
```

## Next Steps

### Immediate
1. **Run nf-test suite** - Verify both test cases pass
2. **Review test output** - Ensure masking behavior matches expected results
3. **If tests fail:**
   - Compare masked output with bbmask output on same input
   - Check if entropy calculation differs
   - Verify k-mer extraction is correct

### Follow-up
1. **Performance testing** - Run on actual large ONT data to verify memory savings
2. **Update CHANGELOG.md** - Document the change
3. **Consider documentation** - Update any user-facing docs mentioning BBTools
4. **Merge to stable branch** - After tests pass and review complete

## Test Data Reference
- **Test file:** `test-data/toy-data/test-random-low-complexity.fastq`
- **Total reads:** 16
- **Low complexity reads:** First 5 (should be masked)
- **High complexity reads:** Remaining 11 (should not be masked)

## Questions or Issues?

If tests fail or behavior differs from expected:
1. Check the entropy calculation in `mask_fastq/src/main.rs:22-41`
2. Verify k-mer extraction logic in `mask_fastq/src/main.rs:44-51`
3. Compare with bbmask.sh output on a small test file
4. Review the unit tests which all pass: `cargo test` in mask_fastq directory

## Git Status
```
Current branch: mike-issue323-test
Main branch: master

Untracked files:
  modules/local/maskRead/mask_fastq/

Modified files:
  modules/local/maskRead/main.nf
```

**Important:** The mask_fastq directory is currently untracked. You may want to add it to the repo or add the binary to .gitignore depending on project conventions.
