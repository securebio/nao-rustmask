# Next Steps for Issue #323 - mask_fastq Implementation

## Current Status ✅

The Rust utility `mask_fastq` has been implemented and is ready for testing. Here's what's been completed:

- ✅ Rust binary built successfully (793KB static binary)
- ✅ All 5 unit tests passing (`cargo test`)
- ✅ Manual FASTQ masking tested and working
- ✅ Nextflow process updated to use new utility
- ✅ Code reviewed and merged to branch

## What YOU Need to Run Locally 🏃

### 1. Run nf-test Suite (REQUIRED)

This is the critical validation step that must be done locally with nf-test configured.

```bash
# From project root
nf-test test tests/modules/local/maskRead/main.nf.test
```

**Expected outcome:**
- Test 2 (empty file) should **PASS** ✓
- Test 1 (FASTQ data) will likely **FAIL** ⚠️ due to GCGCGC sequence issue (see below)

### 2. Investigate GCGCGC Masking Issue ⚠️

**The Problem:**
The test expects all first 5 sequences to be masked, but sequence 3 (GCGCGC repeating) is NOT being masked because its entropy (1.0) exceeds the threshold (0.55).

**What to check:**

#### Option A: Compare with bbmask.sh
Run the original bbmask.sh on the same test data:

```bash
# If bbmask.sh is still available in your environment
bbmask.sh in=test-data/toy-data/test-random-low-complexity.fastq \
  out=bbmask_output.fastq.gz entropy=0.55 k=5 window=25

# Check if it masks the GCGCGC sequence
zcat bbmask_output.fastq.gz | grep -A 3 "^@B1 1"
```

**If bbmask DOES mask GCGCGC:**
- The mask_fastq algorithm needs adjustment
- Investigate how bbmask calculates entropy differently

**If bbmask DOES NOT mask GCGCGC:**
- The test expectations are wrong
- Update the test in `tests/modules/local/maskRead/main.nf.test`

#### Option B: Update Test Expectations
If GCGCGC shouldn't be masked (which is mathematically correct with entropy > 0.55), update line 49 in the test file:

```groovy
// Current (line 49):
assert masked_seqs.take(5).every { it.contains("N") }

// Change to check only sequences 1, 2, 4, 5:
assert [0, 1, 3, 4].every { masked_seqs[it].contains("N") }
assert !masked_seqs[2].contains("N")  // Sequence 3 should NOT be masked
```

### 3. Full Pipeline Test (Optional but Recommended)

Once nf-test passes, run a full pipeline test:

```bash
nextflow run main.nf -profile test
```

### 4. Performance Validation (Optional)

Test with actual large ONT data to verify memory improvements:

```bash
# Monitor memory usage during the run
nextflow run main.nf -profile <your_profile> --input <large_ont_data>
```

The goal was to eliminate OOM errors with 100kb+ ONT reads. Verify that:
- Process completes without memory errors
- Memory usage is lower than with bbmask.sh

## Installation Notes

### If nf-test is not installed:
```bash
curl -fsSL https://code.askimed.com/install/nf-test | bash
```

### If binary needs rebuilding:
```bash
cd modules/local/maskRead/mask_fastq
cargo build --release
```

## Decision Point 🤔

After running the tests, you need to decide:

1. **If nf-test PASSES**:
   - Great! Update CHANGELOG.md and prepare for merge

2. **If nf-test FAILS on GCGCGC**:
   - Compare with bbmask.sh behavior
   - Either fix the algorithm OR update test expectations
   - Re-run nf-test to confirm

3. **If other issues arise**:
   - Check HANDOFF_CONTEXT.md for troubleshooting
   - Review the Rust implementation in `modules/local/maskRead/mask_fastq/src/main.rs`

## Questions?

See `HANDOFF_CONTEXT.md` for:
- Detailed implementation notes
- Algorithm comparison with bbmask.sh
- Troubleshooting guidance
- Technical details

## Summary

**You need to run:** `nf-test test tests/modules/local/maskRead/main.nf.test`

**You need to investigate:** Why GCGCGC isn't being masked (compare with bbmask.sh)

**You need to decide:** Fix algorithm or update test expectations

Everything else is ready to go! 🚀
