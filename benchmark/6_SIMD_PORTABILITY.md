# SIMD Portability Concerns Explained

## What is SIMD and Why is it Architecture-Specific?

**SIMD** (Single Instruction, Multiple Data) refers to special CPU instructions that process multiple data elements simultaneously. However, **different CPU architectures have completely different SIMD instruction sets**, which creates portability challenges.

### The Problem: No Universal SIMD

Unlike regular operations (add, multiply, branches) that work the same across all CPUs, SIMD instructions are **vendor and architecture specific**:

``` bg-bg-200
// Regular code - works everywhere
let result = a + b;  // ✓ Portable

// SIMD code - only works on specific CPUs
use std::arch::x86_64::*;
let result = _mm_add_epi32(a, b);  // ✗ Only x86_64 CPUs!
```

## Major CPU Architecture SIMD Instruction Sets

### 1. **x86_64 / Intel / AMD** (Most servers and desktops)

Multiple generations of SIMD, from oldest to newest:

\| Instruction Set \| Year \| Register Size \| Example \|
\|----------------\|------\|---------------\|---------\|
\| **SSE** (Streaming SIMD Extensions) \| 1999 \| 128-bit (16 bytes) \| `_mm_add_ps()` \|
\| **SSE2** \| 2001 \| 128-bit \| `_mm_add_epi32()` \|
\| **SSE3, SSE4** \| 2004-2006 \| 128-bit \| More operations \|
\| **AVX** \| 2011 \| 256-bit (32 bytes) \| `_mm256_add_ps()` \|
\| **AVX2** \| 2013 \| 256-bit \| `_mm256_add_epi32()` \|
\| **AVX-512** \| 2016 \| 512-bit (64 bytes) \| `_mm512_add_epi32()` \|

**Key point:** Even within x86_64, older CPUs don't support newer instructions!

### 2. **ARM** (Mobile devices, Apple Silicon, AWS Graviton)

\| Instruction Set \| Register Size \| Used In \|
\|----------------\|---------------\|---------\|
\| **NEON** \| 128-bit \| All modern ARM CPUs \|
\| **SVE** (Scalable Vector Extension) \| 128-2048 bit \| High-end ARM servers \|

**Different functions than x86!** `vaddq_s32()` instead of `_mm_add_epi32()`

### 3. **Other Architectures**

- **PowerPC**: AltiVec
- **RISC-V**: RVV (RISC-V Vector Extension)
- **WebAssembly**: SIMD128

**Each has completely different instruction names and behavior!**

## Portability Problems for mask_fastq

### Problem 1: Code Won't Compile

If we write SIMD code for x86_64:

```rust
use std::arch::x86_64::*;

pub fn encode_kmer_simd(bases: &[u8]) -> Option<u16> {
    unsafe {
        let vec = _mm_loadu_si128(bases.as_ptr() as *const __m128i);
        // ... SIMD operations ...
    }
}
```

**What happens on other platforms:**

```sh
# Compiling on ARM (Mac M1/M2, Raspberry Pi, AWS Graviton)
$ cargo build
error: failed to resolve: use of undeclared crate or module `x86_64`
  --> src/lib.rs:1:17
   |
1  | use std::arch::x86_64::*;
   |                 ^^^^^^ use of undeclared crate or module `x86_64`

# Compiling on 32-bit systems, RISC-V, etc.
# Same error - module doesn't exist!
```

### Problem 2: Runtime Crashes

Even on x86_64, old CPUs crash on new instructions:

```rust
// Using AVX2 (from 2013)
unsafe {
    let result = _mm256_add_epi32(a, b);  // Requires AVX2
}
```

**Running on a 2010 CPU (pre-AVX2):**

```sh
Illegal instruction (core dumped)
```

No compiler error, just **instant crash at runtime**!

### Problem 3: Performance Cliffs

Some instructions only perform well on specific CPU generations:

```rust
// This instruction exists in AVX but is SLOW on some CPUs
let result = _mm256_permutevar8x32_epi32(a, b);
```

- AMD Zen1 (2017): ~6 cycles
- AMD Zen2 (2019): ~3 cycles
- Intel Skylake: ~3 cycles
- Intel Haswell: ~10+ cycles

**Same code, 3x performance difference!**

## Solutions and Workarounds

### Solution 1: Conditional Compilation (What We'd Need to Do)

Use Rust's `cfg` attributes to compile different code for different platforms:

```rust
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

#[cfg(target_arch = "aarch64")]  // ARM 64-bit
use std::arch::aarch64::*;

// x86_64 version
#[cfg(target_arch = "x86_64")]
pub fn encode_kmer_simd(bases: &[u8]) -> Option<u16> {
    unsafe {
        // Use SSE/AVX instructions
        let vec = _mm_loadu_si128(bases.as_ptr() as *const __m128i);
        // ... x86 SIMD code ...
    }
}

// ARM version
#[cfg(target_arch = "aarch64")]
pub fn encode_kmer_simd(bases: &[u8]) -> Option<u16> {
    unsafe {
        // Use NEON instructions
        let vec = vld1q_u8(bases.as_ptr());
        // ... ARM SIMD code ...
    }
}

// Fallback for everything else
#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
pub fn encode_kmer_simd(bases: &[u8]) -> Option<u16> {
    // Call the non-SIMD version
    encode_kmer(bases)
}
```

**Maintenance burden:** Now you have to maintain 3+ versions of the same function!

### Solution 2: Runtime CPU Feature Detection

Check at runtime if CPU supports the instruction:

```rust
pub fn encode_kmer_fast(bases: &[u8]) -> Option<u16> {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            return unsafe { encode_kmer_avx2(bases) };
        }
        if is_x86_feature_detected!("sse2") {
            return unsafe { encode_kmer_sse2(bases) };
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        return unsafe { encode_kmer_neon(bases) };
    }

    // Fallback
    encode_kmer(bases)
}
```

**Now you need:** SSE2 version + AVX2 version + ARM version + fallback = **4 implementations**!

### Solution 3: Use Portable SIMD Library

Rust has experimental portable SIMD (`std::simd`) that abstracts differences:

```rust
#![feature(portable_simd)]  // Requires nightly Rust
use std::simd::*;

pub fn encode_kmer_simd(bases: &[u8]) -> Option<u16> {
    let vec = u8x16::from_slice(bases);
    // Same code works on x86, ARM, etc!
}
```

**Problems:**

- Requires **nightly Rust** (unstable)
- Not all SIMD operations are supported
- Still generates architecture-specific code at compile time

### Solution 4: External Crate (What Many Projects Do)

Use battle-tested libraries:

```
[dependencies]
# Automatically handles x86/ARM/fallback
packed_simd = "0.3"
# or
simdeez = "1.0"
```

**Trade-off:** Add dependency, but someone else handles portability

## Real-World Example: Our mask_fastq Code

### Current Situation (No SIMD)

```rust
pub fn encode_kmer(bases: &[u8]) -> Option<u16> {
    let mut encoded: u16 = 0;
    for &base in bases {
        let bits = match base {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => return None,
        };
        encoded = (encoded << 2) | bits;
    }
    Some(encoded)
}
```

**Works on:** ✓ x86_64, ✓ ARM, ✓ RISC-V, ✓ PowerPC, ✓ WebAssembly, ✓ everywhere!

### With SIMD (Minimal Viable Implementation)

```rust
pub fn encode_kmer(bases: &[u8]) -> Option<u16> {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("sse2") {
            return unsafe { encode_kmer_sse2(bases) };
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        return unsafe { encode_kmer_neon(bases) };
    }

    // Fallback to portable version
    encode_kmer_portable(bases)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
unsafe fn encode_kmer_sse2(bases: &[u8]) -> Option<u16> {
    use std::arch::x86_64::*;
    // ... 50 lines of SIMD intrinsics ...
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn encode_kmer_neon(bases: &[u8]) -> Option<u16> {
    use std::arch::aarch64::*;
    // ... 50 lines of DIFFERENT SIMD intrinsics ...
}

fn encode_kmer_portable(bases: &[u8]) -> Option<u16> {
    // Original implementation
}
```

**Code growth:** 20 lines → **150+ lines**

## Testing Challenges

### Testing on Different Architectures


You need to test on **every platform you support**:

```sh
# Test x86_64 with SSE2 (most CPUs)
cargo test --target x86_64-unknown-linux-gnu

# Test x86_64 with AVX2 (newer CPUs)
cargo test --target x86_64-unknown-linux-gnu
# But how do you force AVX2 vs SSE2 in tests?

# Test ARM
cargo test --target aarch64-unknown-linux-gnu

# Test RISC-V
cargo test --target riscv64gc-unknown-linux-gnu

# Test WebAssembly
cargo test --target wasm32-wasi
```

**Problem:** Most developers only have x86_64! How do you test ARM/RISC-V?

### Cross-Platform CI Needed

```yaml
# .github/workflows/test.yml
jobs:
  test:
    strategy:
      matrix:
        arch: [x86_64, aarch64, riscv64]
        os: [ubuntu, macos, windows]
        # Plus testing different CPU features!
```

**Cost:** Longer CI times, more complex setup

## Documentation Burden

Users need to know about architecture differences:

```
## Performance

Performance varies by CPU architecture:

- x86_64 with AVX2: ~6,000 reads/sec
- x86_64 with SSE2: ~5,000 reads/sec
- ARM with NEON: ~5,500 reads/sec
- Other platforms: ~2,000 reads/sec (no SIMD)

To check your CPU features:
- Linux: `cat /proc/cpuinfo | grep flags`
- macOS: `sysctl machdep.cpu.features`
```

## Distribution Complexity

### Binary Distribution

If you distribute pre-built binaries:

```
mask_fastq-linux-x86_64-avx2
mask_fastq-linux-x86_64-sse2
mask_fastq-linux-aarch64
mask_fastq-linux-riscv64
mask_fastq-macos-x86_64
mask_fastq-macos-arm64
mask_fastq-windows-x86_64.exe
```

**7 binaries instead of 1!** (And that's simplified)

### Package Managers

Bioconda/Conda needs to know:

```yaml
# meta.yaml
requirements:
  build:
    - {{ compiler('c') }}
  host:
    - rust >=1.70
  run:
    # How to specify AVX2 requirement?
```

Some package managers can't express "use this version on AVX2 CPUs"!

## Real-World Portability Examples

### Example 1: Bioinformatics Tool (Minimap2)

Minimap2 is a widely-used sequence aligner with SIMD:

```
#ifdef __SSE2__
    // SSE2 code for x86
#elif defined(__ARM_NEON)
    // NEON code for ARM
#else
    // Portable fallback
#endif
```

**Result:**

- Works everywhere
- But maintainer spends time on SIMD variants
- Bug reports: "crashes on my old server" (AVX on non-AVX CPU)

### Example 2: Our Situation

Currently:

```
mask_fastq works on: x86_64, ARM, RISC-V, PowerPC, etc.
Installation: cargo install mask_fastq  ← Works everywhere!
```

With SIMD:

```
mask_fastq works best on: x86_64 with SSE2+
Also works on: ARM (slower SIMD), others (fallback)
Installation: cargo install mask_fastq  ← Still works, but...
  - x86_64 users: Great performance
  - ARM users: Need to verify NEON implementation
  - Other users: Report issues "slower than advertised"
```

## Recommendation for mask_fastq

### Option A: x86_64 Only + Fallback (Pragmatic)

```rust
pub fn encode_kmer(bases: &[u8]) -> Option<u16> {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("sse2") {  // Present on all 64-bit x86 CPUs
            return unsafe { encode_kmer_sse2(bases) };
        }
    }

    // Fallback for ARM, RISC-V, old x86, etc.
    encode_kmer_portable(bases)
}
```

**Pros:**

- Covers 95%+ of servers (x86_64)
- Only 2 implementations to maintain
- SSE2 is guaranteed on x86_64

**Cons:**

- ARM users get slower version
- Still need to test fallback path

### Option B: x86 + ARM (Comprehensive)

Support the two most common architectures:

**Pros:**

- Covers 99%+ of real-world usage
- Good performance everywhere

**Cons:**

- 3 implementations (x86, ARM, fallback)
- Need ARM test infrastructure

### Option C: Use Portable SIMD Crate

``` bg-bg-200
[dependencies]
packed_simd_2 = "0.3"
```

**Pros:**

- One implementation
- Library handles portability

**Cons:**

- External dependency
- May not be as optimized as hand-written

### Option D: No SIMD (Current State)

**Pros:**

- [x] Works everywhere
- [x] One simple implementation
- [x] Easy to test
- [x] Easy to maintain
- [x] Already fast (3.2x speedup from array optimization)

**Cons:**

- Missing potential 1.4x additional speedup

## My Recommendation

For `mask_fastq`, I'd suggest **Option A** if you decide to pursue SIMD:

1.  **Implement SSE2 version for x86_64**
  - SSE2 is guaranteed available (part of x86_64 spec)
  - No runtime detection needed
  - Covers majority of users
2.  **Keep portable fallback**
  - Maintains universal compatibility
  - Only ~10-20% slower than SIMD
  - Required for non-x86 platforms
3.  **Document architecture differences**
    - Be transparent about performance expectations
    - Provide benchmark numbers per platform

**This gives you:**

- 1.4x speedup for most users (x86_64)
- Universal compatibility maintained
- Only 2 implementations to test/maintain
- No complex runtime detection
- No external dependencies

**Code structure:**

```rust
// lib.rs
#[cfg(target_arch = "x86_64")]
mod simd_x86;

pub fn encode_kmer(bases: &[u8]) -> Option<u16> {
    #[cfg(target_arch = "x86_64")]
    {
        unsafe { simd_x86::encode_kmer_sse2(bases) }
    }

    #[cfg(not(target_arch = "x86_64"))]
    {
        encode_kmer_portable(bases)
    }
}
```

Only ~100 extra lines for SSE2, maintains portability, achieves most of the speedup!
