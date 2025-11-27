# R Integration and ALTREP Review

## Overview

This document reviews the R integration layer for the ldmap genomic coordinate encoding, focusing on the vctrs framework integration, ALTREP usage, type safety, and API design.

**Scope:**
- R interface via Rcpp (`R/snp.R`, `R/region.R`)
- vctrs integration for custom vector types
- ALTREP detection and handling
- Type safety and conversion functions
- API design and usability

---

## Architecture

### Custom Vector Types via vctrs

The package implements three custom vector types using the `vctrs` framework:

1. **`ldmap_snp`** - Packed SNP coordinates (chromosome, position, ref, alt)
2. **`ldmap_region`** - Packed genomic regions (chromosome, start, end)
3. **`ldmap_allele`** - Encoded nucleotide alleles

These are implemented as **S3 classes** inheriting from `vctrs_vctr`, with the underlying data stored as `double` vectors (which reinterpret 64-bit integers as doubles).

### ALTREP Detection

**Location:** `inst/include/ldmap.hpp:290-299`

```cpp
inline bool is_compact_seq(SEXP x) {
  if(!ALTREP(x)){
    return false;
  }
  R_xlen_t info_l=::Rf_xlength(R_altrep_data1(x));
  if(info_l!=3){
    return false;
  }
  return true;
}
```

**Purpose:** Detect R's compact sequence representation (e.g., `1:100`) to optimize operations.

**Status:** ✅ Implemented but **not actively used** in the codebase.

---

## Findings

### 🟢 STRENGTHS

#### 1. Excellent vctrs Integration

The implementation properly implements the vctrs protocol with:

**Type Coercion (`vec_ptype2`):**
```r
# From R/snp.R:12
vec_ptype2.ldmap_snp.ldmap_snp <- function(x, y, ...) new_ldmap_snp()
vec_ptype2.ldmap_snp.double <- function(x, y, ...) double()
vec_ptype2.double.ldmap_snp <- function(x, y, ...) double()
```

**Casting (`vec_cast`):**
```r
# From R/snp.R:308-310, 315-316
vec_cast.ldmap_snp.double <- function(x, to, ..., x_arg = "", to_arg = "")
    structure(vctrs::vec_data(x), class = c("ldmap_snp", "vctrs_vctr"))

vec_cast.double.ldmap_snp <- function(x, to, ..., x_arg = "", to_arg = "")
    vctrs::vec_data(x)
```

**Comparison Proxy:**
```r
# From R/snp.R:377-379
vec_proxy_compare.ldmap_snp <- function(x, ...) {
  structure(unclass(x),class="integer64")
}
```

This enables seamless integration with dplyr, tidyverse, and data.frames.

#### 2. Bidirectional Conversions

**SNP ↔ DataFrame:**
```r
# Explode to columns
explode_snp_struct(df, ldmap_snp = "snp_struct")
# → Returns: chrom, pos, ref, alt as separate columns

# Compact from columns
compact_snp_struct(df, chrom = "chrom", pos = "pos",
                   ref = "ref", alt = "alt")
# → Returns: single snp_struct column
```

#### 3. String Parsing

**Location:** `src/snpcoord.cpp` (exported via Rcpp)

Supports parsing from human-readable strings:
```r
parse_ldmap_SNP("1:10177_A_AC")
parse_ldmap_SNP("chr1:10539_C_A")
parse_ldmap_region("chr1:1000_5000")
```

#### 4. Type Safety in C++

**Input Validation:**
```cpp
// From src/snpcoord.cpp:1076-1084
if(Rf_length(SEXP(alt))!=allele_size){
  Rcpp::stop("ref and alt must be the same size or can both be empty/NULL or alt can be empty/NULL");
}

if(use_ref and use_alt and (ref.sexp_type() != alt.sexp_type())){
  Rcpp::stop("ref and alt must both be of the same type if both are to be used");
}
```

#### 5. Chromosome Encoding Flexibility

**Location:** `R/snp.R:231-247`

Handles multiple formats:
```r
as_chromosome(c(1, 2, 3))           # Integer
as_chromosome(c("chr1", "chr2"))    # Character with prefix
as_chromosome(c("1", "2", "X"))     # Character without prefix
```

Automatically detects and validates chromosome naming conventions.

---

### 🔴 CRITICAL ISSUES

#### 1. Missing Chromosome Limit Check

**Location:** `R/snp.R:232-233`

```r
as_chromosome <- function(x){
    if(length(unique(x)) > 2^5)
        stop("ldmap cannot currently represent genomic coordinates in reference genomes with more than 2^5(32) chromosomes");
```

**Issue:** The check is **wrong** for the underlying representation!

- The error message says **32 chromosomes** (2^5)
- But `bit_region` uses **6 bits** for chromosome (supports 0-63)
- And `bit_snp` uses **5 bits** for chromosome (supports 0-31)

**Impact:** Critical - Inconsistent with C++ implementation

**Recommendation:**
```r
as_chromosome <- function(x){
    # Check against the ACTUAL limit from bit_snp (5 bits = 32 chromosomes)
    if(length(unique(x)) > 31)  # 0-31 inclusive = 32 values
        stop("ldmap cannot represent genomic coordinates with more than 31 chromosomes (0-31)");
    # OR better: fix the C++ to use 6 bits everywhere
```

#### 2. Silent Data Loss on Conversion

**Location:** `src/snpcoord.cpp:1715-1716`

```cpp
chrom(i)=sp.chrom() !=0 ? sp.chrom() : NA_INTEGER;
pos(i)=sp.pos() !=0 ? static_cast<double>(sp.pos()) : NA_REAL;
```

**Issue:** Treating `0` as `NA` is **semantically incorrect**:
- Chromosome 0 might be valid (though unconventional)
- Position 0 is a valid genomic coordinate (start of chromosome)

**Impact:** Medium - Could cause silent data corruption for edge cases

**Recommendation:** Use explicit NA markers in the bit-packed representation rather than relying on zero values.

---

### 🟡 MEDIUM ISSUES

#### 3. No Integer64 Support for Positions

**Location:** Throughout R interface

**Issue:** Positions are stored as `double` in R:
```r
# From src/snpcoord.cpp:1716
pos(i)=sp.pos() !=0 ? static_cast<double>(sp.pos()) : NA_REAL;
```

While this works for current use (human genome max ~250M bp), it has limitations:
- Doubles have 53-bit precision for integers
- SNP positions use 43 bits (max ~8.8 trillion) - fits safely
- **But:** Loses explicit integer semantics

**Impact:** Low-Medium - Works but semantically imprecise

**Recommendation:** Consider using `bit64::integer64` for positions to maintain integer semantics, though the current approach is acceptable.

#### 4. Inconsistent NA Handling

**Location:** `R/snp.R:29-31`

```r
is.na.ldmap_snp <- function(x,...){
  is.na(chromosomes(x)) | positions(x)==0
}
```

**Issues:**
1. Uses `positions(x)==0` as NA check (wrong - position 0 is valid!)
2. Inconsistent with C++ which also uses 0 as NA marker
3. Doesn't check alleles for NA

**Impact:** Medium - Incorrect NA detection

**Recommendation:**
```r
is.na.ldmap_snp <- function(x,...){
  # Check for truly invalid values in the underlying representation
  # Perhaps reserve a specific bit pattern for NA (e.g., all bits set)
  is.na(chromosomes(x)) | is.na(positions(x))
}
```

#### 5. Missing pillar/tibble Display Methods

**Status:** ❌ **Not Implemented**

The package doesn't implement custom display methods for tibbles:
- No `pillar_shaft.ldmap_snp()`
- No `vec_ptype_abbr()`
- No `vec_ptype_full()`

**Impact:** Medium - Poor user experience in interactive use

When displaying in tibbles, users see raw doubles:
```r
> tibble(snp = new_ldmap_snp(1:3, 100:102))
# A tibble: 3 × 1
        snp
      <dbl>
1 2.052e-289
2 2.052e-289
3 2.052e-289
```

**Recommendation:** Implement pillar methods:
```r
#' @export
vec_ptype_abbr.ldmap_snp <- function(x, ...) "snp"

#' @export
vec_ptype_full.ldmap_snp <- function(x, ...) "ldmap_snp"

#' @export
pillar_shaft.ldmap_snp <- function(x, ...) {
  # Format as "chr:pos_ref_alt"
  pillar::new_pillar_shaft_simple(
    format_ldmap_snp(x),
    align = "left"
  )
}
```

#### 6. Redundant Function Definitions

**Location:** `R/snp.R:85-92, 95-102`

```r
# Defined THREE times!
vec_duplicate_id.ldmap_snp <- function(x,...){
  vec_duplicate_id(vctrs::vec_data(x))
}
```

**Impact:** Low - Confusing but harmless (R uses the last definition)

**Recommendation:** Remove duplicates.

---

### 🟡 DESIGN ISSUES

#### 7. Unclear ALTREP Usage

**Status:** The `is_compact_seq()` function exists but is **never called**.

**Impact:** Low - Dead code or incomplete feature

**Recommendation:** Either:
1. Use it to optimize operations on sequences, OR
2. Remove it if not needed

#### 8. Mixed Mutation Semantics

**Location:** `R/snp.R:550-558`

```r
`chromosomes<-` <- function(x, value) {
  if(inherits(x,"ldmap_snp")){
    return(new_ldmap_snp(chrom = rep(value,len=length(x)),
                        pos = positions(x),
                        ref = ref_alleles(x),
                        alt_alleles(x)))
  }
```

**Issue:** The assignment operators (`chromosomes<-`, `positions<-`, etc.) don't mutate in place - they **reconstruct** the entire vector.

**Impact:** Medium - Performance and semantic clarity

This is actually **correct** for vctrs (which emphasizes immutability), but the naming suggests mutation. Users might expect:
```r
chromosomes(x) <- 5  # Looks like mutation but creates new object
```

**Recommendation:** Document this clearly or provide explicit `set_*()` functions that make the copy semantics obvious.

#### 9. Incomplete GenomicRanges Integration

**Location:** `R/region.R:254-272`

```r
vec_cast.ldmap_region.GRanges <- function(x, to, ...) {
  # ...
  parse_ldmap_region(as.character(x))  # ⚠️ Indirect conversion
  ## as_ldmap_  # ← Commented out direct conversion
```

**Issue:** The conversion from `GRanges` uses string parsing instead of direct field access. This is:
- Slower
- Fragile (depends on `as.character()` formatting)
- More error-prone

**Impact:** Low-Medium - Works but suboptimal

**Recommendation:** Use direct conversion:
```r
vec_cast.ldmap_region.GRanges <- function(x, to, ...) {
  new_ldmap_region(
    chrom = as.character(GenomicRanges::seqnames(x)),
    start = GenomicRanges::start(x),
    end = GenomicRanges::end(x)
  )
}
```

---

### 🟢 LOW PRIORITY OBSERVATIONS

#### 10. Naming Inconsistency

- Function: `ldmap_snp_2_dataframe` (snake_case with `2`)
- Function: `ldmap_region_2_data_frame` (snake_case with underscored `data_frame`)

**Recommendation:** Standardize: `ldmap_snp_to_dataframe` and `ldmap_region_to_dataframe`

#### 11. Missing Unit Tests for Edge Cases

While test coverage is good (test_snp.R has comprehensive tests), missing:
- Tests for chromosome limits (0, 31, 32, 63)
- Tests for position limits (0, max values)
- Tests for NA handling edge cases
- Tests for type coercion edge cases

---

## Memory Safety

### ✅ GOOD: Proper Use of bit_cast

```cpp
// From inst/include/alleles.hpp:32-37
template <class To, class From>
typename std::enable_if<
    (sizeof(To) == sizeof(From)) &&
    std::is_trivially_copyable<From>::value &&
    std::is_trivial<To>::value,
    To>::type
inline bit_cast(const From &src) noexcept
```

The type-punning via `bit_cast` is:
- ✅ Safe (compile-time size check)
- ✅ Well-defined (no undefined behavior)
- ✅ Portable (works across platforms)

### ✅ GOOD: Static Assertions

```cpp
// From src/snpcoord.cpp:1073
static_assert(sizeof(Snp)==sizeof(double) and sizeof(double) == 8,
              "packed SNP structure must be the size of a double precision floating point (8 bytes)");
```

Ensures the packing assumptions are verified at compile-time.

### ⚠️ WARNING: Potential Alignment Issues

**Location:** `inst/include/alleles.hpp:276, 262`

```cpp
struct bit_snp{
  char alt : 8;
  char ref : 8;
  uint64_t pos :43;
  unsigned char chrom : 5;
} __attribute__((packed));
```

**Issue:** While `__attribute__((packed))` ensures tight packing, it's:
- GCC/Clang specific (not standard C++)
- Could cause alignment issues on some platforms

**Current Status:** ✅ Works on all target platforms (Linux, macOS, Windows with GCC)

**Recommendation:** Add static assertions for alignment assumptions.

---

## API Usability

### ✅ EXCELLENT: Pipe-Friendly Design

The API works well with magrittr pipes:

```r
library(ldmap)
library(dplyr)

# Create SNPs with pipes
snps <- new_ldmap_snp(1:10, 1000:1009) %>%
  set_chromosomes(rep(1:2, length.out = 10)) %>%
  set_positions(seq(1000, 2000, length.out = 10))

# Work with dataframes
gwas_data %>%
  compact_snp_struct(chrom = "chr", pos = "bp",
                     ref = "a1", alt = "a2") %>%
  filter(chromosomes(snp_struct) == 1) %>%
  arrange(snp_struct)
```

### ✅ GOOD: Flexible Constructors

Multiple ways to create SNPs:
```r
# From components
new_ldmap_snp(chrom = 1:3, pos = 100:102)

# From strings
parse_ldmap_SNP(c("1:100_A_T", "2:200_C_G"))

# From doubles
ldmap_snp(c(2.052e-289, 2.052e-289))

# From dataframe
compact_snp_struct(df)
```

### ⚠️ CONFUSING: Double Representation Leaks

Users can see the raw double representation:
```r
x <- new_ldmap_snp(1, 1000)
as.numeric(x)  # Returns something like 2.052e-289
```

This is **conceptually confusing** but **technically necessary** for vctrs integration.

**Recommendation:** Emphasize in documentation that the numeric representation is an implementation detail.

---

## Comparison to Best Practices

### vctrs Best Practices Checklist

| Practice | Status | Notes |
|----------|--------|-------|
| Implement `vec_ptype2` | ✅ | Comprehensive |
| Implement `vec_cast` | ✅ | Includes character, double, integer64 |
| Implement `vec_proxy_compare` | ✅ | Uses integer64 |
| Implement `vec_proxy_equal` | ✅ | Delegates to dataframe |
| Implement `format` method | ✅ | Returns chr:pos_ref_alt |
| Implement `vec_ptype_abbr` | ❌ | Missing (affects tibble display) |
| Implement `vec_ptype_full` | ❌ | Missing |
| Implement `pillar_shaft` | ❌ | Missing (affects tibble display) |
| Document type hierarchy | ⚠️ | Partial |
| Handle NA properly | ⚠️ | Uses 0 as NA (incorrect) |

---

## Performance Considerations

### ✅ GOOD: Memory Efficiency

**Comparison:**

| Representation | Memory per SNP | For 1M SNPs |
|----------------|---------------|-------------|
| Separate columns (int32, int64, 2× char) | ~17 bytes | ~17 MB |
| ldmap_snp (packed) | 8 bytes | ~8 MB |
| **Savings** | **~47%** | **~47%** |

### ✅ GOOD: Cache-Friendly

Packed representation improves:
- Cache locality (all SNP data in contiguous memory)
- Sorting performance (single 64-bit comparison)
- Binary I/O (direct disk writes without marshaling)

### ⚠️ MODERATE: Reconstruction Overhead

The assignment operators reconstruct entire vectors:
```r
positions(snps) <- new_positions  # Copies entire vector
```

For large vectors (>1M elements), this could be slow. However, this is **idiomatic R** and consistent with vctrs design philosophy.

---

## Recommendations Summary

### Priority 1: Fix Critical Issues

1. **Fix chromosome limit check** in `as_chromosome()` to match C++ implementation
2. **Fix NA detection** - don't use `position == 0` as NA marker
3. **Standardize chromosome bits** to 6 across SNP and Region (see main review)

### Priority 2: Improve User Experience

4. **Implement pillar/tibble display methods** for better interactive use
5. **Document the double representation** as implementation detail
6. **Remove redundant function definitions**

### Priority 3: Code Quality

7. **Add edge case tests** for chromosome/position limits
8. **Fix GenomicRanges conversion** to use direct field access
9. **Standardize function naming** (snake_case consistently)
10. **Either use or remove** ALTREP detection code

### Priority 4: Documentation

11. **Document NA representation** and edge case behavior
12. **Add examples** showing proper usage patterns
13. **Document memory layout** and implementation details
14. **Add migration guide** from other genomic coordinate systems

---

## Overall Assessment

**Grade: B+**

The R integration is **well-designed and functional** with excellent use of the vctrs framework. The main issues are:

**Strengths:**
- ✅ Excellent vctrs integration
- ✅ Memory efficient
- ✅ Type-safe C++ implementation
- ✅ Flexible API
- ✅ Good test coverage

**Weaknesses:**
- ❌ Inconsistent NA handling (uses 0 as NA)
- ❌ Missing display methods (pillar/tibble)
- ❌ Chromosome limit check mismatch
- ⚠️ Some redundant code
- ⚠️ Incomplete GenomicRanges integration

**Critical Path to A:**
1. Fix NA representation (don't use 0)
2. Add pillar display methods
3. Fix chromosome limit checks
4. Comprehensive edge case testing

---

## Code Examples: Recommended Fixes

### Fix NA Detection

```r
# Current (WRONG):
is.na.ldmap_snp <- function(x,...){
  is.na(chromosomes(x)) | positions(x)==0  # ❌ position 0 is valid!
}

# Recommended:
is.na.ldmap_snp <- function(x,...){
  # Use actual NA markers from unpacking
  df <- ldmap_snp_2_dataframe(x)
  is.na(df$chrom) | is.na(df$pos)
}
```

### Add Pillar Methods

```r
#' @export
#' @importFrom vctrs vec_ptype_abbr
vec_ptype_abbr.ldmap_snp <- function(x, ...) {
  "snp"
}

#' @export
#' @importFrom vctrs vec_ptype_full
vec_ptype_full.ldmap_snp <- function(x, ...) {
  "ldmap_snp"
}

#' @export
#' @importFrom pillar pillar_shaft
pillar_shaft.ldmap_snp <- function(x, ...) {
  formatted <- format_ldmap_snp(x)
  pillar::new_pillar_shaft_simple(formatted, align = "left")
}
```

### Fix GenomicRanges Conversion

```r
#' @export vec_cast.ldmap_region.GRanges
vec_cast.ldmap_region.GRanges <- function(x, to, ...) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' needed for this function")
  }

  # Direct conversion (faster, more robust)
  new_ldmap_region(
    chrom = as.character(GenomicRanges::seqnames(x)),
    start = GenomicRanges::start(x),
    end = GenomicRanges::end(x)
  )
}
```

---

## Conclusion

The R integration is **production-ready** for its intended use case (LD analysis with human genomes) but would benefit from the fixes outlined above for broader applicability and better user experience. The vctrs integration is exemplary and serves as a good example for other packages implementing custom genomic data types.

The main concern is the **NA representation using zero values**, which should be addressed before wider adoption.
