# Project Plan: Extract Common Permutation Utilities

## Executive Summary

**Project Goal**: Extract and consolidate duplicated integral permutation code across Psi4 modules into a reusable utility library with smart, concise APIs.

**Key Innovation**: Auto-inference of DPD indices and labels from input buffers
- Reduces API calls from 5 parameters to 1-2
- Eliminates redundant specification errors
- Automatic spin-case handling

**Impact**:
- Reduce ~200-300 lines of duplicated permutation logic
- Simplify API: `ovov_to_oovv(&I, file)` vs `buf4_sort(&I, file, prqs, idx1, idx2, label)`
- Improve maintainability and consistency across modules
- Establish foundation for future integral sorting refactoring

**Timeline**: 3-4 weeks (1 developer)
**Risk Level**: LOW
**Priority**: HIGH

---

## Table of Contents
1. [Background](#background)
2. [Current State Analysis](#current-state-analysis)
3. [Goals and Objectives](#goals-and-objectives)
4. [Scope](#scope)
5. [Design Proposal](#design-proposal)
6. [Implementation Plan](#implementation-plan)
7. [Testing Strategy](#testing-strategy)
8. [Timeline and Milestones](#timeline-and-milestones)
9. [Risk Assessment](#risk-assessment)
10. [Success Criteria](#success-criteria)
11. [Future Work](#future-work)

---

## Background

### Problem Statement

Integral sorting and transformation operations are duplicated across multiple Psi4 modules:
- **DCT**: 140+ permutation calls across RHF/UHF implementations
- **FNOCC**: Custom permutation logic in 2,437-line sortintegrals.cc
- **PSIMRCC**: Duplicate sorting in multiple files
- **cctransort, ccdensity, ccenergy**: Repetitive buf4_sort patterns

**Key Statistics:**
- 68 files use the `prqs` permutation (chemist → physicist notation)
- 10 buf4_sort calls in DCT RHF vs 57 in DCT UHF (nearly identical logic)
- ~200-300 lines of permutation logic could be consolidated

### Why This Matters

1. **Maintenance burden**: Changes to permutation logic must be replicated across modules
2. **Error-prone**: Inconsistent implementations lead to subtle bugs
3. **Readability**: Repetitive buf4_sort sequences obscure high-level intent
4. **Foundation**: Enables future refactoring of RHF/UHF duplication

---

## Current State Analysis

### Permutation Patterns Used

Based on analysis of `psi4/src/psi4/libdpd/buf4_sort.cc`, the following permutations are common:

| Pattern | Transformation | Meaning | Usage Count* |
|---------|---------------|---------|--------------|
| `prqs` | (pq\|rs) → (pr\|qs) | Chemist → Physicist notation | ~150+ |
| `rspq` | (pq\|rs) → (rs\|pq) | Swap bra/ket pairs | ~80+ |
| `psrq` | (pq\|rs) → (ps\|rq) | Mixed permutation | ~40+ |
| `qprs` | (pq\|rs) → (qp\|rs) | Swap within bra | ~30+ |
| `rqps` | (pq\|rs) → (rq\|ps) | Complex permutation | ~20+ |
| `pqsr` | (pq\|rs) → (pq\|sr) | Swap within ket | ~20+ |

*Approximate counts from grep analysis

### Common Integral Space Transformations

```cpp
// Pattern 1: Chemist to Physicist notation
// (OV|OV) → <OO|VV>
global_dpd_->buf4_sort(&I, file, prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");

// Pattern 2: Create transpose
// <OO|VV> → <VV|OO>
global_dpd_->buf4_sort(&I, file, rspq, ID("[V,V]"), ID("[O,O]"), "MO Ints <VV|OO>");

// Pattern 3: Create complementary orderings for contractions
// <OV|OV> → <VO|OV>
global_dpd_->buf4_sort(&I, file, qprs, ID("[V,O]"), ID("[O,V]"), "MO Ints <VO|OV>");
```

### Files with High Duplication

**DCT Module:**
- `dct_integrals_RHF.cc`: 10 permutation calls
- `dct_integrals_UHF.cc`: 57 permutation calls (same patterns × spin cases)
- `dct_intermediates_UHF.cc`: 20 permutation calls
- `dct_gradient_UHF.cc`: 5 permutation calls
- `dct_df_operations.cc`: 9 permutation calls

**CC Modules:**
- `cctransort/sort_tei_rhf.cc`: Repetitive A/B/C/D/E/F integral sorting
- `cctransort/sort_tei_uhf.cc`: Same patterns × spin cases
- `ccdensity/*.cc`: Scattered permutations throughout

**Analysis Findings:**
1. Nearly identical permutation sequences across modules
2. Only differences: file numbers, labels, spin labels (O/o, V/v)
3. No abstraction layer - direct DPD calls everywhere

---

## Goals and Objectives

### Primary Goals

1. **Reduce Code Duplication**
   - Consolidate ~200-300 lines of repetitive permutation logic
   - Eliminate copy-paste code across modules

2. **Improve Code Clarity**
   - Replace cryptic `buf4_sort` calls with semantic function names
   - Document permutation intent clearly

3. **Establish Foundation**
   - Create reusable library for current and future modules
   - Enable future RHF/UHF template consolidation

### Non-Goals (Out of Scope)

- ❌ Refactor FNOCC's custom out-of-core sorting algorithm
- ❌ Modify performance-critical inner loops
- ❌ Change DPD library internals
- ❌ Unify RHF/UHF implementations (future work)

---

## Scope

### In Scope

1. **Create new utility library** in `psi4/src/psi4/libtrans/`
2. **Wrap common permutation patterns** with semantic functions
3. **Refactor DCT module** to use new utilities (highest duplication)
4. **Refactor selected CC modules** (cctransort as proof-of-concept)
5. **Documentation** for permutation patterns and usage

### Out of Scope

- Complete refactoring of all 68 files (incremental approach)
- Performance optimization (maintain current performance)
- FNOCC module (highly specialized)
- Changing existing APIs (wrapper-only approach)

---

## Design Proposal

### Architecture

```
psi4/src/psi4/libtrans/
├── integral_permutations.h     (public interface)
├── integral_permutations.cc    (implementation)
└── integral_permutations_impl.h (template implementations)
```

### Public Interface

```cpp
// File: psi4/src/psi4/libtrans/integral_permutations.h

#ifndef _PSI_SRC_LIB_LIBTRANS_INTEGRAL_PERMUTATIONS_H_
#define _PSI_SRC_LIB_LIBTRANS_INTEGRAL_PERMUTATIONS_H_

#include "psi4/libdpd/dpd.h"
#include <string>

namespace psi {
namespace libtrans {

/**
 * @brief Utility class for common integral permutation operations
 *
 * This class provides semantic wrappers around DPD buf4_sort operations
 * to make integral transformations more readable and maintainable.
 *
 * **Design Philosophy:**
 * - Functions infer DPD indices and labels from input buffer when possible
 * - Optional parameters allow customization when needed
 * - Standard naming conventions are enforced for consistency
 *
 * All methods are static - no instantiation required.
 */
class IntegralPermutations {
public:
    //===========================================
    // Common Integral Space Transformations
    //===========================================

    /**
     * @brief Transform (OV|OV) → <OO|VV>
     *
     * Common transformation for building 2-electron integrals in
     * occupied-occupied/virtual-virtual space.
     *
     * **Smart Defaults:**
     * - Automatically determines output indices from input buffer
     * - Generates standard label "MO Ints <OO|VV>" (or <oo|vv>, <Oo|Vv> for spin cases)
     * - Handles all spin cases (RHF, UHF-AA, UHF-BB, UHF-AB) automatically
     *
     * @param InBuf Input (OV|OV) buffer (reads DPD parameters to infer output)
     * @param outfilenum Output file number
     * @param label Optional custom label (auto-generated if empty)
     *
     * @code
     * // Simple case - everything inferred
     * IntegralPermutations::ovov_to_oovv(&I, PSIF_LIBTRANS_DPD);
     *
     * // Custom label if needed
     * IntegralPermutations::ovov_to_oovv(&I, PSIF_LIBTRANS_DPD, "Custom <OO|VV>");
     * @endcode
     */
    static void ovov_to_oovv(
        dpdbuf4 *InBuf,
        int outfilenum,
        const std::string &label = ""
    );

    /**
     * @brief Transform (VV|OO) → <OV|OV>
     *
     * Creates <OV|OV> integrals from (VV|OO) chemist's notation.
     * Auto-detects spin case and generates appropriate indices/label.
     *
     * @param InBuf Input (VV|OO) buffer
     * @param outfilenum Output file
     * @param label Optional custom label (default: "MO Ints <OV|OV>")
     */
    static void vvoo_to_ovov(
        dpdbuf4 *InBuf,
        int outfilenum,
        const std::string &label = ""
    );

    /**
     * @brief Create transpose: <OO|VV> → <VV|OO>
     *
     * Swaps bra and ket pairs. Auto-detects indices from input.
     *
     * @param InBuf Input <OO|VV> buffer
     * @param outfilenum Output file
     * @param label Optional custom label (default: "MO Ints <VV|OO>")
     */
    static void transpose_oovv(
        dpdbuf4 *InBuf,
        int outfilenum,
        const std::string &label = ""
    );

    /**
     * @brief Generic transpose: <AB|CD> → <CD|AB>
     *
     * Swaps bra and ket pairs for any integral type.
     *
     * @param InBuf Input buffer
     * @param outfilenum Output file
     * @param label Optional custom label (auto-generated from input if empty)
     */
    static void transpose_bra_ket(
        dpdbuf4 *InBuf,
        int outfilenum,
        const std::string &label = ""
    );

    //===========================================
    // Additional Common Transformations
    //===========================================

    /**
     * @brief Transform (VO|OO) → <OV|OO>
     *
     * Common for OOOV-type integrals. Auto-detects spin case.
     *
     * @param InBuf Input (VO|OO) buffer
     * @param outfilenum Output file
     * @param label Optional custom label
     */
    static void vooo_to_ovoo(
        dpdbuf4 *InBuf,
        int outfilenum,
        const std::string &label = ""
    );

    /**
     * @brief Transform (OV|VV) → <OV|VV>
     *
     * Common for OVVV-type integrals. Auto-detects spin case.
     *
     * @param InBuf Input (OV|VV) buffer
     * @param outfilenum Output file
     * @param label Optional custom label
     */
    static void ovvv_to_ovvv(
        dpdbuf4 *InBuf,
        int outfilenum,
        const std::string &label = ""
    );

    //===========================================
    // Notation Conversions (General)
    //===========================================

    /**
     * @brief Convert chemist's notation (pq|rs) to physicist's notation <pr|qs>
     *
     * Applies the prqs permutation: (pq|rs) → (pr|qs)
     * This is the most common transformation in quantum chemistry codes.
     *
     * Unlike the specific functions above, this requires explicit indices
     * because the space types cannot be inferred from generic (pq|rs).
     *
     * @param InBuf Input buffer in chemist's notation
     * @param outfilenum PSI file number for output
     * @param pq_indices DPD indices for bra (pr)
     * @param rs_indices DPD indices for ket (qs)
     * @param label String label for output buffer
     */
    static void chemist_to_physicist(
        dpdbuf4 *InBuf,
        int outfilenum,
        int pq_indices,
        int rs_indices,
        const std::string &label
    );

    //===========================================
    // Batch Operations
    //===========================================

    /**
     * @brief Perform standard OVOV transformations in batch
     *
     * Common pattern: Generate all permutations needed for OVOV integrals:
     * - (OV|OV) → <OO|VV>
     * - <OO|VV> → <VV|OO> (if requested)
     *
     * All indices and labels are auto-generated.
     *
     * @param InBuf Input (OV|OV) buffer
     * @param outfilenum Output file
     * @param generate_transpose Also create <VV|OO>
     */
    static void standard_ovov_transformations(
        dpdbuf4 *InBuf,
        int outfilenum,
        bool generate_transpose = true
    );

    //===========================================
    // Low-Level Interface (for advanced use)
    //===========================================

    /**
     * @brief Apply arbitrary permutation with full control
     *
     * Direct wrapper around DPD buf4_sort for cases where automatic
     * inference isn't suitable. Most code should use the semantic
     * functions above instead.
     *
     * @param InBuf Input buffer
     * @param outfilenum Output file
     * @param permutation DPD indices enum (prqs, rspq, etc.)
     * @param pq_indices New bra indices
     * @param rs_indices New ket indices
     * @param label Output label
     */
    static void apply_permutation(
        dpdbuf4 *InBuf,
        int outfilenum,
        indices permutation,
        int pq_indices,
        int rs_indices,
        const std::string &label
    );

    //===========================================
    // Utilities
    //===========================================

    /**
     * @brief Get human-readable name for permutation
     *
     * @param permutation DPD indices enum
     * @return String describing the permutation (e.g., "Chemist to Physicist")
     */
    static std::string get_permutation_name(indices permutation);

    /**
     * @brief Validate permutation is supported
     *
     * @param permutation DPD indices enum
     * @return true if permutation is implemented and tested
     */
    static bool is_permutation_supported(indices permutation);

private:
    //===========================================
    // Internal Helpers
    //===========================================

    /**
     * @brief Derive output OO indices from input OV indices
     *
     * Analyzes input buffer's DPD parameters to determine appropriate
     * output indices for OO space, handling spin cases correctly.
     */
    static int derive_oo_indices(dpdbuf4 *InBuf);

    /**
     * @brief Derive output VV indices from input OV indices
     */
    static int derive_vv_indices(dpdbuf4 *InBuf);

    /**
     * @brief Generate standard label for integral type
     *
     * @param bra_space Space designation for bra (e.g., "OO", "oo", "Oo")
     * @param ket_space Space designation for ket (e.g., "VV", "vv", "Vv")
     * @return Standard label like "MO Ints <OO|VV>"
     */
    static std::string generate_standard_label(
        const std::string &bra_space,
        const std::string &ket_space
    );

    /**
     * @brief Detect spin case from DPD indices
     *
     * @return String: "AA", "BB", or "AB" for unrestricted; "R" for restricted
     */
    static std::string detect_spin_case(dpdbuf4 *InBuf);
};

} // namespace libtrans
} // namespace psi

#endif
```

### Implementation Example

```cpp
// File: psi4/src/psi4/libtrans/integral_permutations.cc

#include "integral_permutations.h"
#include "psi4/libdpd/dpd.h"
#include <stdexcept>
#include <sstream>

namespace psi {
namespace libtrans {

//===========================================
// Internal Helper Implementations
//===========================================

std::string IntegralPermutations::detect_spin_case(dpdbuf4 *InBuf) {
    // Examine DPD index labels to determine spin case
    // This is a simplified example - actual implementation would
    // parse the index structure more carefully

    int pqnum = InBuf->params->pqnum;
    int rsnum = InBuf->params->rsnum;

    // Check if indices contain mixed case (uppercase/lowercase)
    // indicating alpha-beta spin case
    // DPD convention: O=alpha-occ, o=beta-occ, V=alpha-vir, v=beta-vir

    // This is pseudocode - actual implementation would use
    // DPD's index lookup tables
    if (/* indices indicate mixed spin */) {
        return "AB";
    } else if (/* indices indicate beta-beta */) {
        return "BB";
    } else {
        return "AA";  // or "R" for restricted
    }
}

int IntegralPermutations::derive_oo_indices(dpdbuf4 *InBuf) {
    // Extract occupied-occupied index from occupied-virtual input
    // Using DPD's index combination tables

    std::string spin_case = detect_spin_case(InBuf);

    if (spin_case == "AA") {
        return ID("[O,O]");
    } else if (spin_case == "BB") {
        return ID("[o,o]");
    } else {  // AB
        return ID("[O,o]");
    }
}

int IntegralPermutations::derive_vv_indices(dpdbuf4 *InBuf) {
    // Extract virtual-virtual index from occupied-virtual input

    std::string spin_case = detect_spin_case(InBuf);

    if (spin_case == "AA") {
        return ID("[V,V]");
    } else if (spin_case == "BB") {
        return ID("[v,v]");
    } else {  // AB
        return ID("[V,v]");
    }
}

std::string IntegralPermutations::generate_standard_label(
    const std::string &bra_space,
    const std::string &ket_space)
{
    std::ostringstream oss;
    oss << "MO Ints <" << bra_space << "|" << ket_space << ">";
    return oss.str();
}

//===========================================
// Public API Implementations
//===========================================

void IntegralPermutations::ovov_to_oovv(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label)
{
    // Validate input
    if (InBuf == nullptr) {
        throw std::runtime_error("IntegralPermutations::ovov_to_oovv: null input buffer");
    }

    // Auto-detect output indices
    int oo_indices = derive_oo_indices(InBuf);
    int vv_indices = derive_vv_indices(InBuf);

    // Generate label if not provided
    std::string final_label = label;
    if (label.empty()) {
        std::string spin_case = detect_spin_case(InBuf);
        std::string bra, ket;

        if (spin_case == "AA") {
            bra = "OO"; ket = "VV";
        } else if (spin_case == "BB") {
            bra = "oo"; ket = "vv";
        } else {
            bra = "Oo"; ket = "Vv";
        }

        final_label = generate_standard_label(bra, ket);
    }

    // Perform transformation using DPD library
    // (OV|OV) → <OO|VV> using prqs permutation
    global_dpd_->buf4_sort(InBuf, outfilenum, prqs, oo_indices, vv_indices, final_label);
}

void IntegralPermutations::transpose_oovv(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label)
{
    // Swap OO and VV to get VV|OO
    // Input should be <OO|VV> format (physicist's notation)

    if (InBuf == nullptr) {
        throw std::runtime_error("IntegralPermutations::transpose_oovv: null input buffer");
    }

    // Read indices from input (they're already in the right format)
    int oo_indices = InBuf->params->pqnum;  // Currently OO
    int vv_indices = InBuf->params->rsnum;  // Currently VV

    // Generate label if not provided
    std::string final_label = label;
    if (label.empty()) {
        // Detect spaces from input and reverse them
        std::string spin_case = detect_spin_case(InBuf);
        // ... generate "MO Ints <VV|OO>" with appropriate capitalization
        final_label = "MO Ints <VV|OO>";  // Simplified
    }

    // Swap bra and ket using rspq permutation
    global_dpd_->buf4_sort(InBuf, outfilenum, rspq, vv_indices, oo_indices, final_label);
}

void IntegralPermutations::transpose_bra_ket(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label)
{
    // Generic transpose for any integral type
    if (InBuf == nullptr) {
        throw std::runtime_error("IntegralPermutations::transpose_bra_ket: null input buffer");
    }

    int pq_indices = InBuf->params->pqnum;
    int rs_indices = InBuf->params->rsnum;

    // Auto-generate label if not provided by reversing input label
    std::string final_label = label;
    if (label.empty()) {
        // Parse input label and reverse bra/ket
        // Simplified: just add "transposed" suffix
        final_label = std::string(InBuf->file.label) + " (transposed)";
    }

    global_dpd_->buf4_sort(InBuf, outfilenum, rspq, rs_indices, pq_indices, final_label);
}

void IntegralPermutations::standard_ovov_transformations(
    dpdbuf4 *InBuf,
    int outfilenum,
    bool generate_transpose)
{
    // Transform (OV|OV) → <OO|VV>
    ovov_to_oovv(InBuf, outfilenum);  // Uses auto-generated label

    if (generate_transpose) {
        // Open the just-created <OO|VV> buffer
        dpdbuf4 OOVVBuf;
        int oo_idx = derive_oo_indices(InBuf);
        int vv_idx = derive_vv_indices(InBuf);

        global_dpd_->buf4_init(&OOVVBuf, outfilenum, InBuf->file.my_irrep,
                              oo_idx, vv_idx, oo_idx, vv_idx, 0,
                              "MO Ints <OO|VV>");

        // Create <VV|OO> transpose
        transpose_oovv(&OOVVBuf, outfilenum);  // Auto-generates label

        global_dpd_->buf4_close(&OOVVBuf);
    }
}

void IntegralPermutations::chemist_to_physicist(
    dpdbuf4 *InBuf,
    int outfilenum,
    int pq_indices,
    int rs_indices,
    const std::string &label)
{
    // Low-level interface - requires explicit indices
    if (InBuf == nullptr) {
        throw std::runtime_error("IntegralPermutations::chemist_to_physicist: null input buffer");
    }

    // Perform transformation using DPD library
    // (pq|rs) → (pr|qs) using prqs permutation
    global_dpd_->buf4_sort(InBuf, outfilenum, prqs, pq_indices, rs_indices, label);
}

void IntegralPermutations::apply_permutation(
    dpdbuf4 *InBuf,
    int outfilenum,
    indices permutation,
    int pq_indices,
    int rs_indices,
    const std::string &label)
{
    // Direct wrapper for advanced use
    if (InBuf == nullptr) {
        throw std::runtime_error("IntegralPermutations::apply_permutation: null input buffer");
    }

    global_dpd_->buf4_sort(InBuf, outfilenum, permutation, pq_indices, rs_indices, label);
}

//===========================================
// Utility Functions
//===========================================

std::string IntegralPermutations::get_permutation_name(indices permutation) {
    switch(permutation) {
        case prqs: return "Chemist to Physicist notation (prqs)";
        case rspq: return "Transpose bra and ket (rspq)";
        case psrq: return "Mixed permutation (ps|rq)";
        case qprs: return "Swap within bra (qprs)";
        case pqsr: return "Swap within ket (pqsr)";
        case rqps: return "Complex permutation (rqps)";
        default: return "Unknown permutation";
    }
}

bool IntegralPermutations::is_permutation_supported(indices permutation) {
    // List of tested and supported permutations
    switch(permutation) {
        case prqs:
        case rspq:
        case psrq:
        case qprs:
        case pqsr:
        case rqps:
            return true;
        default:
            return false;
    }
}

} // namespace libtrans
} // namespace psi
```

### Usage Example (Before/After)

**Before (DCT RHF):**
```cpp
void DCTSolver::sort_OVOV_integrals_RHF() {
    dpdbuf4 I;
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,V]"),
                           "MO Ints <OO|VV>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, psrq, ID("[O,V]"), ID("[O,V]"),
                           "MO Ints SF <OV|OV>:<Ov|oV>");
    global_dpd_->buf4_close(&I);
}
```

**After (DCT RHF with utilities):**
```cpp
void DCTSolver::sort_OVOV_integrals_RHF() {
    dpdbuf4 I;
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");

    // Auto-infers indices and generates standard label!
    IntegralPermutations::ovov_to_oovv(&I, PSIF_LIBTRANS_DPD);

    // Custom label for special case
    IntegralPermutations::apply_permutation(&I, PSIF_LIBTRANS_DPD, psrq,
                                           ID("[O,V]"), ID("[O,V]"),
                                           "MO Ints SF <OV|OV>:<Ov|oV>");

    global_dpd_->buf4_close(&I);
}
```

**Even Better - Common Pattern:**
```cpp
void DCTSolver::sort_OVOV_integrals_RHF() {
    dpdbuf4 I;
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");

    // One function call does it all!
    IntegralPermutations::standard_ovov_transformations(&I, PSIF_LIBTRANS_DPD);

    global_dpd_->buf4_close(&I);
}
```

**UHF Example - Handles All Spin Cases Automatically:**
```cpp
void DCTSolver::sort_OVOV_integrals_UHF() {
    dpdbuf4 I;

    // Alpha-Alpha case
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    IntegralPermutations::ovov_to_oovv(&I, PSIF_LIBTRANS_DPD);  // Generates "MO Ints <OO|VV>"
    global_dpd_->buf4_close(&I);

    // Beta-Beta case
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                           ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
    IntegralPermutations::ovov_to_oovv(&I, PSIF_LIBTRANS_DPD);  // Generates "MO Ints <oo|vv>"
    global_dpd_->buf4_close(&I);

    // Alpha-Beta case
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                           ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
    IntegralPermutations::ovov_to_oovv(&I, PSIF_LIBTRANS_DPD);  // Generates "MO Ints <Oo|Vv>"
    global_dpd_->buf4_close(&I);
}
```

**Benefits:**
- **Concise**: 1 line instead of 5 parameters for common cases
- **Self-documenting**: `ovov_to_oovv` vs cryptic `prqs`
- **Automatic**: Infers indices and labels from input buffer
- **Spin-aware**: Handles RHF/UHF cases automatically
- **Consistent**: Standard naming conventions enforced
- **Flexible**: Can still override with custom labels when needed
- **Safe**: Centralized validation and error checking

---

## Design Rationale

### Why Auto-Inference?

The original design in this plan required explicit specification of output indices and labels:
```cpp
// Verbose version (rejected)
IntegralPermutations::ovov_to_oovv(&I, file, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");
```

This has several problems:

1. **Redundancy**: The function name `ovov_to_oovv` already encodes the transformation
   - Input: OV × OV → Output: OO × VV
   - We're specifying the same information three times!

2. **Error-Prone**: Easy to pass wrong indices
   ```cpp
   // Bug! Indices don't match function name
   IntegralPermutations::ovov_to_oovv(&I, file, ID("[V,V]"), ID("[O,O]"), "Wrong!");
   ```

3. **Maintenance Burden**: Every call site must be updated if conventions change

### The Solution: Smart Defaults

By reading the input buffer's DPD parameters, we can automatically infer:

1. **Spin Case**: Is this α-α, β-β, or α-β?
   - Detected from index structure ([O,V] vs [o,v] vs [O,v])

2. **Output Indices**: Derive from input indices
   - [O,V] input → [O,O] and [V,V] output
   - [o,v] input → [o,o] and [v,v] output
   - [O,v] input → [O,o] and [V,v] output

3. **Standard Labels**: Follow naming conventions
   - "MO Ints <OO|VV>" for α-α
   - "MO Ints <oo|vv>" for β-β
   - "MO Ints <Oo|Vv>" for α-β

### API Design Philosophy

**Principle**: Make the common case simple, keep advanced use possible

```cpp
// 95% of use cases: simple, clean, automatic
IntegralPermutations::ovov_to_oovv(&I, PSIF_LIBTRANS_DPD);

// Custom label when needed (4% of cases)
IntegralPermutations::ovov_to_oovv(&I, PSIF_LIBTRANS_DPD, "My Custom Label");

// Full control for edge cases (1% of cases)
IntegralPermutations::apply_permutation(&I, PSIF_LIBTRANS_DPD, prqs,
                                       custom_idx1, custom_idx2, "Custom");
```

### Comparison to Alternatives

| Approach | API Complexity | Flexibility | Error Risk |
|----------|----------------|-------------|------------|
| **Smart Defaults (Chosen)** | Low (1-2 params) | High | Low |
| Explicit Parameters | High (5 params) | High | High |
| Enum Spin Case | Medium (3 params) | Medium | Medium |
| Template Specialization | Low (1 param) | Medium | Low |

**Why not templates?**
```cpp
// Template approach would look like:
IntegralPermutations::ovov_to_oovv<AlphaAlpha>(&I, file);
```

Pros:
- Compile-time resolution
- Type-safe

Cons:
- Can't handle runtime spin case detection
- More complex implementation
- Harder to extend

**Decision**: Runtime auto-inference is more flexible and just as safe with validation.

### DPD Buffer Parameter Reading

The key insight is that `dpdbuf4` already contains all the information we need:

```cpp
struct dpdbuf4 {
    struct {
        int pqnum;      // Bra index combination ID
        int rsnum;      // Ket index combination ID
        int nirreps;    // Number of irreps
        // ... more metadata
    } *params;

    struct {
        int my_irrep;   // Irrep of this buffer
        char label[80]; // Human-readable label
    } file;
};
```

We can query DPD's index tables to determine:
- What orbital spaces are involved (O, o, V, v)
- How to combine them for output ([O,O], [V,V], etc.)
- What the standard label should be

This is **not guessing** - it's reading metadata that already exists!

---

## Implementation Plan

### Phase 1: Foundation (Week 1)

**Goal**: Create library infrastructure and basic utilities

**Tasks:**
1. Create file structure
   - `psi4/src/psi4/libtrans/integral_permutations.h`
   - `psi4/src/psi4/libtrans/integral_permutations.cc`
   - Update `psi4/src/psi4/libtrans/CMakeLists.txt`

2. Implement core permutation wrappers
   - `chemist_to_physicist()`
   - `physicist_to_chemist()`
   - `transpose_bra_ket()`
   - `apply_permutation()` (generic wrapper)

3. Add utility functions
   - `get_permutation_name()`
   - `is_permutation_supported()`
   - Input validation helpers

4. Write unit tests (see Testing Strategy)

**Deliverables:**
- Compiling library with 4-5 core functions
- Unit tests passing
- Doxygen documentation

### Phase 2: Semantic Wrappers (Week 2)

**Goal**: Add high-level semantic transformations

**Tasks:**
1. Implement integral space transformations
   - `ovov_to_oovv()` and spin variants
   - `oooo_permutations()` (if needed)
   - `vvvv_permutations()` (if needed)

2. Add batch operations
   - `standard_ovov_transformations()`
   - `standard_oovv_transformations()`

3. Extend unit tests for new functions

4. Write usage examples and documentation

**Deliverables:**
- 10-12 semantic wrapper functions
- Comprehensive test coverage
- Usage guide document

### Phase 3: DCT Module Refactoring (Week 2-3)

**Goal**: Refactor DCT as proof-of-concept

**Tasks:**
1. Refactor `dct_integrals_RHF.cc`
   - Replace buf4_sort calls with semantic wrappers
   - Verify tests pass (no functionality change)

2. Refactor `dct_integrals_UHF.cc`
   - Same transformations for UHF
   - Highlight opportunities for future template consolidation

3. Refactor `dct_df_operations.cc`
   - Update density-fitted integral sorting

4. Run DCT test suite
   - Verify all tests pass
   - Performance regression testing

**Deliverables:**
- DCT module using new utilities
- All DCT tests passing
- Performance parity documented

### Phase 4: CC Module Refactoring (Week 3-4)

**Goal**: Extend to CC modules

**Tasks:**
1. Refactor `cctransort/sort_tei_rhf.cc`
   - High duplication, good test case

2. Refactor `cctransort/sort_tei_uhf.cc`
   - Validate spin-aware utilities

3. Selectively refactor ccdensity files
   - Target files with >5 permutation calls
   - Leave low-usage files for future

4. Run CC test suite

**Deliverables:**
- 5-10 CC files refactored
- All CC tests passing
- Migration guide for other modules

### Phase 5: Documentation and Finalization (Week 4)

**Goal**: Complete documentation and prepare for merge

**Tasks:**
1. Write comprehensive documentation
   - Developer guide on using IntegralPermutations
   - Permutation pattern reference
   - Migration guide for other modules

2. Code review
   - Request review from Psi4 maintainers
   - Address feedback

3. Performance validation
   - Benchmark DCT and CC modules
   - Document any performance changes

4. Create PR and merge

**Deliverables:**
- Complete documentation
- Approved PR
- Merged to master

---

## Testing Strategy

### Unit Tests

**File**: `psi4/tests/pytests/test_integral_permutations.py`

**Test Cases:**

1. **Permutation Correctness**
   ```python
   def test_chemist_to_physicist():
       """Verify (pq|rs) → (pr|qs) transformation"""
       # Create small test buffer with known values
       # Apply transformation
       # Verify output matches analytical result
   ```

2. **Spin Case Handling**
   ```python
   def test_mixed_spin_transformation():
       """Verify (OV|ov) → <Oo|Vv> for α-β case"""
   ```

3. **Error Handling**
   ```python
   def test_null_buffer_error():
       """Verify proper error on null input"""
   ```

### Integration Tests

**Use existing module test suites:**

1. **DCT Tests** (`psi4/tests/dct*`)
   - Run all DCT tests before/after refactoring
   - Compare energies to 10 decimal places
   - Verify timings within 5% tolerance

2. **CC Tests** (`psi4/tests/cc*`)
   - Same verification for CC modules

3. **Regression Suite**
   - Run full Psi4 regression suite
   - No changes to any test outputs

### Performance Tests

**Benchmarks:**
```bash
# DCT performance
time psi4 dct1.in  # Before refactoring
time psi4 dct1.in  # After refactoring

# Compare timings
```

**Acceptance Criteria:**
- No performance regression >5%
- Prefer <2% variation (within noise)

### Test Matrix

| Test Type | Coverage | Pass Criteria |
|-----------|----------|---------------|
| Unit Tests | All new functions | 100% pass |
| DCT Integration | All dct* tests | Energies match to 10 decimals |
| CC Integration | All cc* tests | Energies match to 10 decimals |
| Regression | Full suite | No test failures |
| Performance | DCT, CC benchmarks | <5% regression |

---

## Timeline and Milestones

### Gantt Chart

```
Week 1: Foundation
├── Day 1-2: File structure, core wrappers
├── Day 3-4: Utility functions, validation
└── Day 5: Unit tests, documentation

Week 2: Semantic Wrappers & DCT Refactoring
├── Day 1-2: Semantic wrapper functions
├── Day 3-5: DCT RHF/UHF refactoring, testing

Week 3: DCT Completion & CC Start
├── Day 1-2: DCT DF operations, final testing
├── Day 3-5: CC module refactoring (cctransort)

Week 4: CC Completion & Finalization
├── Day 1-2: CC ccdensity refactoring
├── Day 3: Documentation
├── Day 4: Code review, address feedback
└── Day 5: Final testing, merge
```

### Milestones

| Milestone | Week | Deliverable |
|-----------|------|-------------|
| M1: Library Created | 1 | Compiling library with tests |
| M2: DCT Refactored | 2 | DCT using utilities, tests pass |
| M3: CC Refactored | 3 | cctransort using utilities |
| M4: Ready for Merge | 4 | Documentation complete, approved PR |

---

## Risk Assessment

### Risk Matrix

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Performance regression | Low | High | Benchmark early, inline functions |
| Test failures | Medium | High | Incremental testing, frequent validation |
| API mismatch | Low | Medium | Careful DPD API study, wrapper-only approach |
| Scope creep | Medium | Medium | Strict phase boundaries, defer non-critical |
| Incomplete refactoring | Low | Low | Document migration path for future |

### Specific Risks

**Risk 1: Performance Overhead from Function Calls**
- **Probability**: Low
- **Impact**: High (unacceptable if significant)
- **Mitigation**:
  - Use `inline` for small wrappers
  - Enable compiler optimization (-O3)
  - Benchmark early (Week 1)
  - Consider header-only implementation if needed

**Risk 2: DPD API Changes**
- **Probability**: Very Low
- **Impact**: High
- **Mitigation**:
  - Study DPD documentation thoroughly
  - Use only stable, public DPD APIs
  - Wrapper approach isolates changes

**Risk 3: Resistance to Adoption**
- **Probability**: Medium
- **Impact**: Medium
- **Mitigation**:
  - Demonstrate clear benefits in DCT
  - Provide migration guide
  - Don't force adoption - incremental

**Risk 4: Incomplete Testing**
- **Probability**: Low
- **Impact**: High
- **Mitigation**:
  - Use existing regression suite
  - Energy comparison to high precision
  - Multiple test systems

---

## Success Criteria

### Quantitative Metrics

1. **Code Reduction**
   - ✅ Reduce DCT permutation code by >30% (goal: 40-50 lines)
   - ✅ Reduce CC permutation code by >20% (goal: 30-40 lines)
   - ✅ Total reduction: 200+ lines across modules

2. **Test Coverage**
   - ✅ 100% of new functions have unit tests
   - ✅ All DCT tests pass with energies matching to 10 decimals
   - ✅ All CC tests pass with energies matching to 10 decimals
   - ✅ No regressions in full test suite

3. **Performance**
   - ✅ No performance regression >5% in any module
   - ✅ Ideal: <2% variation (within measurement noise)

4. **Documentation**
   - ✅ All public functions have Doxygen comments
   - ✅ Developer guide written
   - ✅ Migration guide for other modules

### Qualitative Goals

1. **Code Readability**
   - Code reviews confirm improved readability
   - Semantic function names vs cryptic permutations

2. **Maintainability**
   - Single point of change for permutation logic
   - Clear separation of concerns

3. **Foundation for Future Work**
   - Enables RHF/UHF template consolidation
   - Pattern for other refactoring efforts

---

## Future Work

### Phase 2 Enhancements (3-6 months)

1. **Template-Based Spin Case Handling**
   - Use utilities to consolidate RHF/UHF in DCT
   - Estimated: 350-400 line reduction

2. **Extend to All CC Modules**
   - Migrate remaining ccdensity, ccenergy files
   - Estimated: 100-150 line reduction

3. **Metadata-Driven Sorting**
   - Define integral spaces in metadata
   - Generate permutations from metadata
   - Long-term vision: eliminate manual permutations

### Potential Improvements

1. **Performance Optimization**
   - If overhead detected, investigate:
     - Header-only implementation
     - Template meta-programming for compile-time resolution
     - Caching of frequently-used permutations

2. **Extended Validation**
   - Add symmetry checking
   - Verify DPD index compatibility
   - Runtime validation mode for development

3. **Integration with IntegralTransform**
   - Coordinate with existing IntegralTransform class
   - Possible integration of utilities into IntegralTransform

---

## Appendices

### A. Permutation Reference

Quick reference for DPD permutation codes:

```
Input: (pq|rs)

prqs → (pr|qs)  # Chemist → Physicist
prsq → (pr|sq)  #
psqr → (ps|qr)  #
psrq → (ps|rq)  #

qprs → (qp|rs)  # Swap in bra
qpsr → (qp|sr)  #
qrps → (qr|ps)  #
qrsp → (qr|sp)  #
qspr → (qs|pr)  #
qsrp → (qs|rp)  #

rpqs → (rp|qs)  #
rpsq → (rp|sq)  #
rqps → (rq|ps)  #
rqsp → (rq|sp)  #
rsqp → (rs|qp)  #
rspq → (rs|pq)  # Swap bra↔ket

spqr → (sp|qr)  #
sprq → (sp|rq)  #
sqpr → (sq|pr)  #
sqrp → (sq|rp)  #
srqp → (sr|qp)  #
srpq → (sr|pq)  #
```

### B. Example Migration

**File**: `dct_integrals_RHF.cc`

**Lines Before**: 373 total, ~40 permutation-related
**Lines After**: ~365 total, ~32 permutation-related
**Reduction**: ~8 lines + improved clarity

---

## Approvals

| Role | Name | Signature | Date |
|------|------|-----------|------|
| Project Lead | | | |
| Psi4 Maintainer | | | |
| Code Reviewer | | | |

---

## Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2025-11-18 | Claude | Initial plan |

---

**END OF PLAN**
