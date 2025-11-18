# Project Plan: Extract Common Permutation Utilities

## Executive Summary

**Project Goal**: Extract and consolidate duplicated integral permutation code across Psi4 modules into a reusable utility library.

**Impact**:
- Reduce ~200-300 lines of duplicated permutation logic
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
 * All methods are static - no instantiation required.
 */
class IntegralPermutations {
public:
    //===========================================
    // Notation Conversions
    //===========================================

    /**
     * @brief Convert chemist's notation (pq|rs) to physicist's notation <pr|qs>
     *
     * Applies the prqs permutation: (pq|rs) → (pr|qs)
     * This is the most common transformation in quantum chemistry codes.
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

    /**
     * @brief Convert physicist's notation <pr|qs> to chemist's notation (pq|rs)
     *
     * Inverse of chemist_to_physicist. Uses prqs permutation (self-inverse for this case).
     *
     * @param InBuf Input buffer in physicist's notation
     * @param outfilenum PSI file number for output
     * @param pq_indices DPD indices for compound index (pq)
     * @param rs_indices DPD indices for compound index (rs)
     * @param label String label for output buffer
     */
    static void physicist_to_chemist(
        dpdbuf4 *InBuf,
        int outfilenum,
        int pq_indices,
        int rs_indices,
        const std::string &label
    );

    //===========================================
    // Common Integral Space Transformations
    //===========================================

    /**
     * @brief Transform (OV|OV) → <OO|VV>
     *
     * Common transformation for building 2-electron integrals in
     * occupied-occupied/virtual-virtual space.
     *
     * @param InBuf Input (OV|OV) buffer
     * @param outfilenum Output file
     * @param oo_indices DPD indices for OO space
     * @param vv_indices DPD indices for VV space
     * @param label Output label
     */
    static void ovov_to_oovv(
        dpdbuf4 *InBuf,
        int outfilenum,
        int oo_indices,
        int vv_indices,
        const std::string &label
    );

    /**
     * @brief Create transpose: <AB|CD> → <CD|AB>
     *
     * Swaps bra and ket pairs using rspq permutation.
     *
     * @param InBuf Input buffer
     * @param outfilenum Output file
     * @param cd_indices DPD indices for new bra (CD)
     * @param ab_indices DPD indices for new ket (AB)
     * @param label Output label
     */
    static void transpose_bra_ket(
        dpdbuf4 *InBuf,
        int outfilenum,
        int cd_indices,
        int ab_indices,
        const std::string &label
    );

    //===========================================
    // Spin-Aware Transformations (UHF/ROHF)
    //===========================================

    /**
     * @brief Transform (OV|ov) → <Oo|Vv> for mixed-spin case
     *
     * Handles α-β spin case transformation.
     *
     * @param InBuf Input buffer
     * @param outfilenum Output file
     * @param oo_indices DPD indices for Oo
     * @param vv_indices DPD indices for Vv
     * @param label Output label
     */
    static void ovov_to_oovv_mixed_spin(
        dpdbuf4 *InBuf,
        int outfilenum,
        int oo_indices,
        int vv_indices,
        const std::string &label
    );

    //===========================================
    // Generic Permutation Interface
    //===========================================

    /**
     * @brief Apply arbitrary permutation with semantic naming
     *
     * Thin wrapper around DPD buf4_sort with better error checking
     * and documentation.
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
    // Batch Operations
    //===========================================

    /**
     * @brief Perform standard OVOV transformations in batch
     *
     * Common pattern: Generate all permutations needed for OVOV integrals:
     * - (OV|OV) → <OO|VV>
     * - <OO|VV> → <VV|OO>
     * - Additional requested permutations
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

namespace psi {
namespace libtrans {

void IntegralPermutations::chemist_to_physicist(
    dpdbuf4 *InBuf,
    int outfilenum,
    int pq_indices,
    int rs_indices,
    const std::string &label)
{
    // Validate input
    if (InBuf == nullptr) {
        throw std::runtime_error("IntegralPermutations::chemist_to_physicist: null input buffer");
    }

    // Perform transformation using DPD library
    // (pq|rs) → (pr|qs) using prqs permutation
    global_dpd_->buf4_sort(InBuf, outfilenum, prqs, pq_indices, rs_indices, label);
}

void IntegralPermutations::ovov_to_oovv(
    dpdbuf4 *InBuf,
    int outfilenum,
    int oo_indices,
    int vv_indices,
    const std::string &label)
{
    // This is semantically equivalent to chemist_to_physicist
    // but provides a more descriptive name for this common operation
    chemist_to_physicist(InBuf, outfilenum, oo_indices, vv_indices, label);
}

void IntegralPermutations::transpose_bra_ket(
    dpdbuf4 *InBuf,
    int outfilenum,
    int cd_indices,
    int ab_indices,
    const std::string &label)
{
    global_dpd_->buf4_sort(InBuf, outfilenum, rspq, cd_indices, ab_indices, label);
}

void IntegralPermutations::standard_ovov_transformations(
    dpdbuf4 *InBuf,
    int outfilenum,
    bool generate_transpose)
{
    // Extract DPD indices from input buffer
    // This is a simplified example - actual implementation would need more detail

    // Transform (OV|OV) → <OO|VV>
    ovov_to_oovv(InBuf, outfilenum,
                 /* determine OO indices */,
                 /* determine VV indices */,
                 "MO Ints <OO|VV>");

    if (generate_transpose) {
        // Open the just-created <OO|VV> buffer
        dpdbuf4 OOVVBuf;
        // ... initialize buffer ...

        // Create <VV|OO> transpose
        transpose_bra_ket(&OOVVBuf, outfilenum,
                         /* VV indices */,
                         /* OO indices */,
                         "MO Ints <VV|OO>");

        global_dpd_->buf4_close(&OOVVBuf);
    }
}

std::string IntegralPermutations::get_permutation_name(indices permutation) {
    switch(permutation) {
        case prqs: return "Chemist to Physicist notation";
        case rspq: return "Transpose bra and ket";
        case psrq: return "Mixed permutation (ps|rq)";
        case qprs: return "Swap within bra";
        case pqsr: return "Swap within ket";
        // ... more cases ...
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

    // Semantic, self-documenting function names
    IntegralPermutations::ovov_to_oovv(&I, PSIF_LIBTRANS_DPD,
                                       ID("[O,O]"), ID("[V,V]"),
                                       "MO Ints <OO|VV>");

    IntegralPermutations::apply_permutation(&I, PSIF_LIBTRANS_DPD, psrq,
                                           ID("[O,V]"), ID("[O,V]"),
                                           "MO Ints SF <OV|OV>:<Ov|oV>");

    global_dpd_->buf4_close(&I);
}
```

**Benefits:**
- Clear intent: `ovov_to_oovv` vs cryptic `prqs`
- Centralized: One place to update if DPD API changes
- Documented: Function documentation explains transformation
- Validated: Input checking and error messages

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
