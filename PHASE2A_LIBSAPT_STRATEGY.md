# Phase 2a: libsapt_solver Migration Strategy

**Date:** 2025-11-18
**Module:** libsapt_solver
**Scope:** 967 block_matrix/init_matrix calls (49% of total legacy usage)
**Status:** IN PROGRESS

---

## Executive Summary

The libsapt_solver module is the largest user of legacy matrix utilities in Psi4. This document outlines the strategic approach to migrating all 967 calls to modern SharedMatrix.

**Challenge:** Massive scope with highly complex quantum chemistry code
**Approach:** Phased migration from simple to complex files
**Estimated Effort:** 20-30 hours total

---

## File Inventory & Complexity Analysis

### Complexity Categories

Based on block_matrix/init_matrix call counts:

#### Tier 1: Simple (< 10 calls) - **5 files, 22 calls**
| File | Calls | Priority |
|------|-------|----------|
| elst12.cc | 3 | ⭐ Start here |
| elst13.cc | 3 | ⭐ Start here |
| ind-disp30.cc | 4 | ⭐ |
| utils.cc | 5 | ⭐ |
| disp20.cc | 7 | ⭐ |

**Estimated time:** 30-45 minutes total

#### Tier 2: Low Complexity (10-30 calls) - **6 files, 118 calls**
| File | Calls | Notes |
|------|-------|-------|
| sapt.cc | 12 | Main class file |
| exch11.cc | 12 | |
| disp30.cc | 16 | |
| exch10.cc | 25 | |
| disp22t.cc | 26 | |
| exch-ind30.cc | 27 | |

**Estimated time:** 1-2 hours total

#### Tier 3: Medium Complexity (30-60 calls) - **5 files, 196 calls**
| File | Calls | Notes |
|------|-------|-------|
| sapt2.cc | 33 | |
| disp22sdq.cc | 33 | |
| ind20.cc | 38 | |
| exch-disp30.cc | 38 | |
| sapt0.cc | 54 | Core SAPT0 implementation |

**Estimated time:** 2-4 hours total

#### Tier 4: High Complexity (60-110 calls) - **4 files, 373 calls**
| File | Calls | Notes |
|------|-------|-------|
| amplitudes.cc | 64 | Amplitude calculation routines |
| exch-ind20.cc | 76 | |
| exch-disp20.cc | 82 | |
| exch-ind-disp30.cc | 108 | |

**Estimated time:** 4-6 hours total

#### Tier 5: Very High Complexity (> 110 calls) - **2 files, 263 calls**
| File | Calls | Notes |
|------|-------|-------|
| disp2ccd.cc | 125 | Coupled-cluster dispersion |
| exch12.cc | 138 | **Most complex file** |

**Estimated time:** 4-6 hours total

---

## Migration Patterns Identified

### Pattern 1: Temporary Matrix for PSIO Read (Most Common)

**Before:**
```cpp
double **pAA = block_matrix(aoccA, aoccA);
psio_->read_entry(ampfile, pAAlabel, (char *)pAA[0], sizeof(double) * aoccA * aoccA);
// ... use pAA ...
free_block(pAA);
```

**After:**
```cpp
auto pAA = std::make_shared<Matrix>("pAA", aoccA, aoccA);
psio_->read_entry(ampfile, pAAlabel, (char *)pAA->get_pointer(), sizeof(double) * aoccA * aoccA);
// ... use pAA->pointer() or pAA->get_pointer() ...
// Automatic cleanup
```

### Pattern 2: Temporary Matrix for BLAS Operations

**Before:**
```cpp
double **X = block_matrix(m, n);
C_DGEMM('N', 'N', m, n, k, 1.0, A[0], lda, B[0], ldb, 0.0, X[0], ldc);
// ... use X ...
free_block(X);
```

**After:**
```cpp
auto X = std::make_shared<Matrix>("X", m, n);
C_DGEMM('N', 'N', m, n, k, 1.0, A[0], lda, B[0], ldb, 0.0, X->get_pointer(), ldc);
// ... use X->pointer() ...
// Automatic cleanup
```

### Pattern 3: Matrix Used in Tight Loop

**Before:**
```cpp
for (int h = 0; h < nirrep; h++) {
    double **T = block_matrix(rows[h], cols[h]);
    // ... operations ...
    free_block(T);
}
```

**After:**
```cpp
for (int h = 0; h < nirrep; h++) {
    auto T = std::make_shared<Matrix>("T", rows[h], cols[h]);
    // ... operations with T->pointer() ...
    // Automatic cleanup at end of each iteration
}
```

---

## Migration Approach

### Phase 1: Simple Files (Tier 1)
**Goal:** Establish pattern, build confidence
**Files:** elst12.cc, elst13.cc, ind-disp30.cc, utils.cc, disp20.cc
**Calls:** 22
**Time:** 45 minutes

**Steps:**
1. Migrate elst13.cc first (3 calls, already analyzed)
2. Verify pattern works
3. Apply to remaining Tier 1 files
4. Commit batch

### Phase 2: Low Complexity (Tier 2)
**Goal:** Handle moderate volume
**Files:** 6 files
**Calls:** 118
**Time:** 1-2 hours

**Steps:**
1. Migrate files in order of complexity
2. Commit every 2-3 files
3. Test incrementally

### Phase 3: Medium Complexity (Tier 3)
**Goal:** Core SAPT functionality
**Files:** 5 files including sapt0.cc
**Calls:** 196
**Time:** 2-4 hours

**Steps:**
1. Extra care with sapt0.cc (most used)
2. Test sapt0 calculations specifically
3. Commit per file for easy rollback

### Phase 4: High Complexity (Tier 4)
**Goal:** Complex amplitude routines
**Files:** 4 files
**Calls:** 373
**Time:** 4-6 hours

**Steps:**
1. Analyze each file's unique patterns
2. Migrate conservatively
3. Extensive testing

### Phase 5: Very High Complexity (Tier 5)
**Goal:** Most complex files
**Files:** exch12.cc, disp2ccd.cc
**Calls:** 263
**Time:** 4-6 hours

**Steps:**
1. Deep analysis before migration
2. Consider splitting into multiple commits
3. Comprehensive testing

---

## Risk Mitigation

### Testing Strategy

**Per-Tier Testing:**
- Tier 1: Quick compilation check
- Tier 2-3: Run basic SAPT tests
- Tier 4-5: Full SAPT test suite

**Regression Tests:**
```bash
pytest tests/pytests/test_sapt*
pytest tests/sapt*
```

### Rollback Plan

Each tier committed separately:
- Tier 1: Single commit (low risk)
- Tier 2: Commit every 2-3 files
- Tier 3+: Commit per file

### Known Challenges

1. **PSIO Read/Write:** Requires raw pointers
   - Solution: Use `->get_pointer()` for flat arrays

2. **BLAS/LAPACK:** Requires contiguous memory
   - Solution: Use `->pointer()[0]` or `->get_pointer()`

3. **Large file size:** exch12.cc is 67KB
   - Solution: Work methodically, use search/replace carefully

---

## Success Criteria

- ✅ All 967 calls migrated to SharedMatrix
- ✅ Full SAPT test suite passes
- ✅ No performance regression (< 1%)
- ✅ No memory leaks (valgrind clean)
- ✅ Code compiles without warnings

---

## Progress Tracking

| Tier | Files | Calls | Status | Commits |
|------|-------|-------|--------|---------|
| Tier 1 | 5 | 22 | ⏳ Pending | - |
| Tier 2 | 6 | 118 | ⏳ Pending | - |
| Tier 3 | 5 | 196 | ⏳ Pending | - |
| Tier 4 | 4 | 373 | ⏳ Pending | - |
| Tier 5 | 2 | 263 | ⏳ Pending | - |
| **Total** | **22** | **972** | **0%** | **0** |

*Note: Total is 972 vs reported 967 due to grep count variations*

---

## Timeline Estimate

**Optimistic:** 15 hours (3 days part-time)
**Realistic:** 20-25 hours (4-5 days part-time)
**Pessimistic:** 30 hours (6 days part-time)

**Assumes:**
- 1 developer
- Good test coverage available
- No major blockers discovered

---

## Next Steps

1. ✅ Complete this strategy document
2. ⏳ Begin Tier 1 migration (elst13.cc first)
3. ⏳ Validate pattern with testing
4. ⏳ Proceed tier by tier

---

**Document Status:** APPROVED - Ready for implementation
**Start Date:** 2025-11-18
**Target Completion:** 2025-11-25 (1 week)
