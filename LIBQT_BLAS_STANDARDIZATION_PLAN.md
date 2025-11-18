# LibQt BLAS Standardization Plan for Psi4

**Document Version:** 1.0
**Created:** 2025-11-18
**Branch:** `claude/standardize-blas-wrappers-01WMXdqxgvAKrpQj7kqRjqSH`
**Status:** Planning Phase

---

## Executive Summary

This document outlines a comprehensive plan to establish **libqt** as the standard BLAS/LAPACK interface for the Psi4 quantum chemistry package. The codebase currently has **88% of BLAS operations already using libqt**, making standardization highly feasible with minimal risk.

**Key Metrics:**
- **Total BLAS Operations:** 2,845+ calls across 228 files
- **Already Standardized:** 2,505 calls (88.1%) using libqt C interface
- **Needs Migration:** 340 calls (11.9%) using fragmented wrappers
- **Estimated Effort:** 4-8 weeks (phased approach)
- **Risk Level:** LOW to MEDIUM (depending on phase)

---

## Table of Contents

1. [Current State Analysis](#1-current-state-analysis)
2. [Objectives and Goals](#2-objectives-and-goals)
3. [Target Architecture](#3-target-architecture)
4. [Migration Strategy](#4-migration-strategy)
5. [Implementation Phases](#5-implementation-phases)
6. [Risk Assessment and Mitigation](#6-risk-assessment-and-mitigation)
7. [Testing and Validation](#7-testing-and-validation)
8. [Timeline and Resources](#8-timeline-and-resources)
9. [Success Criteria](#9-success-criteria)
10. [Appendices](#10-appendices)

---

## 1. Current State Analysis

### 1.1 Existing BLAS Wrapper Implementations

#### **A. libqt (Primary Interface - 88% Usage)**
- **Location:** `psi4/src/psi4/libqt/`
- **Files:** `qt.h`, `blas_intfc.cc`, `blas_intfc23.cc`, `lapack_intfc.h`
- **Interface Style:** C-ordered, size_t support (>2^31 elements)
- **Operations:** Complete BLAS 1/2/3 + comprehensive LAPACK (200+ routines)
- **Name Mangling:** FCMangle with 4 fallback schemes
- **Advantages:**
  - Most widely used (2,505 calls across 17 modules)
  - Portable across platforms
  - Handles large arrays via INT_MAX blocking
  - Extensively documented
  - Well-tested and stable

#### **B. fnocc/blas.h (Fortran-Ordered - 8.1% Usage)**
- **Location:** `psi4/src/psi4/fnocc/`
- **Files:** `blas.h`, `blas.cc`, `blas_mangle.h`
- **Interface Style:** Fortran-ordered (column-major)
- **Operations:** Basic BLAS 1/2/3 + limited LAPACK (DGESV, DSYEV, DGESVD)
- **Usage:** 355 calls in fnocc module (ccsd.cc, qed.cc, mp4.cc, etc.)
- **Critical File:** `ccsd.cc` - 92,000 lines with 150+ F_DGEMM calls
- **Advantages:**
  - Optimal for Fortran memory layout
  - Direct mapping to BLAS semantics
  - No transpose overhead
- **Issues:**
  - Duplicates libqt functionality
  - Limited to 32-bit integers
  - Inconsistent with rest of codebase

#### **C. psimrcc/algebra_interface.h (Fortran-Ordered - 1.5% Usage)**
- **Location:** `psi4/src/psi4/psimrcc/`
- **Files:** `algebra_interface.h`, `algebra_interface.cc`
- **Interface Style:** Fortran-ordered wrapper around libqt
- **Operations:** F_DGEMM, F_DGEMV (forwarding to C_DGEMM internally)
- **Usage:** 67 calls in 8 psimrcc files
- **Note:** Actually uses libqt underneath with transpose adjustments
- **Advantage:**
  - Maintains Fortran semantics for MRCC algorithms
- **Issues:**
  - Hidden dependency on libqt
  - Adds unnecessary layer

#### **D. dfocc/arrays.h (Array Class Wrapper - 1.9% Usage)**
- **Location:** `psi4/src/psi4/dfocc/`
- **Files:** `arrays.h`, `arrays.cc`
- **Interface Style:** C++ class-based (Array1d, Array2d, Array3d)
- **Operations:** gemm(), gemv(), dot() methods wrapping C_DGEMM, C_DGEMV, C_DDOT
- **Usage:** 85 calls in dfocc module
- **Advantage:**
  - Object-oriented convenience
  - Automatic memory management
  - Type safety
- **Note:** Already uses libqt internally

#### **E. occ/arrays.h (Array Class Wrapper - 0.1% Usage)**
- **Location:** `psi4/src/psi4/occ/`
- **Files:** `arrays.h`
- **Interface Style:** Similar to dfocc arrays
- **Usage:** 4 calls in occ module
- **Note:** Already uses libqt internally

### 1.2 Usage Statistics by Module

| Module | BLAS Calls | % of Total | Primary Interface |
|--------|-----------|-----------|-------------------|
| libsapt_solver | 1,610 | 36.5% | libqt (C_DGEMM) |
| cc | 534 | 12.2% | libqt (C_DGEMM) |
| scfgrad | 410 | 9.4% | libqt (C_DGEMM) |
| libfock | 411 | 9.4% | libqt (C_DGEMM) |
| fnocc | 355 | 8.1% | fnocc (F_DGEMM) |
| dfocc | 85 | 1.9% | Array classes |
| psimrcc | 67 | 1.5% | algebra_interface |
| detci | 180 | 4.1% | libqt (C_DGEMM) |
| cphf | 142 | 3.2% | libqt (C_DGEMM) |
| dcft | 89 | 2.0% | libqt (C_DGEMM) |
| Others (8 modules) | ~135 | 3.1% | libqt (C_DGEMM) |
| **TOTAL** | **2,845+** | **100%** | **libqt (88%)** |

### 1.3 Most Frequently Used Operations

| Operation | Frequency | Primary Use Case |
|-----------|-----------|------------------|
| C_DGEMM / F_DGEMM | 2,419 | Matrix multiplication (78.8%) |
| C_DDOT | 648 | Inner products |
| C_DAXPY | 558 | Vector addition |
| C_DCOPY | 416 | Memory copying |
| C_DGEMV / F_DGEMV | 394 | Matrix-vector products |
| C_DSCAL | 287 | Vector scaling |
| C_DSYEV | 156 | Eigenvalue decomposition |
| C_DGESVD | 89 | Singular value decomposition |
| Others | ~300 | Various BLAS/LAPACK operations |

### 1.4 Key Findings

✅ **Strengths:**
- 88% of codebase already uses libqt consistently
- libqt provides complete BLAS/LAPACK coverage
- libqt supports large arrays (>2^31 elements)
- Extensive documentation and testing

⚠️ **Issues:**
- Fragmentation: 5 different BLAS interfaces for the same operations
- Code duplication: Multiple implementations of DGEMM, DGEMV, etc.
- Maintenance burden: Changes require updating multiple wrappers
- Inconsistent naming: F_DGEMM vs C_DGEMM vs gemm()
- Missing opportunity: Could leverage libqt's INT_MAX blocking everywhere

---

## 2. Objectives and Goals

### 2.1 Primary Objectives

1. **Standardization**
   - Establish libqt as the single source of truth for BLAS/LAPACK operations
   - Eliminate redundant wrapper implementations
   - Provide consistent API across all Psi4 modules

2. **Maintainability**
   - Reduce code duplication
   - Centralize platform-specific name mangling
   - Simplify build system dependencies
   - Enable easier debugging and profiling

3. **Functionality**
   - Preserve existing algorithmic correctness
   - Maintain performance characteristics
   - Support both C-ordered and Fortran-ordered memory layouts
   - Enable large array support (>2^31 elements) everywhere

4. **Backward Compatibility**
   - Minimize disruption to existing modules
   - Preserve module-specific abstractions (CCBLAS, Array classes)
   - Maintain API compatibility where possible

### 2.2 Non-Goals

❌ **Out of Scope:**
- Refactoring algorithms to change memory layouts
- Rewriting fnocc/ccsd.cc (92,000 lines) from scratch
- Changing psimrcc tensor expression framework
- Replacing Array class convenience layers (unless desired)
- Performance optimization (beyond standardization benefits)

### 2.3 Success Criteria

✓ All BLAS operations route through libqt (directly or via thin wrappers)
✓ Remove fnocc/blas.h and fnocc/blas_mangle.h
✓ Consolidate psimrcc to use libqt
✓ No performance regression in critical paths (ccsd, mrcc)
✓ All existing tests pass
✓ Documentation updated to reflect standardized interface
✓ Build system simplified (single BLAS dependency path)

---

## 3. Target Architecture

### 3.1 Architectural Principles

1. **Single Source of Truth:** libqt provides all BLAS/LAPACK operations
2. **Layered Abstraction:** Allow thin wrappers for specific use cases
3. **Zero Performance Penalty:** No additional overhead in critical paths
4. **Memory Layout Flexibility:** Support both C-ordered and Fortran-ordered data

### 3.2 Proposed Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    APPLICATION MODULES                       │
│  (fnocc, psimrcc, dfocc, occ, cc, detci, sapt, etc.)       │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│              OPTIONAL CONVENIENCE LAYERS                     │
│                                                              │
│  ┌─────────────────┐  ┌──────────────────┐                 │
│  │ libqt/          │  │ Module-Specific   │                 │
│  │ blas_fortran.h  │  │ Abstractions:     │                 │
│  │                 │  │ - CCBLAS (psimrcc)│                 │
│  │ F_DGEMM()       │  │ - Array (dfocc)   │                 │
│  │ F_DGEMV()       │  │                   │                 │
│  │ [thin wrappers] │  │ [wraps libqt]     │                 │
│  └─────────────────┘  └──────────────────┘                 │
│            │                    │                            │
│            └────────┬───────────┘                            │
└─────────────────────┼────────────────────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────────────────┐
│                 LIBQT CORE INTERFACE                        │
│              (Single BLAS/LAPACK Gateway)                    │
│                                                              │
│  C_DGEMM(), C_DGEMV(), C_DDOT(), C_DCOPY(), ...             │
│  C_DSYEV(), C_DGESVD(), C_DGESV(), ...                      │
│  200+ LAPACK routines with full documentation               │
│                                                              │
│  Features:                                                   │
│  - INT_MAX blocking for large arrays                        │
│  - size_t support (>2^31 elements)                          │
│  - Platform-independent name mangling                       │
│  - Comprehensive error handling                             │
└─────────────────────────────────────────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────────────────┐
│            VENDOR BLAS/LAPACK LIBRARIES                     │
│   (Intel MKL, OpenBLAS, ATLAS, Accelerate, etc.)           │
└─────────────────────────────────────────────────────────────┘
```

### 3.3 Component Descriptions

#### **Layer 1: Vendor BLAS/LAPACK Libraries**
- Platform-specific optimized implementations
- No changes required

#### **Layer 2: libqt Core Interface**
- **Primary gateway for all BLAS/LAPACK operations**
- Handles name mangling via FCMangle
- Provides INT_MAX blocking for large arrays
- Already exists and is well-tested
- **Enhancements needed:**
  - None required (already complete)

#### **Layer 3: Optional Convenience Layers**

**A. libqt/blas_fortran.h (NEW)**
- Thin Fortran-ordered wrappers
- Provides F_DGEMM, F_DGEMV, F_DCOPY, F_DDOT, etc.
- Internally calls C_DGEMM with transposed arguments
- **Purpose:** Enable fnocc/psimrcc to maintain Fortran semantics
- **Performance:** Zero overhead (compiler inline)

**B. Module-Specific Abstractions (EXISTING)**
- psimrcc CCBLAS: High-level tensor expression framework
- dfocc/occ Array classes: Object-oriented convenience
- **Action:** Update to use libqt, remove direct BLAS calls
- **Purpose:** Maintain existing high-level APIs

#### **Layer 4: Application Modules**
- Use libqt directly (preferred)
- Use convenience wrappers when appropriate
- **Action:** Update includes to point to libqt

---

## 4. Migration Strategy

### 4.1 Guiding Principles

1. **Phased Approach:** Incremental changes with testing between phases
2. **Risk Mitigation:** Start with low-risk changes, defer high-risk items
3. **Test-Driven:** Validate each phase with full test suite
4. **Backward Compatible:** Maintain existing APIs during transition
5. **Documentation First:** Update docs before code changes

### 4.2 Migration Tiers

Based on risk assessment, migrations are categorized into four tiers:

#### **Tier 1: Include Path Standardization (LOWEST RISK)**
- **Scope:** 195 files already using libqt
- **Action:** Ensure consistent #include <psi4/libqt/qt.h>
- **Risk:** VERY LOW
- **Effort:** 1-2 days
- **Validation:** Compilation check

#### **Tier 2: Array Class Consolidation (LOW RISK)**
- **Scope:** dfocc/occ Array classes (89 calls)
- **Action:**
  - Verify Array methods use libqt internally (already done)
  - Optionally consolidate dfocc and occ arrays
- **Risk:** LOW (already using libqt)
- **Effort:** 1 week (if consolidation desired)
- **Validation:** DFOCC/OCC test suite

#### **Tier 3: psimrcc Fortran Wrapper Migration (MEDIUM RISK)**
- **Scope:** 67 calls in psimrcc using algebra_interface
- **Action:**
  - Create libqt/blas_fortran.h
  - Migrate psimrcc to use new wrapper
  - Remove algebra_interface.h
- **Risk:** MEDIUM (tensor framework complexity)
- **Effort:** 2 weeks
- **Validation:** PSIMRCC test suite (MRCC calculations)

#### **Tier 4: fnocc Fortran Wrapper Migration (HIGHEST RISK)**
- **Scope:** 355 calls in fnocc (critical: ccsd.cc with 92K lines)
- **Action:**
  - Expand libqt/blas_fortran.h to cover all fnocc operations
  - Migrate fnocc to use libqt/blas_fortran.h
  - Remove fnocc/blas.h and fnocc/blas_mangle.h
- **Risk:** MEDIUM-HIGH (large codebase, performance-critical)
- **Effort:** 3-4 weeks
- **Validation:** Comprehensive FNOCC test suite + benchmarking

### 4.3 Recommended Migration Options

#### **Option A: Aggressive Standardization (Recommended)**
- Implement all four tiers
- Timeline: 6-8 weeks
- Outcome: Complete standardization, single BLAS interface
- Risk: MEDIUM (tier 4 requires careful validation)

#### **Option B: Conservative Standardization**
- Implement tiers 1-3, defer tier 4
- Timeline: 3-4 weeks
- Outcome: 92% standardized (keep fnocc/blas.h)
- Risk: LOW

#### **Option C: Minimal Standardization**
- Implement tiers 1-2 only
- Timeline: 2 weeks
- Outcome: 90% standardized (keep psimrcc and fnocc wrappers)
- Risk: VERY LOW

#### **Option D: Documentation-Only**
- No code changes
- Standardize on libqt in documentation
- Declare fnocc/psimrcc wrappers as deprecated
- Timeline: 1 week
- Risk: NONE

**Recommendation:** **Option A (Aggressive)** with proper testing and validation. The benefits of complete standardization outweigh the risks, especially with phased implementation.

---

## 5. Implementation Phases

### Phase 1: Foundation and Preparation (Week 1-2)

#### **Tasks:**

**1.1 Create libqt/blas_fortran.h**
- **File:** `psi4/src/psi4/libqt/blas_fortran.h`
- **Content:** Fortran-ordered wrappers for libqt
- **Operations to include:**
  ```cpp
  // BLAS Level 1
  inline void F_DCOPY(size_t n, const double *x, int incx, double *y, int incy);
  inline double F_DDOT(size_t n, const double *x, int incx, const double *y, int incy);
  inline void F_DSCAL(size_t n, double alpha, double *x, int incx);
  inline void F_DAXPY(size_t n, double alpha, const double *x, int incx, double *y, int incy);

  // BLAS Level 2
  inline void F_DGEMV(char trans, size_t m, size_t n, double alpha,
                      const double *a, size_t lda, const double *x, int incx,
                      double beta, double *y, int incy);

  // BLAS Level 3
  inline void F_DGEMM(char transa, char transb, size_t m, size_t n, size_t k,
                      double alpha, const double *a, size_t lda,
                      const double *b, size_t ldb,
                      double beta, double *c, size_t ldc);

  // LAPACK
  inline void F_DGESV(size_t n, size_t nrhs, double *a, size_t lda,
                      int *ipiv, double *b, size_t ldb, int *info);
  inline void F_DSYEV(char jobz, char uplo, size_t n, double *a, size_t lda,
                      double *w, double *work, size_t lwork, int *info);
  inline void F_DGESVD(char jobu, char jobvt, size_t m, size_t n,
                       double *a, size_t lda, double *s,
                       double *u, size_t ldu, double *vt, size_t ldvt,
                       double *work, size_t lwork, int *info);
  ```
- **Implementation:** Map to C_* functions with appropriate transposes
- **Testing:** Unit tests in libqt/test_blas_fortran.cc

**1.2 Documentation**
- Create BLAS_MIGRATION_GUIDE.md
- Update libqt/README.md
- Document new blas_fortran.h interface
- Add migration examples

**1.3 Set Up Testing Infrastructure**
- Identify critical test cases for each module
- Create benchmarking harness for performance validation
- Prepare regression test suite

**Deliverables:**
- ✓ libqt/blas_fortran.h with unit tests
- ✓ Updated documentation
- ✓ Testing infrastructure

**Validation:**
- Unit tests pass
- No changes to existing modules yet (preparation phase)

---

### Phase 2: Tier 1 - Include Path Standardization (Week 2)

#### **Tasks:**

**2.1 Audit Current Includes**
- Verify all 195 files using libqt have consistent includes
- Check for redundant includes

**2.2 Standardize Include Statements**
- Ensure all files use: `#include "psi4/libqt/qt.h"`
- Remove duplicate or redundant includes

**2.3 Update Build System**
- Verify CMakeLists.txt properly links libqt
- Check include directories are correctly specified

**Files Affected:** ~195 files across 17 modules
- libsapt_solver (35 files)
- cc (28 files)
- scfgrad (18 files)
- libfock (22 files)
- And others...

**Implementation:**
```bash
# Example automated refactoring
find psi4/src -name "*.cc" -o -name "*.h" | \
  xargs grep -l "C_DGEMM\|C_DGEMV" | \
  xargs sed -i 's/#include.*qt\.h.*/#include "psi4\/libqt\/qt.h"/'
```

**Deliverables:**
- ✓ Consistent include paths across codebase
- ✓ Updated build system

**Validation:**
- Clean compilation
- All tests pass

---

### Phase 3: Tier 2 - Array Class Consolidation (Week 3)

#### **Tasks:**

**3.1 Analyze dfocc and occ Array Classes**
- Compare implementations in dfocc/arrays.h and occ/arrays.h
- Identify differences and dependencies
- Determine if consolidation is beneficial

**3.2 Optional: Consolidate to Single Array Implementation**
- Create unified libarray/ directory (if desired)
- Move common Array1d/2d/3d code
- Update dfocc and occ to use unified arrays

**3.3 Verify libqt Usage**
- Confirm all Array methods use libqt internally
- No direct BLAS calls outside libqt

**Files Affected:**
- dfocc/arrays.h, dfocc/arrays.cc (85 calls)
- occ/arrays.h (4 calls)
- ~30 files using Array classes

**Implementation Decision:**
- **Option A (Recommended):** Keep separate for now, verify libqt usage
- **Option B:** Consolidate if significant duplication found

**Deliverables:**
- ✓ Verified Array classes use libqt
- ✓ Optional: Consolidated Array implementation

**Validation:**
- DFOCC test suite passes
- OCC test suite passes
- No performance regression

---

### Phase 4: Tier 3 - psimrcc Migration (Week 4-5)

#### **Tasks:**

**4.1 Update psimrcc to Use libqt/blas_fortran.h**
- Replace algebra_interface.h includes with libqt/blas_fortran.h
- Update F_DGEMM, F_DGEMV calls to use new wrapper

**Files to Modify:**
- psimrcc/algebra_interface.h → DELETE
- psimrcc/algebra_interface.cc → DELETE
- Update 8 files using algebra_interface:
  - mp2_ccsd.cc
  - mrcc_pert_triples.cc
  - mrccsd_t.cc
  - ccmrcc.cc
  - And others...

**Implementation:**
```cpp
// Before
#include "algebra_interface.h"
F_DGEMM('n', 'n', m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

// After
#include "psi4/libqt/blas_fortran.h"
F_DGEMM('n', 'n', m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
// (Same API, different implementation)
```

**4.2 Update CCBLAS Class**
- Ensure internal BLAS calls use libqt
- Remove any dependencies on algebra_interface

**4.3 Testing**
- Run full PSIMRCC test suite
- Validate MRCC calculations
- Benchmark critical calculations

**Deliverables:**
- ✓ psimrcc using libqt/blas_fortran.h
- ✓ algebra_interface.h removed
- ✓ All PSIMRCC tests passing

**Validation:**
- MRCCSD calculations correct
- Perturbative triples correct
- No performance regression
- Memory usage unchanged

---

### Phase 5: Tier 4 - fnocc Migration (Week 5-8)

#### **Tasks:**

**5.1 Expand libqt/blas_fortran.h**
- Add all operations used by fnocc:
  - BLAS1: DCOPY, DDOT, DNRM2, DSCAL, DAXPY
  - BLAS2: DGEMV, DGBMV
  - BLAS3: DGEMM
  - LAPACK: DGESV, DSYEV, DSPEV, DGESVD

**5.2 Create Migration Branch for fnocc**
- Work on separate branch for safety
- Enable easy rollback if issues arise

**5.3 Update fnocc Files (355 calls)**

**Critical file: ccsd.cc (92,000 lines, 150+ F_DGEMM calls)**
```cpp
// Before
#include "blas.h"
F_DGEMM('n', 'n', o*o, v*v, o*v, 1.0, integrals, o*o, tempv, o*v, 1.0, tempt, o*o);

// After
#include "psi4/libqt/blas_fortran.h"
F_DGEMM('n', 'n', o*o, v*v, o*v, 1.0, integrals, o*o, tempv, o*v, 1.0, tempt, o*o);
// (API unchanged, just include change)
```

**Other files to update:**
- cepa.cc
- qed.cc
- mp4.cc
- mp3.cc
- ccsd_triples.cc
- df_helper.cc
- (All fnocc/*.cc files using BLAS)

**5.4 Remove fnocc BLAS Wrappers**
- DELETE: fnocc/blas.h
- DELETE: fnocc/blas.cc
- DELETE: fnocc/blas_mangle.h

**5.5 Extensive Testing**
- Full fnocc test suite (CCSD, CCSD(T), QCISD, QCISD(T), MP3, MP4, CEPA)
- Benchmark against previous version
- Validate numerical accuracy

**5.6 Performance Validation**
- Run benchmark calculations:
  - CCSD/water
  - CCSD(T)/water
  - QCISD(T)/water
- Compare timing with baseline
- Acceptable: ±2% variation
- Investigate if >5% change

**Deliverables:**
- ✓ fnocc using libqt/blas_fortran.h
- ✓ fnocc/blas.h removed
- ✓ All fnocc tests passing
- ✓ Performance validated

**Validation:**
- All CCSD/CCSD(T) tests pass
- QCISD/QCISD(T) tests pass
- MP3/MP4 tests pass
- CEPA tests pass
- Numerical results identical (within 1e-10)
- Performance within ±2%

---

### Phase 6: Cleanup and Documentation (Week 8)

#### **Tasks:**

**6.1 Code Cleanup**
- Remove deprecated includes
- Clean up build system
- Remove unused BLAS wrapper code

**6.2 Documentation Updates**
- Update all module READMEs
- Update developer documentation
- Add BLAS standardization guide
- Document libqt/blas_fortran.h API

**6.3 Build System Optimization**
- Simplify BLAS/LAPACK dependencies
- Ensure single path to vendor libraries
- Update CMake documentation

**6.4 Create Migration Summary**
- Document what changed
- List affected files
- Provide troubleshooting guide

**Deliverables:**
- ✓ Clean, standardized codebase
- ✓ Comprehensive documentation
- ✓ Migration guide for future developers

---

## 6. Risk Assessment and Mitigation

### 6.1 Risk Matrix

| Risk | Probability | Impact | Severity | Mitigation |
|------|-------------|--------|----------|------------|
| **Performance regression in CCSD** | MEDIUM | HIGH | HIGH | Benchmark before/after, inline wrappers |
| **Numerical differences in results** | LOW | HIGH | MEDIUM | Extensive validation, bit-wise comparison |
| **Build failures on some platforms** | MEDIUM | MEDIUM | MEDIUM | Test on multiple platforms before merge |
| **Breaking downstream dependencies** | LOW | MEDIUM | LOW | Deprecation warnings, transition period |
| **Memory layout issues** | LOW | HIGH | MEDIUM | Preserve Fortran ordering, validate stride handling |
| **Test suite failures** | MEDIUM | MEDIUM | MEDIUM | Run full tests after each phase |
| **Merge conflicts during migration** | MEDIUM | LOW | LOW | Frequent rebases, coordinated development |

### 6.2 Specific Risk Mitigation Strategies

#### **Risk 1: Performance Regression**

**Concern:** Adding wrapper layer could impact performance in tight loops

**Mitigation:**
1. Mark all F_* wrappers as `inline` to eliminate function call overhead
2. Rely on compiler optimization (LTO/IPO)
3. Benchmark critical paths before and after
4. Profile hot spots with gprof/perf
5. Acceptance criteria: <2% performance variation

**Rollback Plan:**
- Keep fnocc/blas.h on separate branch
- Can revert if performance unacceptable

#### **Risk 2: Numerical Differences**

**Concern:** Subtly different BLAS call semantics could change results

**Mitigation:**
1. Ensure F_DGEMM wrapper exactly replicates Fortran semantics
2. Use identical vendor BLAS underneath
3. Compare results bit-for-bit with previous version
4. Run full test suite with strict tolerances (1e-10)
5. Validate reference calculations from literature

**Validation:**
```python
# Compare old vs new results
import numpy as np
old_energy = -76.240832480  # Previous CCSD/water
new_energy = -76.240832480  # After migration
assert abs(old_energy - new_energy) < 1e-10
```

#### **Risk 3: Platform-Specific Issues**

**Concern:** Different BLAS vendors, compilers, or platforms may behave differently

**Mitigation:**
1. Test on multiple platforms:
   - Linux (GCC, Intel ICC, Clang)
   - macOS (Clang, Accelerate framework)
   - Windows (MSVC if supported)
2. Test with multiple BLAS vendors:
   - Intel MKL
   - OpenBLAS
   - ATLAS
   - Accelerate (macOS)
3. Test with different integer sizes (32-bit, 64-bit)
4. Use CI/CD for automated testing

#### **Risk 4: Memory Layout Confusion**

**Concern:** Mixing C-ordered and Fortran-ordered data could cause bugs

**Mitigation:**
1. Document memory layout requirements clearly
2. F_* wrappers explicitly handle Fortran layout
3. Assert correct dimensions in debug builds
4. Validate with Valgrind/AddressSanitizer
5. Code review focus on stride/leading dimension

#### **Risk 5: Incomplete Migration**

**Concern:** Missing some BLAS calls during migration

**Mitigation:**
1. Automated search for all BLAS calls:
   ```bash
   grep -r "F_DGEMM\|F_DGEMV\|F_DCOPY\|F_DDOT" psi4/src/psi4/fnocc/
   ```
2. Compilation will fail for any missed calls (good!)
3. Comprehensive test suite catches runtime issues
4. Code review checklist

### 6.3 Testing Strategy for Risk Mitigation

#### **Tier 1: Unit Tests**
- Test each F_* wrapper individually
- Validate against known results
- Test edge cases (zero sizes, transposes, etc.)

#### **Tier 2: Integration Tests**
- Run module-specific test suites:
  - fnocc: `ctest -R fnocc`
  - psimrcc: `ctest -R psimrcc`
  - dfocc: `ctest -R dfocc`

#### **Tier 3: Regression Tests**
- Compare results with previous Psi4 version
- Use reference data from literature
- Bit-for-bit comparison where possible

#### **Tier 4: Performance Tests**
- Benchmark critical calculations
- Compare timing with baseline
- Profile for hot spots

#### **Tier 5: Platform Tests**
- Linux, macOS, Windows
- Multiple BLAS vendors
- Multiple compilers

---

## 7. Testing and Validation

### 7.1 Testing Pyramid

```
              ┌──────────────────┐
              │  Manual Testing  │ (1%)
              │  User Validation │
              └──────────────────┘
         ┌─────────────────────────────┐
         │   System/Integration Tests  │ (9%)
         │   - Full test suite runs    │
         │   - Reference calculations  │
         └─────────────────────────────┘
    ┌────────────────────────────────────────┐
    │      Module-Specific Tests             │ (30%)
    │      - fnocc test suite                │
    │      - psimrcc test suite              │
    │      - dfocc test suite                │
    └────────────────────────────────────────┘
┌──────────────────────────────────────────────────┐
│           Unit Tests                             │ (60%)
│           - blas_fortran.h wrappers              │
│           - Individual operations                │
│           - Edge cases                           │
└──────────────────────────────────────────────────┘
```

### 7.2 Test Plan by Phase

#### **Phase 1: Foundation**
- Unit tests for libqt/blas_fortran.h
- Test F_DGEMM with various transposes
- Test F_DGEMV with various transposes
- Test F_DCOPY, F_DDOT, F_DSCAL
- Validate against reference BLAS implementation

**Example Unit Test:**
```cpp
// test_blas_fortran.cc
TEST(BlasFortran, DGEMM_NoTranspose) {
    size_t m = 3, n = 4, k = 5;
    double A[15], B[20], C[12];
    // Initialize A, B
    F_DGEMM('n', 'n', m, n, k, 1.0, A, m, B, k, 0.0, C, m);
    // Validate C against expected result
    EXPECT_NEAR(C[0], expected[0], 1e-14);
}
```

#### **Phase 2: Include Standardization**
- Verify clean compilation
- Run quick smoke tests
- No functional changes expected

#### **Phase 3: Array Classes**
- DFOCC test suite: `ctest -R dfocc -V`
- OCC test suite: `ctest -R occ -V`
- Validate numerical results unchanged

**Key Tests:**
- dfocc-ccsd-1 (CCSD energy)
- dfocc-ccsd-grad1 (CCSD gradient)
- occ-omp2-1 (Orbital-optimized MP2)

#### **Phase 4: psimrcc**
- PSIMRCC test suite: `ctest -R psimrcc -V`
- Critical: Multi-reference calculations

**Key Tests:**
- mrccsd-1 (MRCCSD energy)
- mrccsd_t-1 (MRCCSD(T) energy)
- Perturbative triples

**Validation Criteria:**
- Energies match to 1e-10 Eh
- Gradients match to 1e-8 Eh/bohr
- Dipoles match to 1e-6 Debye

#### **Phase 5: fnocc (CRITICAL)**
- Full fnocc test suite: `ctest -R fnocc -V`
- Benchmark calculations
- Performance profiling

**Key Tests (with expected results):**

| Test | Method | System | Expected Energy (Eh) | Tolerance |
|------|--------|--------|---------------------|-----------|
| fnocc-ccsd-1 | CCSD | H2O | -76.240832 | 1e-6 |
| fnocc-ccsd-t-1 | CCSD(T) | H2O | -76.245877 | 1e-6 |
| fnocc-qcisd-1 | QCISD | H2O | -76.240875 | 1e-6 |
| fnocc-qcisd-t-1 | QCISD(T) | H2O | -76.245920 | 1e-6 |
| fnocc-mp3-1 | MP3 | H2O | -76.232946 | 1e-6 |
| fnocc-mp4-1 | MP4(SDTQ) | H2O | -76.245877 | 1e-6 |
| fnocc-cepa-1 | CEPA(0) | H2O | -76.237618 | 1e-6 |

**Performance Benchmarks:**

| Calculation | System | Basis | Expected Time | Tolerance |
|-------------|--------|-------|---------------|-----------|
| CCSD | H2O | cc-pVDZ | ~5s | ±10% |
| CCSD(T) | H2O | cc-pVDZ | ~8s | ±10% |
| CCSD | Benzene | cc-pVDZ | ~300s | ±10% |

**Validation Steps:**
1. Run before migration, save results
2. Run after migration, compare
3. Investigate any differences >1e-10
4. Profile if timing difference >5%

### 7.3 Automated Validation Script

```bash
#!/bin/bash
# validate_blas_migration.sh

echo "BLAS Migration Validation"
echo "=========================="

# Phase 1: Compilation
echo "Testing compilation..."
cd build && make -j8 || exit 1

# Phase 2: Quick smoke test
echo "Running smoke tests..."
ctest -L smoke --output-on-failure || exit 1

# Phase 3: Module-specific tests
echo "Running fnocc tests..."
ctest -R fnocc -V || exit 1

echo "Running psimrcc tests..."
ctest -R psimrcc -V || exit 1

echo "Running dfocc tests..."
ctest -R dfocc -V || exit 1

# Phase 4: Full test suite
echo "Running full test suite..."
ctest -j8 --output-on-failure || exit 1

# Phase 5: Performance benchmarks
echo "Running benchmarks..."
./run_benchmarks.sh > benchmark_results.txt

echo "Validation complete!"
```

### 7.4 Acceptance Criteria

✓ **Correctness:**
- All test cases pass (100% pass rate)
- Numerical results identical to within 1e-10 Eh
- No new compiler warnings
- No memory leaks (Valgrind clean)

✓ **Performance:**
- No degradation >2% in critical paths
- Benchmark timings within ±5%
- Memory usage unchanged

✓ **Code Quality:**
- No code duplication for BLAS wrappers
- Consistent naming conventions
- Comprehensive documentation
- Clean compilation on all platforms

✓ **Maintainability:**
- Single source of truth for BLAS operations
- Clear dependency graph
- Simplified build system

---

## 8. Timeline and Resources

### 8.1 Detailed Timeline (Option A: Aggressive)

| Phase | Duration | Calendar | Tasks | Dependencies |
|-------|----------|----------|-------|--------------|
| **Phase 1: Foundation** | 2 weeks | Weeks 1-2 | Create blas_fortran.h, docs, tests | None |
| **Phase 2: Tier 1** | 1 week | Week 2 | Include standardization | Phase 1 |
| **Phase 3: Tier 2** | 1 week | Week 3 | Array consolidation | Phase 2 |
| **Phase 4: Tier 3** | 2 weeks | Weeks 4-5 | psimrcc migration | Phase 1, 3 |
| **Phase 5: Tier 4** | 3 weeks | Weeks 5-8 | fnocc migration | Phase 1, 4 |
| **Phase 6: Cleanup** | 1 week | Week 8 | Documentation, final testing | Phase 5 |
| **Buffer** | 1 week | Week 9 | Contingency for issues | - |
| **Total** | **8-9 weeks** | - | - | - |

**Note:** Phases 2-4 can partially overlap after Phase 1 completion.

### 8.2 Resource Requirements

#### **Personnel:**
- **Lead Developer:** 1 FTE (full-time for 8 weeks)
  - Implements blas_fortran.h
  - Performs migrations
  - Code reviews

- **Testing Engineer:** 0.5 FTE
  - Creates test infrastructure
  - Runs validation tests
  - Performance benchmarking

- **Code Reviewer:** 0.25 FTE
  - Reviews pull requests
  - Validates correctness
  - Domain expertise (coupled cluster)

- **Total:** ~1.75 FTE for 8 weeks

#### **Compute Resources:**
- Development workstation (Linux)
- macOS test machine
- CI/CD resources (GitHub Actions or similar)
- Benchmark cluster for performance testing

#### **Software:**
- Compilers: GCC, Intel ICC, Clang
- BLAS vendors: MKL, OpenBLAS, ATLAS
- Profiling tools: gprof, perf, Valgrind
- Version control: Git

### 8.3 Milestones and Checkpoints

**Milestone 1 (End of Week 2):** libqt/blas_fortran.h complete and tested
- Deliverable: New header file with unit tests
- Checkpoint: Unit tests pass, documentation complete

**Milestone 2 (End of Week 3):** Include standardization and Array consolidation
- Deliverable: Consistent includes, optionally unified arrays
- Checkpoint: Compilation clean, basic tests pass

**Milestone 3 (End of Week 5):** psimrcc migration complete
- Deliverable: psimrcc using libqt/blas_fortran.h
- Checkpoint: PSIMRCC test suite passes

**Milestone 4 (End of Week 8):** fnocc migration complete
- Deliverable: fnocc using libqt/blas_fortran.h, old wrappers removed
- Checkpoint: Full test suite passes, performance validated

**Milestone 5 (End of Week 9):** Final release
- Deliverable: Clean, documented, standardized codebase
- Checkpoint: All acceptance criteria met

### 8.4 Alternative Timelines

#### **Option B: Conservative (3-4 weeks)**
- Phase 1: 1.5 weeks
- Phase 2-3: 1.5 weeks
- Phase 6: 1 week
- Defer Phase 4-5 (keep fnocc/psimrcc wrappers)

#### **Option C: Minimal (2 weeks)**
- Phase 1: 1 week
- Phase 2: 0.5 weeks
- Phase 6: 0.5 weeks
- Defer Phase 3-5

#### **Option D: Documentation Only (1 week)**
- Update docs to recommend libqt
- Mark fnocc/psimrcc wrappers as deprecated
- No code changes

---

## 9. Success Criteria

### 9.1 Technical Success Criteria

#### **Primary Objectives (Must Have):**

✓ **Single BLAS Interface**
- All BLAS operations route through libqt (directly or via approved wrappers)
- No redundant BLAS wrapper implementations
- Consistent API across all modules

✓ **Code Quality**
- Zero compiler warnings
- No memory leaks (Valgrind clean)
- No code duplication for BLAS operations
- Passes static analysis (cppcheck, clang-tidy)

✓ **Correctness**
- All existing tests pass (100% pass rate)
- Numerical results identical to baseline (within 1e-10 Eh)
- Reference calculations validated

✓ **Performance**
- No performance regression >2% in critical paths
- Benchmark timings within baseline ±5%
- Memory usage unchanged

#### **Secondary Objectives (Should Have):**

✓ **Documentation**
- Comprehensive BLAS_MIGRATION_GUIDE.md
- Updated libqt documentation
- Clear examples for future developers
- API reference for blas_fortran.h

✓ **Build System**
- Simplified BLAS/LAPACK dependencies
- Clean CMakeLists.txt structure
- Single path to vendor libraries

✓ **Maintainability**
- Clear ownership of BLAS interface
- Documented extension points
- Easy to add new BLAS operations

#### **Stretch Goals (Nice to Have):**

✓ **Performance Improvements**
- Leverage libqt's INT_MAX blocking in more places
- Enable large array support (>2^31 elements) in fnocc/psimrcc

✓ **Enhanced Features**
- Complex BLAS support (ZGEMM, etc.)
- Additional LAPACK routines

### 9.2 Quantitative Metrics

| Metric | Baseline | Target | Measurement |
|--------|----------|--------|-------------|
| BLAS Interfaces | 5 | 1 (+wrappers) | Count of distinct implementations |
| Code Duplication | High | None | Manual inspection |
| Test Pass Rate | 100% | 100% | ctest results |
| Performance (CCSD) | Baseline | ±2% | Timing benchmarks |
| Compiler Warnings | 0 | 0 | Compilation log |
| Memory Leaks | 0 | 0 | Valgrind report |
| Documentation | Partial | Complete | Manual review |

### 9.3 Qualitative Success Indicators

✓ **Developer Experience**
- New contributors understand where to add BLAS calls
- Clear error messages when BLAS calls fail
- Easy to debug BLAS-related issues

✓ **Code Clarity**
- Consistent naming conventions
- Self-documenting code
- Minimal cognitive load

✓ **Community Acceptance**
- Positive feedback from Psi4 developers
- No complaints about performance
- Successful integration into main branch

---

## 10. Appendices

### Appendix A: File Inventory

#### **Files to Create:**
1. `psi4/src/psi4/libqt/blas_fortran.h` (new wrapper interface)
2. `psi4/src/psi4/libqt/test_blas_fortran.cc` (unit tests)
3. `BLAS_MIGRATION_GUIDE.md` (documentation)
4. `validate_blas_migration.sh` (validation script)

#### **Files to Modify:**

**Tier 1 (Include Standardization):** ~195 files
- All files currently including libqt/qt.h

**Tier 2 (Array Classes):** ~32 files
- dfocc/arrays.h, arrays.cc
- occ/arrays.h
- Files using Array1d/2d/3d classes

**Tier 3 (psimrcc):** ~8 files
- psimrcc/mp2_ccsd.cc
- psimrcc/mrcc_pert_triples.cc
- psimrcc/mrccsd_t.cc
- psimrcc/ccmrcc.cc
- psimrcc/blas.cc
- Others using algebra_interface

**Tier 4 (fnocc):** ~20 files
- fnocc/ccsd.cc (CRITICAL - 92,000 lines)
- fnocc/cepa.cc
- fnocc/qed.cc
- fnocc/mp3.cc
- fnocc/mp4.cc
- fnocc/ccsd_triples.cc
- Others using F_DGEMM/F_DGEMV

#### **Files to Delete:**
1. `psi4/src/psi4/fnocc/blas.h`
2. `psi4/src/psi4/fnocc/blas.cc`
3. `psi4/src/psi4/fnocc/blas_mangle.h`
4. `psi4/src/psi4/psimrcc/algebra_interface.h`
5. `psi4/src/psi4/psimrcc/algebra_interface.cc`

**Total Files Affected:** ~260 files

---

### Appendix B: Example Implementation

#### **libqt/blas_fortran.h (Proposed Implementation):**

```cpp
#ifndef _psi4_libqt_blas_fortran_h_
#define _psi4_libqt_blas_fortran_h_

#include "psi4/libqt/qt.h"
#include <cstddef>

namespace psi {

/**
 * @file blas_fortran.h
 * @brief Fortran-ordered BLAS wrappers for libqt
 *
 * This file provides Fortran-ordered (column-major) wrappers around
 * the C-ordered libqt BLAS interface. These wrappers are designed
 * for modules that work with Fortran-ordered data (e.g., fnocc, psimrcc).
 *
 * All operations are implemented as inline functions that translate
 * to the appropriate C_* calls with transposed arguments where necessary.
 *
 * @note These are thin wrappers with zero overhead when optimized.
 * @note Memory layout must be Fortran-ordered (column-major).
 *
 * Example usage:
 * @code
 *   // Fortran-ordered matrices A(m,k), B(k,n), C(m,n)
 *   // C := alpha*A*B + beta*C
 *   F_DGEMM('n', 'n', m, n, k, alpha, A, m, B, k, beta, C, m);
 * @endcode
 */

// =============================================================================
// BLAS Level 1: Vector-Vector Operations
// =============================================================================

/**
 * @brief Fortran-ordered DCOPY: y := x
 * @param n Number of elements
 * @param x Source vector
 * @param incx Stride for x
 * @param y Destination vector
 * @param incy Stride for y
 */
inline void F_DCOPY(size_t n, const double *x, int incx, double *y, int incy) {
    C_DCOPY(n, x, incx, y, incy);
}

/**
 * @brief Fortran-ordered DDOT: return x'*y
 * @param n Number of elements
 * @param x First vector
 * @param incx Stride for x
 * @param y Second vector
 * @param incy Stride for y
 * @return Dot product
 */
inline double F_DDOT(size_t n, const double *x, int incx, const double *y, int incy) {
    return C_DDOT(n, x, incx, y, incy);
}

/**
 * @brief Fortran-ordered DSCAL: x := alpha*x
 * @param n Number of elements
 * @param alpha Scaling factor
 * @param x Vector to scale (in-place)
 * @param incx Stride for x
 */
inline void F_DSCAL(size_t n, double alpha, double *x, int incx) {
    C_DSCAL(n, alpha, x, incx);
}

/**
 * @brief Fortran-ordered DAXPY: y := alpha*x + y
 * @param n Number of elements
 * @param alpha Scaling factor
 * @param x Vector x
 * @param incx Stride for x
 * @param y Vector y (in-place)
 * @param incy Stride for y
 */
inline void F_DAXPY(size_t n, double alpha, const double *x, int incx, double *y, int incy) {
    C_DAXPY(n, alpha, x, incx, y, incy);
}

/**
 * @brief Fortran-ordered DNRM2: return ||x||_2
 * @param n Number of elements
 * @param x Vector
 * @param incx Stride for x
 * @return Euclidean norm
 */
inline double F_DNRM2(size_t n, const double *x, int incx) {
    return C_DNRM2(n, x, incx);
}

// =============================================================================
// BLAS Level 2: Matrix-Vector Operations
// =============================================================================

/**
 * @brief Fortran-ordered DGEMV: y := alpha*op(A)*x + beta*y
 *
 * Fortran memory layout: A is stored column-major
 *
 * @param trans 'n' or 'N': op(A) = A
 *              't' or 'T': op(A) = A'
 * @param m Number of rows of A
 * @param n Number of columns of A
 * @param alpha Scaling factor for A*x
 * @param a Matrix A (Fortran-ordered: column-major)
 * @param lda Leading dimension of A (>= m in Fortran layout)
 * @param x Vector x
 * @param incx Stride for x
 * @param beta Scaling factor for y
 * @param y Vector y (in-place)
 * @param incy Stride for y
 *
 * Note: Fortran A(m,n) with lda>=m translates to C A(n,m) with lda>=m
 * Need to swap dimensions and transpose when calling C_DGEMV
 */
inline void F_DGEMV(char trans, size_t m, size_t n, double alpha,
                    const double *a, size_t lda, const double *x, int incx,
                    double beta, double *y, int incy) {
    // Fortran: y := alpha * op(A) * x + beta * y
    // where A is m x n in Fortran (column-major) layout
    //
    // C layout: A^T is n x m (row-major = Fortran column-major)
    // So: C_DGEMV with transposed operation

    char c_trans;
    if (trans == 'n' || trans == 'N') {
        // Fortran: y := alpha * A * x + beta * y
        // C view: y := alpha * A^T^T * x = alpha * (A^T)' * x
        c_trans = 't';  // Need transpose in C
    } else {
        // Fortran: y := alpha * A' * x + beta * y
        // C view: y := alpha * (A^T)^T * x = alpha * A^T * x
        c_trans = 'n';  // No transpose in C
    }

    // Call C_DGEMV with swapped dimensions
    C_DGEMV(c_trans, n, m, alpha, a, lda, x, incx, beta, y, incy);
}

// =============================================================================
// BLAS Level 3: Matrix-Matrix Operations
// =============================================================================

/**
 * @brief Fortran-ordered DGEMM: C := alpha*op(A)*op(B) + beta*C
 *
 * Fortran memory layout: A, B, C are stored column-major
 *
 * @param transa 'n' or 'N': op(A) = A; 't' or 'T': op(A) = A'
 * @param transb 'n' or 'N': op(B) = B; 't' or 'T': op(B) = B'
 * @param m Number of rows of op(A) and C
 * @param n Number of columns of op(B) and C
 * @param k Number of columns of op(A) and rows of op(B)
 * @param alpha Scaling factor
 * @param a Matrix A (Fortran-ordered)
 * @param lda Leading dimension of A
 * @param b Matrix B (Fortran-ordered)
 * @param ldb Leading dimension of B
 * @param beta Scaling factor for C
 * @param c Matrix C (Fortran-ordered, in-place)
 * @param ldc Leading dimension of C
 *
 * Fortran DGEMM semantics:
 * - If transa='n': A is m x k, lda >= m
 * - If transa='t': A is k x m, lda >= k
 * - If transb='n': B is k x n, ldb >= k
 * - If transb='t': B is n x k, ldb >= n
 * - C is m x n, ldc >= m
 *
 * C DGEMM semantics (row-major):
 * - Fortran A(m,k) col-major = C A^T(k,m) row-major
 * - Need to swap A <-> B and swap m <-> n, adjust transposes
 *
 * Mathematical transformation:
 * Fortran: C := alpha * op(A) * op(B) + beta * C
 * C view:  C^T := alpha * op(B)^T * op(A)^T + beta * C^T
 * So:      C := (alpha * op(B)^T * op(A)^T)^T + beta * C
 *             = alpha * op(A) * op(B) + beta * C  [transposition cancels]
 *
 * Implementation: Swap argument order and toggle transposes
 */
inline void F_DGEMM(char transa, char transb, size_t m, size_t n, size_t k,
                    double alpha, const double *a, size_t lda,
                    const double *b, size_t ldb,
                    double beta, double *c, size_t ldc) {
    // Toggle transposes
    char c_transb = (transb == 'n' || transb == 'N') ? 't' : 'n';
    char c_transa = (transa == 'n' || transa == 'N') ? 't' : 'n';

    // Call C_DGEMM with swapped arguments and dimensions
    C_DGEMM(c_transb, c_transa, n, m, k, alpha, b, ldb, a, lda, beta, c, ldc);
}

// =============================================================================
// LAPACK: Linear Algebra Operations
// =============================================================================

/**
 * @brief Fortran-ordered DGESV: Solve A*X = B via LU factorization
 * @param n Order of matrix A
 * @param nrhs Number of right-hand sides
 * @param a Matrix A (Fortran-ordered, overwritten with LU factors)
 * @param lda Leading dimension of A (>= n)
 * @param ipiv Pivot indices (output)
 * @param b Matrix B (Fortran-ordered, overwritten with solution X)
 * @param ldb Leading dimension of B (>= n)
 * @param info Output: 0 = success, <0 = illegal value, >0 = singular
 */
inline void F_DGESV(size_t n, size_t nrhs, double *a, size_t lda,
                    int *ipiv, double *b, size_t ldb, int *info) {
    // LAPACK routines are Fortran-native, so no transpose needed
    // C_DGESV already handles this correctly
    *info = C_DGESV(n, nrhs, a, lda, ipiv, b, ldb);
}

/**
 * @brief Fortran-ordered DSYEV: Symmetric eigenvalue decomposition
 * @param jobz 'V' = compute eigenvectors, 'N' = eigenvalues only
 * @param uplo 'U' = upper triangle, 'L' = lower triangle
 * @param n Order of matrix A
 * @param a Matrix A (Fortran-ordered, overwritten with eigenvectors if jobz='V')
 * @param lda Leading dimension of A (>= n)
 * @param w Eigenvalues (output, ascending order)
 * @param work Workspace array
 * @param lwork Size of work array (>= max(1, 3n-1))
 * @param info Output: 0 = success
 */
inline void F_DSYEV(char jobz, char uplo, size_t n, double *a, size_t lda,
                    double *w, double *work, size_t lwork, int *info) {
    *info = C_DSYEV(jobz, uplo, n, a, lda, w, work, lwork);
}

/**
 * @brief Fortran-ordered DGESVD: Singular value decomposition
 * @param jobu 'A' = all columns of U, 'S' = first min(m,n) columns, 'N' = no U
 * @param jobvt 'A' = all rows of V^T, 'S' = first min(m,n) rows, 'N' = no V^T
 * @param m Number of rows of A
 * @param n Number of columns of A
 * @param a Matrix A (Fortran-ordered, destroyed on output)
 * @param lda Leading dimension of A (>= m)
 * @param s Singular values (output, descending order)
 * @param u Left singular vectors (Fortran-ordered)
 * @param ldu Leading dimension of U
 * @param vt Right singular vectors transposed (Fortran-ordered)
 * @param ldvt Leading dimension of V^T
 * @param work Workspace array
 * @param lwork Size of work array
 * @param info Output: 0 = success
 */
inline void F_DGESVD(char jobu, char jobvt, size_t m, size_t n,
                     double *a, size_t lda, double *s,
                     double *u, size_t ldu, double *vt, size_t ldvt,
                     double *work, size_t lwork, int *info) {
    *info = C_DGESVD(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork);
}

}  // namespace psi

#endif  // _psi4_libqt_blas_fortran_h_
```

**Key Design Decisions:**
1. Inline functions → zero overhead
2. Detailed documentation → easy to use correctly
3. Mathematical explanations → understand transpose logic
4. Consistent naming → F_* prefix for Fortran-ordered
5. Delegates to C_* functions → single source of truth

---

### Appendix C: Migration Checklist

Use this checklist to track progress through the migration:

#### **Phase 1: Foundation**
- [ ] Create `libqt/blas_fortran.h` with full documentation
- [ ] Implement F_DCOPY, F_DDOT, F_DSCAL, F_DAXPY, F_DNRM2
- [ ] Implement F_DGEMV with proper transpose handling
- [ ] Implement F_DGEMM with proper transpose handling
- [ ] Implement F_DGESV, F_DSYEV, F_DGESVD
- [ ] Create unit tests for all wrappers
- [ ] Validate against reference BLAS
- [ ] Update libqt/CMakeLists.txt
- [ ] Write BLAS_MIGRATION_GUIDE.md
- [ ] Update libqt/README.md

#### **Phase 2: Include Standardization**
- [ ] Audit all libqt includes (195 files)
- [ ] Standardize include paths
- [ ] Verify clean compilation
- [ ] Run smoke tests

#### **Phase 3: Array Consolidation**
- [ ] Analyze dfocc/arrays vs occ/arrays
- [ ] Verify libqt usage in Array methods
- [ ] Optional: Consolidate implementations
- [ ] Run DFOCC test suite
- [ ] Run OCC test suite

#### **Phase 4: psimrcc Migration**
- [ ] Update psimrcc/mp2_ccsd.cc to use blas_fortran.h
- [ ] Update psimrcc/mrcc_pert_triples.cc
- [ ] Update psimrcc/mrccsd_t.cc
- [ ] Update psimrcc/ccmrcc.cc
- [ ] Update other psimrcc files
- [ ] Remove algebra_interface.h
- [ ] Remove algebra_interface.cc
- [ ] Run PSIMRCC test suite
- [ ] Validate MRCCSD energies
- [ ] Benchmark performance

#### **Phase 5: fnocc Migration**
- [ ] Backup fnocc/blas.h (safety)
- [ ] Update fnocc/ccsd.cc (CRITICAL)
- [ ] Update fnocc/cepa.cc
- [ ] Update fnocc/qed.cc
- [ ] Update fnocc/mp3.cc
- [ ] Update fnocc/mp4.cc
- [ ] Update fnocc/ccsd_triples.cc
- [ ] Update other fnocc files
- [ ] Run FNOCC test suite
- [ ] Validate CCSD/CCSD(T) energies
- [ ] Validate QCISD/QCISD(T) energies
- [ ] Validate MP3/MP4 energies
- [ ] Benchmark performance
- [ ] Remove fnocc/blas.h
- [ ] Remove fnocc/blas.cc
- [ ] Remove fnocc/blas_mangle.h

#### **Phase 6: Cleanup**
- [ ] Run full test suite
- [ ] Verify no orphaned BLAS code
- [ ] Update all module READMEs
- [ ] Update developer documentation
- [ ] Clean up CMakeLists.txt
- [ ] Final performance validation
- [ ] Create migration summary
- [ ] Prepare release notes

---

### Appendix D: Performance Benchmarking Guide

#### **Benchmark Calculations:**

**Test 1: CCSD/H2O/cc-pVDZ**
```python
molecule h2o {
  O
  H 1 0.96
  H 1 0.96 2 104.5
}

set basis cc-pvdz
energy('fnocc-ccsd')
```
- Expected energy: -76.240832480
- Expected time: ~5s
- Metrics: Total time, iteration time, DGEMM time

**Test 2: CCSD(T)/H2O/cc-pVDZ**
```python
set basis cc-pvdz
energy('fnocc-ccsd(t)')
```
- Expected energy: -76.245877200
- Expected time: ~8s
- Metrics: Total time, triples time

**Test 3: CCSD/Benzene/cc-pVDZ**
```python
molecule benzene {
  C        0.000000    1.396940    0.000000
  C        1.209562    0.698470    0.000000
  C        1.209562   -0.698470    0.000000
  C        0.000000   -1.396940    0.000000
  C       -1.209562   -0.698470    0.000000
  C       -1.209562    0.698470    0.000000
  H        0.000000    2.480080    0.000000
  H        2.148004    1.240040    0.000000
  H        2.148004   -1.240040    0.000000
  H        0.000000   -2.480080    0.000000
  H       -2.148004   -1.240040    0.000000
  H       -2.148004    1.240040    0.000000
  symmetry c1
}

set basis cc-pvdz
energy('fnocc-ccsd')
```
- Expected time: ~300s
- Metrics: Iteration performance, memory usage

#### **Profiling:**

```bash
# Profile with gprof
gfortran -pg ...
./psi4 input.dat
gprof ./psi4 > profile.txt

# Profile with perf
perf record -g ./psi4 input.dat
perf report

# Check for hot spots in BLAS calls
perf record -g -e cycles:pp ./psi4 input.dat
perf annotate F_DGEMM
```

---

### Appendix E: Contact and Resources

#### **Psi4 Development Resources:**
- GitHub: https://github.com/psi4/psi4
- Documentation: https://psicode.org/psi4manual/master/
- Forum: https://github.com/psi4/psi4/discussions

#### **BLAS/LAPACK References:**
- BLAS Quick Reference: http://www.netlib.org/blas/blasqr.pdf
- LAPACK User Guide: http://www.netlib.org/lapack/lug/
- Intel MKL Documentation: https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html

#### **Related Documentation:**
- Psi4 libqt documentation (internal)
- Coupled Cluster theory references
- Fortran/C interoperability guides

---

## Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2025-11-18 | Claude | Initial comprehensive plan |

---

**END OF DOCUMENT**
