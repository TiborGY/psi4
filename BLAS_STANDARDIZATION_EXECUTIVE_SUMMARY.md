# BLAS Standardization Executive Summary

**Status:** Planning Complete - Ready for Implementation
**Priority:** Medium-High
**Estimated Effort:** 6-8 weeks
**Risk:** Low-Medium

---

## The Situation

Psi4 currently has **5 different BLAS wrapper implementations** doing essentially the same thing:

1. **libqt** (88% of codebase) - Universal C-ordered interface ✅
2. **fnocc/blas.h** (8%) - Fortran-ordered thin wrappers
3. **psimrcc/algebra_interface.h** (1.5%) - Fortran wrappers around libqt
4. **dfocc/arrays.h** (1.9%) - Array class wrappers
5. **occ/arrays.h** (0.1%) - Array class wrappers

**Good news:** The codebase is already 88% standardized on libqt!

**The problem:** Fragmentation creates maintenance burden, code duplication, and missed opportunities.

---

## The Solution

**Establish libqt as the single standard BLAS interface** while preserving existing high-level abstractions.

### Proposed Architecture

```
┌──────────────────────────────┐
│   Application Modules        │
│  (fnocc, psimrcc, dfocc...)  │
└──────────────────────────────┘
               ↓
┌──────────────────────────────┐
│  Optional Convenience Layer  │
│  • libqt/blas_fortran.h      │ ← NEW: Thin Fortran wrappers
│  • CCBLAS (keep)             │
│  • Array classes (keep)      │
└──────────────────────────────┘
               ↓
┌──────────────────────────────┐
│    libqt Core Interface      │ ← Single source of truth
│  (already 88% standardized)  │
└──────────────────────────────┘
               ↓
┌──────────────────────────────┐
│  Vendor BLAS/LAPACK          │
│  (MKL, OpenBLAS, ATLAS...)   │
└──────────────────────────────┘
```

---

## What Gets Created

### New File: `libqt/blas_fortran.h`

A thin Fortran-ordered wrapper layer providing:
- `F_DGEMM()` - Matrix-matrix multiply (Fortran semantics)
- `F_DGEMV()` - Matrix-vector multiply
- `F_DCOPY()`, `F_DDOT()`, `F_DSCAL()`, etc.
- `F_DGESV()`, `F_DSYEV()`, `F_DGESVD()` - LAPACK operations

**Key features:**
- Zero overhead (inline functions)
- Maps to libqt C_* functions internally
- Handles Fortran/C layout translation automatically
- Same API as current fnocc/blas.h

---

## What Gets Removed

1. `fnocc/blas.h` → Use `libqt/blas_fortran.h` instead
2. `fnocc/blas.cc` → No longer needed
3. `fnocc/blas_mangle.h` → libqt handles this
4. `psimrcc/algebra_interface.h` → Use `libqt/blas_fortran.h` instead
5. `psimrcc/algebra_interface.cc` → No longer needed

**What stays:**
- psimrcc CCBLAS (high-level tensor framework)
- dfocc/occ Array classes (convenience layer)
- All algorithm code (no changes to ccsd.cc logic!)

---

## Implementation Strategy

### 4-Tier Phased Approach

#### **Tier 1: Include Standardization (1 week, LOW RISK)**
- Ensure consistent `#include "psi4/libqt/qt.h"` across 195 files
- Already using libqt, just standardize format
- **Impact:** Clean, consistent includes

#### **Tier 2: Array Consolidation (1 week, LOW RISK)**
- Verify dfocc/occ arrays use libqt internally (already do)
- Optional: Consolidate dfocc and occ array implementations
- **Impact:** No code duplication in array classes

#### **Tier 3: psimrcc Migration (2 weeks, MEDIUM RISK)**
- Create `libqt/blas_fortran.h`
- Migrate 67 calls in 8 psimrcc files
- Remove `algebra_interface.h`
- **Impact:** psimrcc standardized on libqt

#### **Tier 4: fnocc Migration (3 weeks, MEDIUM-HIGH RISK)**
- Migrate 355 calls in ~20 fnocc files
- **Critical:** ccsd.cc (92,000 lines with 150+ F_DGEMM calls)
- Remove `fnocc/blas.h`
- **Impact:** Complete standardization

### Timeline Options

| Option | Scope | Timeline | Risk | Standardization |
|--------|-------|----------|------|-----------------|
| **A - Aggressive (Recommended)** | All 4 tiers | 6-8 weeks | Medium | 100% |
| **B - Conservative** | Tiers 1-3 only | 3-4 weeks | Low | 92% |
| **C - Minimal** | Tiers 1-2 only | 2 weeks | Very Low | 90% |
| **D - Doc Only** | Documentation | 1 week | None | 88% (status quo) |

---

## Benefits

### Immediate
✅ **Single source of truth** for BLAS operations
✅ **Eliminate code duplication** (5 implementations → 1 core + wrappers)
✅ **Simplified maintenance** (one place to fix bugs)
✅ **Consistent behavior** across all modules

### Long-term
✅ **Enable large array support** (>2^31 elements) everywhere via libqt's INT_MAX blocking
✅ **Easier platform porting** (centralized name mangling)
✅ **Better debugging** (single BLAS interface to profile/trace)
✅ **Clearer codebase** for new contributors

---

## Risks and Mitigation

### Risk 1: Performance Regression in CCSD
- **Likelihood:** Low
- **Impact:** High
- **Mitigation:**
  - Make all wrappers `inline` (zero overhead)
  - Benchmark before/after
  - Acceptance: <2% variation
  - Rollback plan: Keep fnocc/blas.h on branch

### Risk 2: Numerical Differences
- **Likelihood:** Very Low
- **Impact:** High
- **Mitigation:**
  - F_DGEMM wrapper exactly replicates Fortran semantics
  - Extensive validation (all tests must pass)
  - Bit-for-bit comparison with baseline

### Risk 3: Large Code Change in ccsd.cc
- **Likelihood:** Medium
- **Impact:** Medium
- **Mitigation:**
  - API stays identical (just include change)
  - Extensive testing with full fnocc suite
  - Phase 5 is separate, can defer if needed

---

## Success Criteria

### Must Have (Phase Completion)
✓ All tests pass (100% pass rate)
✓ Numerical results identical to baseline (within 1e-10 Eh)
✓ No performance regression >2%
✓ Zero compiler warnings
✓ Single BLAS interface (libqt core)

### Should Have (Quality)
✓ Comprehensive documentation
✓ Clean build system
✓ Simplified dependencies

### Nice to Have (Stretch Goals)
✓ Performance improvements from libqt features
✓ Enhanced BLAS/LAPACK coverage

---

## Recommended Next Steps

### Week 1-2: Foundation (Start Here!)

**1. Create libqt/blas_fortran.h**
```cpp
// Implement thin Fortran-ordered wrappers
inline void F_DGEMM(...) { /* transpose and call C_DGEMM */ }
inline void F_DGEMV(...) { /* transpose and call C_DGEMV */ }
// etc.
```

**2. Unit Tests**
```cpp
// Test wrapper correctness
TEST(BlasFortran, DGEMM_NoTranspose) { ... }
TEST(BlasFortran, DGEMM_TransposeA) { ... }
```

**3. Documentation**
- Create BLAS_MIGRATION_GUIDE.md
- Update libqt/README.md
- Document blas_fortran.h API

**Deliverable:** Working libqt/blas_fortran.h with tests

---

### Week 2-3: Low-Risk Migrations

**4. Tier 1: Include Standardization**
- Ensure consistent libqt includes
- Verify compilation

**5. Tier 2: Array Verification**
- Confirm arrays use libqt internally
- Run DFOCC/OCC tests

**Deliverable:** Clean, standardized includes

---

### Week 4-5: psimrcc Migration

**6. Migrate psimrcc**
```cpp
// Before
#include "algebra_interface.h"

// After
#include "psi4/libqt/blas_fortran.h"
```

**7. Testing**
- Run PSIMRCC test suite
- Validate MRCCSD calculations
- Benchmark performance

**Deliverable:** psimrcc standardized on libqt

---

### Week 6-8: fnocc Migration (Critical Phase)

**8. Migrate fnocc**
```cpp
// Before
#include "blas.h"

// After
#include "psi4/libqt/blas_fortran.h"
```

**9. Extensive Testing**
- Full fnocc test suite
- CCSD/CCSD(T) validation
- Performance benchmarks

**10. Cleanup**
- Remove fnocc/blas.h
- Remove psimrcc/algebra_interface.h
- Final documentation

**Deliverable:** Complete standardization

---

## Decision Points

### Decision 1: Which Option to Pursue?

**Recommendation: Option A (Aggressive - All 4 Tiers)**

**Rationale:**
- Benefits outweigh risks
- 88% already standardized (momentum)
- Proper testing mitigates risks
- Complete solution vs partial

**Alternative:** Start with Option B (Tiers 1-3), defer Tier 4 if issues arise

---

### Decision 2: Array Class Consolidation?

**Options:**
- **A:** Keep dfocc and occ arrays separate (less work)
- **B:** Consolidate to single implementation (cleaner)

**Recommendation:** Option A (keep separate) unless significant duplication found

---

### Decision 3: Timeline?

**Options:**
- **Aggressive:** 6-8 weeks (1 FTE)
- **Conservative:** 12-16 weeks (0.5 FTE, part-time)
- **Phased:** Complete Tier 1-2 immediately, defer 3-4 to later release

**Recommendation:** Aggressive if resources available, otherwise phased

---

## Resource Requirements

**Personnel:**
- 1 FTE developer (8 weeks)
- 0.5 FTE testing engineer (concurrent)
- 0.25 FTE code reviewer (concurrent)

**Compute:**
- Development workstation (Linux)
- Test machines (Linux, macOS)
- CI/CD resources

**Timeline:**
- **Optimistic:** 6 weeks
- **Realistic:** 8 weeks
- **Pessimistic:** 10 weeks (with buffer)

---

## Quick Start Commands

### 1. Create the new wrapper file
```bash
cd psi4/src/psi4/libqt
# Create blas_fortran.h (see full plan for implementation)
```

### 2. Add unit tests
```bash
# Create test_blas_fortran.cc
# Run: ctest -R test_blas_fortran
```

### 3. Migrate psimrcc (example)
```bash
cd psi4/src/psi4/psimrcc
# Replace #include "algebra_interface.h" with:
# #include "psi4/libqt/blas_fortran.h"
```

### 4. Test
```bash
cd build
ctest -R psimrcc -V
```

### 5. Migrate fnocc (when ready)
```bash
cd psi4/src/psi4/fnocc
# Replace #include "blas.h" with:
# #include "psi4/libqt/blas_fortran.h"
```

### 6. Validate
```bash
ctest -R fnocc -V
# Run benchmarks
```

---

## Documentation Provided

1. **LIBQT_BLAS_STANDARDIZATION_PLAN.md** (this branch)
   - Comprehensive 70-page plan
   - Detailed implementation guide
   - Risk assessment and mitigation
   - Testing strategy
   - Complete example code

2. **BLAS_STANDARDIZATION_EXECUTIVE_SUMMARY.md** (this file)
   - Quick overview
   - Decision framework
   - Immediate next steps

3. **Previous Analysis** (from exploration)
   - BLAS_USAGE_MAP.txt
   - DETAILED_MODULE_BREAKDOWN.txt
   - MIGRATION_ACTION_ITEMS.txt

---

## Questions to Answer

**Before starting:**

1. **Do we have resources for 6-8 weeks?**
   - If yes → Proceed with Option A (Aggressive)
   - If no → Consider Option B (Conservative) or Option C (Minimal)

2. **What's the timeline priority?**
   - If urgent → Option C (2 weeks, 90% standardization)
   - If normal → Option A (8 weeks, complete)

3. **Who owns this effort?**
   - Assign lead developer
   - Assign code reviewer with coupled cluster expertise
   - Assign testing resources

4. **What's the merge strategy?**
   - Single PR or multiple PRs per tier?
   - Feature branch or direct to main?
   - Release coordination

**During implementation:**

5. **Are tests passing at each tier?**
   - Gate progress on test results
   - Don't proceed to next tier if failures

6. **Is performance acceptable?**
   - Benchmark at each tier
   - Investigate if >2% variation

7. **Are there unexpected issues?**
   - Adjust plan as needed
   - Consider rollback if major problems

---

## Contact and Next Steps

**Ready to proceed?**

1. Review this summary and the comprehensive plan
2. Make decision on which option to pursue (A/B/C/D)
3. Allocate resources (personnel, time, compute)
4. Set up development branch
5. Begin Phase 1: Create libqt/blas_fortran.h

**Questions or concerns?**
- Review full plan: LIBQT_BLAS_STANDARDIZATION_PLAN.md
- Check implementation examples in Appendix B
- Review risk mitigation in Section 6

**Approval to start:**
- [ ] Option selected (A/B/C/D)
- [ ] Resources allocated
- [ ] Timeline approved
- [ ] Development branch created
- [ ] Ready to begin Phase 1

---

**Bottom line:** This is a **feasible, low-risk standardization effort** that will **significantly improve code maintainability** with **minimal disruption** to existing functionality. The phased approach allows us to **validate at each step** and **rollback if needed**.

**Recommendation: Proceed with Option A (Aggressive).** Start with Phase 1 (Foundation) and validate thoroughly before moving to subsequent phases.

---

**Document Status:** Ready for Review and Approval
**Next Action:** Decision on which option to pursue
**Expected Start:** Upon approval
