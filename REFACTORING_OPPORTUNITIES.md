# Code Duplication Analysis: Refactoring Opportunities in PSI4

**Analysis Date:** 2025-11-18
**Scope:** psi4/src/psi4 C++ codebase
**Total Files Analyzed:** 800+

## Executive Summary

This document identifies significant code duplication across the PSI4 quantum chemistry package. After comprehensive analysis of the C++ source code, we've identified **20,000-40,000 lines of duplicated code** across multiple categories. The duplications range from identical copy-paste code to similar algorithmic patterns that could share a common implementation.

### Impact Assessment
- **Maintenance burden:** Changes must be replicated across multiple files
- **Bug propagation:** Fixes may not be applied uniformly
- **Code quality:** Inconsistent implementations of same algorithms
- **Testing overhead:** Similar code requires duplicate test coverage
- **Learning curve:** New contributors must understand multiple implementations of same concepts

---

## Critical Priority Duplications

### 1. BLAS/LAPACK Wrapper Duplication (CRITICAL)
**Impact:** Code maintenance nightmare

**Files with identical content:**
- `psi4/src/psi4/fnocc/blas.h` (lines 29-165)
- `psi4/src/psi4/libqt/blas_intfc23_mangle.h` (lines 29-164)
- `psi4/src/psi4/psimrcc/algebra_interface_mangle.h` (lines 29-117)
- `psi4/src/psi4/mcscf/algebra_interface_mangle.h` (lines 29-117)

**What's duplicated:**
- Identical F_DGEMM, F_DGEMV, F_DSYEV, F_DGESVD macros
- Identical FC_SYMBOL handling logic (#if FC_SYMBOL == 2/1/3/4)
- Identical extern "C" declarations for BLAS/LAPACK routines

**Recommendation:**
→ **Create single header:** `psi4/src/psi4/libblas/blas_interface.h`
→ **Migration effort:** Low (header replacement only)
→ **Risk:** Low (compile-time verification)

---

### 2. MOInfo Structure Duplication (CRITICAL)
**Impact:** ~7 identical structure definitions

**Files:**
- `psi4/src/psi4/cc/ccdensity/MOInfo.h`
- `psi4/src/psi4/cc/ccenergy/MOInfo.h`
- `psi4/src/psi4/cc/cceom/MOInfo.h`
- `psi4/src/psi4/cc/cclambda/MOInfo.h`
- `psi4/src/psi4/cc/cchbar/MOInfo.h`
- `psi4/src/psi4/cc/ccresponse/MOInfo.h`
- `psi4/src/psi4/cc/cctriples/MOInfo.h`

**Already exists:** `psi4/src/psi4/libmoinfo/moinfo.h` (unused by CC modules!)

**What's duplicated (in all 7 files):**
- nirreps, nmo, nso, nao
- orbspi, clsdpi, openpi, uoccpi, frdocc, fruocc
- occpi, aoccpi, boccpi, virtpi, avirtpi, bvirtpi
- occ_sym, aocc_sym, bocc_sym, vir_sym, avir_sym, bvir_sym
- occ_off, vir_off arrays
- enuc, escf, eref, ecc fields

**Recommendation:**
→ **Migrate all CC modules to:** `libmoinfo::MOInfo`
→ **Migration effort:** Medium (requires updating get_moinfo.cc files)
→ **Risk:** Low (data structure replacement)

---

### 3. Tensor/Array Class Duplication (CRITICAL)
**Impact:** Could eliminate 50-100 files of duplication

**Duplicate implementations:**
- `psi4/src/psi4/dfocc/tensors.cc` - Tensor1d/2d/3d classes
- `psi4/src/psi4/dfocc/arrays.cc` - Array1d/2d/3d classes
- `psi4/src/psi4/occ/arrays.cc` - Array1d/2d/3d classes

**What's duplicated (identical logic):**
- Memory allocation/deallocation (memalloc, release, init)
- Zero initialization (zero() method using memset)
- Printing functions (print() methods)
- Vector operations (set, add, subtract, get)
- Mathematical operations (rms, dot, gemv)
- Matrix operations (to_lower_triangle, vector_dot)

**Size:**
- DFOCC: 123 files, 108,000 lines
- OCC: 54 files, similar structure
- Estimated duplication: 5,000-10,000 lines

**Recommendation:**
→ **Create:** `psi4/src/psi4/libtensor/` unified tensor library
→ **Use templates** for different dimensionalities and spin cases
→ **Migration effort:** High (extensive refactoring)
→ **Risk:** Medium (requires thorough testing)

---

### 4. Multiple CPHF/CPKS Solver Implementations (CRITICAL)
**Impact:** 4 completely independent response solvers

**Implementations:**
1. `psi4/src/psi4/libfock/apps.h:141` - `class RCPHF : public RBase`
2. `psi4/src/psi4/libfock/hamiltonian.h:142` - `class CPHFRHamiltonian`
3. `psi4/src/psi4/fisapt/fisapt.h:213-267` - `class CPHF_FISAPT`
4. `psi4/src/psi4/libsapt_solver/usapt0.h:338` - `class CPKS_USAPT0`

**What's duplicated:**
- CPHF/CPKS iterative solution procedure
- Preconditioning strategies
- Product formation (Hamiltonian-vector products)
- Convergence criteria

**Recommendation:**
→ **Consolidate to:** Single CPHF solver in libfock
→ **Provide:** Specialized interfaces for FISAPT/SAPT
→ **Migration effort:** High
→ **Risk:** Medium

---

### 5. Multiple MP2 Implementations (CRITICAL)
**Impact:** 9+ separate MP2 energy implementations

**Files:**
- `psi4/src/psi4/cc/ccenergy/mp2_energy.cc`
- `psi4/src/psi4/dfocc/cc_energy.cc`
- `psi4/src/psi4/occ/cc_energy.cc`
- `psi4/src/psi4/dfmp2/mp2.cc`
- `psi4/src/psi4/fnocc/mp2.cc`
- `psi4/src/psi4/f12/mp2.cc`
- `psi4/src/psi4/dlpno/mp2.cc`
- `psi4/src/psi4/dct/dct_mp2_RHF.cc`
- `psi4/src/psi4/dct/dct_mp2_UHF.cc`

**What's duplicated:**
- Core MP2 energy expression: E_MP2 = Σ t_ij^ab * <ij||ab>
- Amplitude formation: t_ij^ab = <ij|ab> / (ε_i + ε_j - ε_a - ε_b)
- Energy accumulation and spin scaling

**Recommendation:**
→ **Create base class:** `MP2Base` with core algorithm
→ **Specialize for:** DF-MP2, Local-MP2, F12-MP2, etc.
→ **Migration effort:** High
→ **Risk:** Medium

---

### 6. Multiple Triples (T) Implementations (CRITICAL)
**Impact:** 9+ separate (T) correction implementations

**Files:**
- `psi4/src/psi4/cc/cctriples/triples.cc`
- `psi4/src/psi4/fnocc/triples.cc`
- `psi4/src/psi4/fnocc/lowmemory_triples.cc`
- `psi4/src/psi4/dfocc/ccsd_triples.cc`
- `psi4/src/psi4/dfocc/uccsd_triples_hm.cc`
- `psi4/src/psi4/dct/dct_triples.cc`
- `psi4/src/psi4/psimrcc/mrcc_pert_triples.cc`
- Plus separate spin-case files (AAA, AAB, ABB, BBB)

**What's duplicated:**
- T3 amplitude formation via contractions
- (T) energy expression evaluation
- Permutation operators

**Recommendation:**
→ **Unify around:** Single (T) kernel with integral backend abstraction
→ **Template by:** Spin case and memory strategy
→ **Migration effort:** Very High
→ **Risk:** Medium-High

---

### 7. Multiple DIIS Implementations (CRITICAL)
**Impact:** 8+ DIIS implementations when libdiis exists

**Files:**
1. `psi4/src/psi4/libdiis/diismanager.{cc,h}` - **General implementation (underutilized!)**
2. `psi4/src/psi4/cc/ccenergy/diis_RHF.cc`
3. `psi4/src/psi4/cc/ccenergy/diis_ROHF.cc`
4. `psi4/src/psi4/cc/ccenergy/diis_UHF.cc`
5. `psi4/src/psi4/cc/cclambda/diis.cc`
6. `psi4/src/psi4/cc/ccresponse/diis.cc`
7. `psi4/src/psi4/dfocc/diis.cc`
8. `psi4/src/psi4/fnocc/diis.cc`
9. `psi4/src/psi4/psimrcc/blas_diis.cc`
10. `psi4/src/psi4/mcscf/scf_diis.cc`

**What's duplicated:**
- Error vector construction
- B-matrix formation and solving
- Extrapolation
- Storage management

**Recommendation:**
→ **Mandate use of:** `libdiis/DIISManager`
→ **Enhance if needed:** Add missing features to DIISManager
→ **Remove:** All duplicate implementations
→ **Migration effort:** Medium
→ **Risk:** Low-Medium

---

## High Priority Duplications

### 8. Integral Sorting and Transformation
**Files:** ~2,000-3,000 lines of duplication

- `psi4/src/psi4/dct/dct_integrals_RHF.cc`
- `psi4/src/psi4/dct/dct_integrals_UHF.cc`
- `psi4/src/psi4/fnocc/sortintegrals.cc` (2,198 lines!)
- `psi4/src/psi4/psimrcc/sort_out_of_core.cc`
- `psi4/src/psi4/psimrcc/sort_mrpt2.cc`

**Duplication:** Separate implementations for OOOO, OOVV, OVOV, OVVV, VVVV integrals
**Recommendation:** Create unified `IntegralSorter` class with template methods

---

### 9. get_moinfo() Functions
**Files:** 9 nearly identical implementations

All CC modules have `get_moinfo.cc` with identical initialization logic
**Recommendation:** Create base `get_moinfo()` in `libccbase/`

---

### 10. Davidson Iterative Diagonalization
**Files:** 3+ independent implementations

- `psi4/src/psi4/detci/diag_h.cc` (lines 143-278)
- `psi4/src/psi4/detci/mitrush_iter.cc`
- `psi4/src/psi4/cc/cceom/diag.cc`

**Recommendation:** Create unified `DavidsonSolver` class in libqt

---

### 11. MO Integral Transformation
**Files:** 4+ separate transformation engines

- `psi4/src/psi4/libtrans/integraltransform.h` - **Main implementation**
- `psi4/src/psi4/dct/half_transform.cc` - DCT custom
- `psi4/src/psi4/psimrcc/transform.h` - PSIMRCC custom
- `psi4/src/psi4/fnocc/df_t1_transformation.cc` - FNOCC custom

**Recommendation:** Standardize on `libtrans::IntegralTransform`

---

### 12. Orbital Orthogonalization
**Files:** 6+ implementations of same algorithms

- `psi4/src/psi4/libmints/orthog.{cc,h}` - Symmetric, Canonical, Cholesky
- `psi4/src/psi4/libqt/schmidt.cc` - Gram-Schmidt
- `psi4/src/psi4/mcscf/scf_S_inverse_sqrt.cc` - Overlap inversion
- `psi4/src/psi4/libmints/matrix.cc` - schmidt_orthog_columns(), canonical_orthogonalization()
- `psi4/src/psi4/detci/civect.cc` - Custom Gram-Schmidt
- `psi4/src/psi4/libscf_solver/rohf.cc` - prepare_canonical_orthogonalization()

**Recommendation:** Consolidate to `libmints/orthog` methods

---

### 13. DCT Module RHF/UHF Duplication
**Files:** 22 paired files doing nearly identical operations

Pattern: Nearly every DCT file has `_RHF.cc` and `_UHF.cc` versions:
- dct_compute_{RHF,UHF}.cc
- dct_scf_{RHF,UHF}.cc
- dct_oo_{RHF,UHF}.cc
- dct_integrals_{RHF,UHF}.cc
- dct_intermediates_{RHF,UHF}.cc
- dct_energy_{RHF,UHF}.cc
- dct_gradient_{RHF,UHF}.cc
- dct_lambda_{RHF,UHF}.cc
- dct_tau_{RHF,UHF}.cc
- dct_density_{RHF,UHF}.cc
- dct_mp2_{RHF,UHF}.cc

**Recommendation:** Template or policy-based design for spin cases

---

### 14. Amplitude Update Code (T1/T2)
**Files:** 38 files, ~3,000-5,000 lines duplication

Multiple implementations in:
- cc/ccenergy, cc/cclambda, cc/cceom
- dfocc: separate files for each method (ccd, ccsd, lccd, etc.)
- occ: similar duplication
- psimrcc: separate update routines

**Recommendation:** Template amplitude update by excitation level

---

### 15. Natural Orbital Formation
**Files:** 16 files reference natural orbitals

- `psi4/src/psi4/fnocc/frozen_natural_orbitals.cc` - Full implementation
- `psi4/src/psi4/libscf_solver/sad.cc` - "Natural Orbitals via UHF"
- Referenced in dfocc, dfmp2, occ modules

**Recommendation:** Centralize in libmints or libtrans

---

## Medium Priority Duplications

### 16. SCF Orbital Construction (form_C)
**Files:** libscf_solver/{rhf,uhf,rohf,cuhf}.cc

Near-identical level shifting and Fock diagonalization
**Recommendation:** Extract common methods to HF base class

---

### 17. Fock Matrix Building
**Files:** 58 files with fock_build patterns

Duplicate J/K matrix construction logic across modules
**Recommendation:** Consolidate with libfock JK infrastructure

---

### 18. Density Matrix Formation
**Files:** dct/dct_intermediates_{RHF,UHF}.cc, libtrans/integraltransform_tpdm.cc

Similar density formation from orbitals
**Recommendation:** Unified density builder with spin-case templates

---

### 19. Diagonalization Wrappers
**Files:** 500-800 lines across multiple modules

Multiple classes implement their own wrappers to LAPACK
**Recommendation:** Standardize on Matrix::diagonalize()

---

### 20. Convergence Checking Patterns
**Files:** 25+ implementations, ~1,250-2,500 lines

Near-identical while loops checking convergence
**Recommendation:** Create ConvergenceChecker utility class

---

### 21. Iteration Printing and Diagnostics
**Files:** 50+ iterative methods, ~1,500-2,500 lines

Nearly identical iteration header printing and diagnostics
**Recommendation:** Create IterationLogger utility class

---

### 22. Exception Handling
**Files:** 217+ files, 1,341 PSIEXCEPTION throws

Inconsistent error messages for similar errors
**Recommendation:** Validation utilities with standard error messages

---

### 23. Wavefunction Initialization
**Files:** 107+ files, ~5,000-15,000 lines estimated

Every wavefunction class has nearly identical initialization boilerplate
**Recommendation:** CRTP or mixins for common initialization

---

### 24. Memory Allocation Patterns
**Files:** fnocc extensively uses malloc/free

MemoryManager exists but not consistently used
**Recommendation:** Mandate MemoryManager or smart pointers

---

### 25. Derivative/Gradient Computation
**Files:** 1,500-2,500 lines duplication

- scfgrad/response.cc, scfgrad/jk_grad.cc
- dfmp2/corr_grad.cc
- libmints/mintshelper.cc

Duplicate derivative integral setup
**Recommendation:** Template derivative computation by integral type

---

### 26. Schwarz Screening
**Files:** Multiple implementations, ~200-400 lines

Different screening approaches in libmints, libfock, lib3index
**Recommendation:** Centralize screening logic

---

### 27. Grid Generation and DFT Integration
**Files:** ~400-600 lines duplication

libfock/cubature.cc, libfock/v.cc, libmints/numinthelper.cc
**Recommendation:** Single grid generation factory

---

### 28. JK Builder Implementations
**Files:** 8+ JK classes, ~500-1,000 lines duplication

Some initialization patterns duplicated across JK builders
**Recommendation:** Extract more common code to JK base class

---

### 29. String Manipulation Functions
**Files:** libpsi4util/stl_string.cc vs local implementations

- Central implementation: `to_lower()`, `to_upper()`, `split()`
- Duplicated in: libmints/basisset.cc, liboptions/liboptions.cc

**Recommendation:** All modules use libpsi4util implementation

---

### 30. Timing/Profiling Code
**Files:** 215+ files manually manage timer_on()/timer_off()

**Recommendation:** RAII-based ScopedTimer class

---

### 31. PSIO Operations
**Files:** 168 files, manual file I/O patterns

Repeated psio_open/close/read/write with manual error checking
**Recommendation:** RAII-based PSIOFile class

---

### 32. get_params() Functions
**Files:** 15+ cc subdirectories

Similar parameter reading across CC modules
**Recommendation:** Base parameter reader class

---

### 33. print_header() Functions
**Files:** 30+ implementations

Similar banner/header formatting across modules
**Recommendation:** HeaderPrinter utility class

---

### 34. Three-Index/DF Operations
**Files:** 167 files reference 3-index operations

B-matrix formation duplicated across modules
**Recommendation:** Consolidate to lib3index/dfhelper

---

### 35. SAPT vs FISAPT Code
**Files:** Major duplication between SAPT and FISAPT

Both compute similar interaction energy components
**Recommendation:** Review shared code opportunities

---

## Summary Statistics

### Code Duplication by Category

| Category | Estimated Lines | Priority | Affected Modules |
|----------|----------------|----------|------------------|
| Tensor/Array classes | 5,000-10,000 | Critical | dfocc, occ |
| Wavefunction init | 5,000-15,000 | Medium | All |
| CC energy computation | 5,000-8,000 | High | dfocc, occ, cc |
| Amplitude updates | 3,000-5,000 | High | cc variants |
| Integral sorting | 2,000-3,000 | High | dct, fnocc, psimrcc |
| Derivative computation | 1,500-2,500 | Medium | scfgrad, dfmp2 |
| Convergence checking | 1,250-2,500 | Medium | All iterative |
| DIIS implementations | 1,500-2,000 | Critical | Multiple |
| Iteration diagnostics | 1,500-2,500 | Medium | All |

### Top Priority Actions

1. **Immediate:**
   - Consolidate BLAS/LAPACK headers (4 files → 1)
   - Migrate CC modules to libmoinfo (eliminate 7 duplicate headers)
   - Mandate use of libdiis/DIISManager

2. **Short-term:**
   - Consolidate CPHF/CPKS solvers
   - Unify MP2 implementations around base class
   - Create IntegralSorter class
   - Extract common get_moinfo() to libccbase

3. **Long-term:**
   - Refactor tensor/array library
   - Template DCT for spin cases
   - Unify triples implementations
   - Davidson solver consolidation

### Estimated Impact

**Total Duplicated Code:** 20,000-40,000 lines
**Potential Reduction:** 15,000-30,000 lines (50-75% reduction)
**Modules Most Affected:**
- DFOCC/OCC (highest duplication)
- CC modules (structural duplication)
- DCT (RHF/UHF splits)
- Response theory (scattered implementations)

---

## Refactoring Strategies

### 1. Template Metaprogramming
Use templates for spin-case variants (RHF/UHF/ROHF)

**Example:**
```cpp
template<SpinType ST>
class DCTEnergy {
    void compute_energy();  // Single implementation for all spin cases
};
```

### 2. Policy-Based Design
For algorithm variants

**Example:**
```cpp
template<typename IntegralBackend, typename MemoryPolicy>
class TriplesCorrection {
    // Single (T) implementation, multiple backends
};
```

### 3. Facade Pattern
For existing underutilized libraries

**Example:**
```cpp
// Instead of reimplementing DIIS:
DIISManager diis;  // Use existing libdiis
```

### 4. CRTP/Mixins
For common wavefunction boilerplate

**Example:**
```cpp
template<typename Derived>
class WavefunctionInit {
    void common_init() { /* shared initialization */ }
};
```

### 5. Factory Patterns
For creating similar objects

**Example:**
```cpp
class IntegralSorterFactory {
    static unique_ptr<IntegralSorter> create(SortType type);
};
```

---

## Recommendations for New Code

1. **Before implementing new feature:**
   - Search for existing implementations
   - Check libqt, libmints, libpsi4util for utilities

2. **Prefer composition over duplication:**
   - Use existing libraries (DIISManager, IntegralTransform, etc.)
   - Extend through inheritance, not copy-paste

3. **Template for spin cases:**
   - Don't create separate RHF/UHF files
   - Template or use polymorphism

4. **Follow DRY principle:**
   - Extract common patterns to utilities
   - Document why code is duplicated if unavoidable

5. **Use RAII:**
   - For timers, file I/O, memory management
   - Reduces boilerplate dramatically

---

## Migration Guide

### Phase 1: Low-Hanging Fruit (Weeks 1-4)
- Consolidate BLAS headers → `libblas/blas_interface.h`
- Migrate CC modules to libmoinfo
- Create HeaderPrinter utility
- Create IterationLogger utility
- Create ConvergenceChecker utility

### Phase 2: Structural Improvements (Months 2-3)
- Mandate DIISManager usage
- Create IntegralSorter class
- Consolidate get_moinfo() functions
- Standardize orthogonalization methods
- Create RAII wrappers (Timer, PSIOFile)

### Phase 3: Major Refactoring (Months 4-6)
- Consolidate CPHF/CPKS solvers
- Unify MP2 implementations
- Create common tensor library design
- Davidson solver unification

### Phase 4: Algorithm Unification (Months 7-12)
- Template DCT for spin cases
- Unify triples implementations
- Consolidate MO transformation
- Refactor amplitude update code

---

## Testing Strategy

For each refactoring:

1. **Preserve existing tests:** All current test cases must pass
2. **Add regression tests:** Ensure numerical results unchanged
3. **Performance testing:** Verify no performance regression
4. **Code coverage:** Ensure refactored code is well-tested
5. **Documentation:** Update documentation for API changes

---

## Conclusion

The PSI4 codebase contains significant duplication resulting from 30 years of development by dozens of contributors. While some duplication may be intentional (e.g., specialized implementations for performance), much of it appears to be inadvertent.

**Key Findings:**
- 20,000-40,000 lines of duplicated code identified
- Critical duplications in core infrastructure (BLAS wrappers, MOInfo, DIIS)
- Significant algorithmic duplication (MP2, triples, CPHF)
- Structural duplication (RHF/UHF splits, wavefunction init)

**Benefits of Refactoring:**
- Reduced maintenance burden
- Improved code quality and consistency
- Easier onboarding for new developers
- Fewer bugs through consolidation
- Better test coverage

**Recommendation:**
Proceed with phased refactoring starting with critical, low-risk consolidations (BLAS headers, MOInfo structures, DIIS usage) and gradually tackling larger structural improvements.
