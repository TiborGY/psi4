# DIIS Implementation Consolidation Analysis

## Executive Summary

This analysis investigates the feasibility of consolidating DIIS (Direct Inversion in the Iterative Subspace) implementations across the Psi4 codebase. The investigation reveals that **most modules are already using the centralized `libdiis/DIISManager`**, contrary to the initial assumption of widespread fragmentation.

### Key Findings:

1. **occ, dfocc, and dct modules already use libdiis/DIISManager** - No migration needed
2. **psimrcc has a legitimate custom implementation** due to tight integration with CCBLAS framework
3. **~112 lines of dead code identified** in dfocc/diis.cc that can be removed
4. **Estimated consolidation impact: ~110 lines** (much less than the originally estimated 500-1000 lines)

---

## Module-by-Module Analysis

### 1. libdiis/DIISManager (Central Implementation)

**Location:** `psi4/src/psi4/libdiis/`

**Architecture:** Python-backed C++ wrapper
- C++ front-end: `diismanager.h` and `diismanager.cc` (73 lines)
- Python back-end: `psi4/driver/procrouting/diis.py` (437 lines)

**Capabilities:**
- Storage policies: InCore, OnDisk
- Removal policies: LargestError, OldestAdded
- DIIS engines: DIIS, ADIIS, EDIIS (with hybrid capabilities)
- Supported data types: Matrix, Vector, dpdbuf4, dpdfile2, BlockedTensor, float
- Automatic caching of dot products for efficiency
- Variadic templates for flexible multi-vector extrapolation

**Key Insight:** This is a robust, feature-rich implementation that handles all common DIIS use cases.

---

### 2. occ Module

**Status:** ✅ **Already using libdiis/DIISManager**

**DIIS Functions:**
- `oo_diis()` - orbital optimization DIIS (psi4/src/psi4/occ/occ_iterations.cc:501)
- `cepa_diis()` - CEPA amplitude DIIS (psi4/src/psi4/occ/cepa_iterations.cc:242)

**Implementation Details:**
```cpp
DIISManager orbital_diis(maxdiis_, "Orbital Optimized DIIS",
                         DIISManager::RemovalPolicy::LargestError,
                         DIISManager::StoragePolicy::OnDisk);
```

**Usage Pattern:**
- Works with DPD buffers (dpdbuf4) for T2 amplitudes
- Uses Vector objects for orbital gradients and rotations
- Handles both RHF and UHF references
- Supports multiple wavefunction types (OMP2, OMP3, OCEPA, OREMP)

**Verdict:** No changes needed. Implementation is clean and efficient.

---

### 3. dfocc Module

**Status:** ✅ **Already using libdiis/DIISManager** + 🔴 **Dead code to remove**

**DIIS Managers:**
- `orbitalDIIS` - orbital optimization (psi4/src/psi4/dfocc/dfocc.h:656)
- `ccsdDiisManager` - CCSD amplitude DIIS (psi4/src/psi4/dfocc/dfocc.h:654)

**Implementation Details:**
```cpp
orbitalDIIS = std::make_shared<DIISManager>(
    cc_maxdiis_, "Orbital Optimized DIIS",
    DIISManager::RemovalPolicy::LargestError,
    DIISManager::StoragePolicy::OnDisk);
```

**Files Using DIISManager:**
- `occ_iterations.cc` - orbital optimization iterations
- `ccsd_iterations.cc` - CCSD iterations
- `ccd_iterations.cc` - CCD iterations
- `lccd_iterations.cc` - LCCD iterations
- `remp_iterations.cc` - REMP iterations
- Plus corresponding `*_t2_amps.cc` files

**Dead Code Identified:**
- `psi4/src/psi4/dfocc/diis.cc` (112 lines) - Custom DIIS function no longer called
- Function signature: `void DFOCC::diis(int dimvec, SharedTensor2d &vecs, ...)`
- **Evidence of obsolescence:**
  - No references found in any dfocc source files
  - Comment in `update_mo.cc:60` states: "replaced by coupled DIIS in oo_diis() above"
  - All DIIS operations now use orbitalDIIS or ccsdDiisManager

**Recommendation:** Remove `diis.cc` and corresponding declaration in `dfocc.h:148-149`.

---

### 4. dct Module

**Status:** ✅ **Already using libdiis/DIISManager**

**DIIS Usage Locations:**
- `dct_oo_RHF.cc:66` - Orbital optimization
- `dct_oo_UHF.cc:61` - Orbital optimization (UHF)
- `dct_compute_RHF.cc:131` - SCF iterations
- `dct_compute_UHF.cc:178, 267, 344` - Multiple DIIS managers for amplitudes, orbitals
- `dct_tau_UHF.cc:722` - Tau DIIS
- `dct_gradient_UHF.cc:109, 1118, 1913` - Gradient response iterations
- `dct_qc.cc:73` - Quadratically convergent iterations

**Implementation Pattern:**
```cpp
DIISManager diisManager(maxdiis_, "DCT DIIS vectors");
diisManager.set_error_vector_size(orbital_gradient_a_.get(), orbital_gradient_b_.get(), &Laa, &Lab, &Lbb);
diisManager.set_vector_size(Xtotal_a_.get(), Xtotal_b_.get(), &Laa, &Lab, &Lbb);
```

**Usage Pattern:**
- Extensive use across different DCT variants
- Mixes Matrix/Vector objects with DPD buffers
- Uses both InCore and OnDisk storage policies depending on context
- Conditional DIIS based on convergence thresholds

**Verdict:** No changes needed. Comprehensive and appropriate usage.

---

### 5. psimrcc Module

**Status:** 🟡 **Custom implementation - migration not recommended**

**Location:** `psi4/src/psi4/psimrcc/blas_diis.cc` (207 lines)

**Custom Implementation Details:**

**Data Structure:**
```cpp
std::vector<std::pair<std::string, std::string>> diis_matrices;  // (amps, delta_amps)
```

**Key Functions:**
- `diis_add()` - Add matrix pairs to DIIS tracking
- `diis_save_t_amps()` - Save amplitudes to PSIO
- `diis()` - Main DIIS extrapolation

**Unique Features:**
1. **Tight CCBLAS Integration:**
   - Works with CCBLAS matrix infrastructure (CCMatIrTmp)
   - Uses string-based matrix naming and lookup
   - Integrated with PSIMRCC's multi-reference framework

2. **Custom Storage:**
   - Uses PSIF_PSIMRCC_INTEGRALS file
   - Custom labeling scheme: `{name}_DIIS_{h}_{step}`
   - Direct F_DDOT calls for dot products

3. **Specialized Logic:**
   - Two DIIS modes: DiisEachCycle, DiisCC
   - Handles singularities explicitly
   - Custom cycle-based indexing

**Why Migration is Difficult:**

1. **Architectural Mismatch:**
   - libdiis expects specific Psi4 object types (Matrix, Vector, dpdbuf4)
   - psimrcc uses CCBLAS matrix infrastructure with string-based lookup
   - Adapting would require significant refactoring of CCBLAS framework

2. **Feature Requirements:**
   - Cycle-based vector management (modulo indexing)
   - String-based matrix identification
   - Integration with existing PSIMRCC PSIO infrastructure

3. **Code Complexity:**
   - PSIMRCC is a complex multi-reference module
   - DIIS is deeply integrated with the iteration logic
   - Refactoring risk is high for marginal benefit

**Recommendation:** Keep custom implementation. The 207 lines are justified by:
- Architectural necessity
- Module-specific requirements
- Risk vs. reward of migration

---

## Consolidation Opportunities

### Immediate Actions (Low Risk, High Value)

#### 1. Remove Dead Code in dfocc
**Files to remove:**
- `psi4/src/psi4/dfocc/diis.cc` (112 lines)

**Changes needed:**
- Remove from `psi4/src/psi4/dfocc/CMakeLists.txt`
- Remove declaration from `psi4/src/psi4/dfocc/dfocc.h:148-149`:
  ```cpp
  void diis(int dimvec, SharedTensor2d &vecs, SharedTensor2d &errvecs,
            SharedTensor1d &vec_new, SharedTensor1d &errvec_new);
  ```

**Impact:** Clean up 112 lines of unused code, reduce maintenance burden

**Risk:** Very low - code is confirmed unused

---

### Future Considerations (Optional Enhancements)

#### 1. Create Shared DIIS Helper Functions for occ/dfocc

While both modules already use libdiis, there is some pattern repetition in:
- Setting up error vector and vector sizes
- Converting between Tensor and Matrix/Vector formats
- DIIS extrapolation patterns

**Potential Shared Utilities:**
```cpp
namespace psi {
namespace occutils {
    // Helper to set up DIIS for standard T2 amplitudes
    void setup_t2_diis(DIISManager& mgr, int naoccA, int navirA, ...);

    // Helper to add T2 amplitudes to DIIS
    void add_t2_to_diis(DIISManager& mgr, dpdbuf4* T, dpdbuf4* R, ...);
}
}
```

**Estimated savings:** ~50-100 lines across modules
**Risk:** Medium - refactoring working code
**Priority:** Low - current code is functional

---

## Detailed Metrics

### Code Volume by Component

| Component | Lines of Code | Status | Action |
|-----------|--------------|--------|--------|
| libdiis (C++) | 73 | Active | Keep |
| libdiis (Python) | 437 | Active | Keep |
| occ DIIS usage | ~150 | Active, using libdiis | Keep |
| dfocc DIIS usage | ~200 | Active, using libdiis | Keep |
| dfocc/diis.cc | 112 | **DEAD CODE** | **Remove** |
| dct DIIS usage | ~250 | Active, using libdiis | Keep |
| psimrcc DIIS | 207 | Active, custom | Keep |

### Consolidation Impact

**Original estimate:** 500-1000 lines of redundant code
**Actual redundancy:** 112 lines (dfocc/diis.cc only)

**Reduction in estimate:** The original assumption that occ, dfocc, and dct had redundant DIIS implementations was incorrect. These modules are already consolidated around libdiis.

---

## Recommendations Summary

### High Priority (Implement Now)

1. ✅ **Remove dfocc/diis.cc dead code**
   - Impact: -112 lines
   - Risk: Very low
   - Effort: 30 minutes

### Low Priority (Consider for Future)

2. 🔶 **Create shared DIIS helper utilities** (optional)
   - Impact: ~50-100 line reduction
   - Risk: Medium
   - Effort: 2-4 hours
   - Benefit: Modest - current code works fine

3. 🔶 **Document libdiis best practices** (optional)
   - Create examples showing common DIIS patterns
   - Add inline documentation for typical use cases
   - Would help future developers

### Not Recommended

4. ❌ **Migrate psimrcc to libdiis**
   - Impact: Would require major refactoring
   - Risk: High
   - Benefit: Minimal - current implementation is appropriate

---

## Technical Details

### libdiis/DIISManager API Summary

**Initialization:**
```cpp
DIISManager mgr(max_vecs, "Label", RemovalPolicy, StoragePolicy);
```

**Setup:**
```cpp
mgr.set_error_vector_size(error1, error2, ...);  // Variadic
mgr.set_vector_size(vec1, vec2, ...);            // Variadic
```

**Iteration:**
```cpp
bool added = mgr.add_entry(error1, error2, ..., vec1, vec2, ...);
if (mgr.subspace_size() >= min_vecs) {
    mgr.extrapolate(vec1, vec2, ...);
}
```

**Cleanup:**
```cpp
mgr.delete_diis_file();  // For OnDisk storage
```

### Supported Data Types

The Python backend (via pybind11) supports:
- `core.Matrix` - 2D matrices
- `core.Vector` - 1D vectors
- `core.dpdbuf4` - 4-index DPD buffers
- `core.dpdfile2` - 2-index DPD files
- `ambit.BlockedTensor` - Ambit tensors
- `float` - Scalar values

---

## Conclusion

The DIIS implementation landscape in Psi4 is **much better than initially assumed**:

1. **Three major modules (occ, dfocc, dct) are already using the centralized libdiis**
2. **Only one module (psimrcc) has a justified custom implementation** due to architectural requirements
3. **One piece of dead code (dfocc/diis.cc)** can be removed for immediate cleanup
4. **No major consolidation work is needed** - the codebase is already well-organized

**Overall Assessment:** The original concern about DIIS fragmentation was overstated. The situation is largely resolved, with only minor cleanup needed.

---

## Implementation Plan

### Phase 1: Immediate Cleanup (Recommended)

**Task:** Remove dead code from dfocc

**Steps:**
1. Remove `psi4/src/psi4/dfocc/diis.cc`
2. Remove function declaration from `psi4/src/psi4/dfocc/dfocc.h` (lines 148-149)
3. Update `psi4/src/psi4/dfocc/CMakeLists.txt` to remove diis.cc from sources
4. Build and test dfocc module
5. Run dfocc regression tests to confirm no breakage

**Estimated time:** 30-60 minutes
**Risk:** Very low
**Tests to run:** All dfocc tests in test suite

### Phase 2: Optional Enhancements (Future Work)

**Task:** Create shared DIIS utility functions (if deemed valuable)

**Decision point:** Evaluate if the ~50-100 line savings justify the refactoring effort

---

## Files Referenced in This Analysis

### Core DIIS Infrastructure
- `psi4/src/psi4/libdiis/diismanager.h`
- `psi4/src/psi4/libdiis/diismanager.cc`
- `psi4/driver/procrouting/diis.py`

### occ Module
- `psi4/src/psi4/occ/occwave.h`
- `psi4/src/psi4/occ/occ_iterations.cc`
- `psi4/src/psi4/occ/cepa_iterations.cc`

### dfocc Module
- `psi4/src/psi4/dfocc/dfocc.h`
- `psi4/src/psi4/dfocc/diis.cc` ⚠️ DEAD CODE
- `psi4/src/psi4/dfocc/occ_iterations.cc`
- `psi4/src/psi4/dfocc/ccsd_iterations.cc`
- `psi4/src/psi4/dfocc/ccd_iterations.cc`
- `psi4/src/psi4/dfocc/lccd_iterations.cc`
- `psi4/src/psi4/dfocc/remp_iterations.cc`
- `psi4/src/psi4/dfocc/update_mo.cc`

### dct Module
- `psi4/src/psi4/dct/dct.h`
- `psi4/src/psi4/dct/dct_oo_RHF.cc`
- `psi4/src/psi4/dct/dct_oo_UHF.cc`
- `psi4/src/psi4/dct/dct_compute_RHF.cc`
- `psi4/src/psi4/dct/dct_compute_UHF.cc`
- `psi4/src/psi4/dct/dct_tau_UHF.cc`
- `psi4/src/psi4/dct/dct_gradient_UHF.cc`
- `psi4/src/psi4/dct/dct_qc.cc`

### psimrcc Module
- `psi4/src/psi4/psimrcc/blas.h`
- `psi4/src/psi4/psimrcc/blas_diis.cc`

---

*Analysis completed: 2025-11-18*
*Analyst: Claude (Sonnet 4.5)*
