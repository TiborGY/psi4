# cclambda DIIS Consolidation Analysis

**Date**: 2025-11-18
**Investigator**: Claude (Anthropic AI)
**Context**: Following successful ccenergy DIIS migration to libdiis

---

## Executive Summary

The cclambda module contains a custom DIIS implementation (~768 lines across three reference types) that is **structurally IDENTICAL** to the ccenergy DIIS code we just successfully migrated. The only difference is that cclambda operates on Lambda amplitudes (L) instead of T amplitudes.

**Key Question**: Can cclambda DIIS be migrated to libdiis like ccenergy?

**Short Answer**: **YES - Extremely High Feasibility** - The proven ccenergy methodology can be applied directly with minimal modifications.

**Recommendation**: **HIGH PRIORITY** - Immediate migration using established pattern from ccenergy.

---

## Current Implementation

### File Structure
- **Location**: `psi4/src/psi4/cc/cclambda/diis.cc`
- **Size**: 858 lines total (includes license header and comments)
- **Actual DIIS code**: ~768 lines of custom implementation
- **Class**: `CCLambdaWavefunction` (inherits from `CCEnergyWavefunction`)

### DIIS Function Signature

```cpp
void CCLambdaWavefunction::diis(int iter, int L_irr)
```

- `iter`: Iteration number
- `L_irr`: Irrep of the target state (for excited state calculations)

### Reference Type Coverage

**RHF**: Lines 81-294 (~213 lines)
- L1 amplitudes: `LIA` (1 component)
- L2 amplitudes: `LIjAb` (1 component)
- **Total**: 2 amplitude components

**ROHF**: Lines 295-570 (~275 lines)
- L1 amplitudes: `LIA`, `Lia` (2 components)
- L2 amplitudes: `LIJAB`, `Lijab`, `LIjAb` (3 components)
- **Total**: 5 amplitude components

**UHF**: Lines 571-851 (~280 lines)
- L1 amplitudes: `LIA`, `Lia` (2 components)
- L2 amplitudes: `LIJAB`, `Lijab`, `LIjAb` (3 components)
- **Total**: 5 amplitude components

---

## Code Structure Analysis

### RHF Implementation (Lines 81-294)

**Step 1: Vector Length Calculation** (Lines 82-90)
```cpp
global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
for (h = 0; h < nirreps; h++) {
    vector_length += L1.params->rowtot[h] * L1.params->coltot[h ^ L_irr];
    vector_length += L2.params->rowtot[h] * L2.params->coltot[h ^ L_irr];
}
```

**Step 2: Error Vector Computation** (Lines 96-136)
- Manual flattening of DPD buffers to 1D array
- Computes `R = L_new - L_old` for both L1 and L2
- Stores to PSIO file `PSIF_CC_DIIS_ERR`

**Step 3: Amplitude Storage** (Lines 142-166)
- Flattens current L amplitudes to 1D array
- Stores to PSIO file `PSIF_CC_DIIS_AMP`

**Step 4: B Matrix Construction** (Lines 178-226)
- Reads error vectors from disk
- Computes dot products: `B[p][q] = error[p] · error[q]`
- Scales by maximum value for numerical stability
- Adds constraint row/column: `B[i][nvector] = -1`

**Step 5: Linear System Solving** (Lines 237-243)
- Solves `B·c = [-1]` using LAPACK `C_DGESV`
- Gets DIIS coefficients in array `c`

**Step 6: Extrapolation** (Lines 245-258)
- Reads old amplitude vectors from disk
- Computes linear combination: `L_new = Σ c[p] * L_old[p]`

**Step 7: Unpack to DPD** (Lines 260-281)
- Manually unpacks 1D array back to DPD file2/buf4 structures
- Writes updated amplitudes to DPD files

**Step 8: Spin-Adaptation** (Lines 269, 284-287)
- RHF-specific: Copies `LIjAb` → `LIJAB` and `Lijab`
- Legacy code for spin-adaptation

### ROHF Implementation (Lines 295-570)

Similar structure to RHF but with **5 amplitude components**:
- L1: `LIA`, `Lia`
- L2: `LIJAB`, `Lijab`, `LIjAb`

Key differences:
- Larger vector length calculation (lines 297-309)
- All 5 components manually flattened to error vector (lines 314-393)
- All 5 components stored as amplitudes (lines 399-447)
- Uses `flin()` instead of `C_DGESV` (line 501)
- Manual unpacking of all 5 components (lines 515-563)

### UHF Implementation (Lines 571-851)

Identical structure to ROHF with **5 amplitude components**:
- L1: `LIA`, `Lia`
- L2: `LIJAB`, `Lijab`, `LIjAb`

Note: UHF L1 components use different DPD index spaces (lines 574-575, 613-626)

---

## Comparison with ccenergy

| Aspect | ccenergy (Original) | ccenergy (libdiis) | cclambda (Current) |
|--------|---------------------|--------------------|--------------------|
| **Amplitude Type** | T (tau) amplitudes | T amplitudes | L (lambda) amplitudes |
| **Data Structure** | DPD buffers | DPD buffers | DPD buffers |
| **Algorithm** | Identical DIIS | Identical DIIS | Identical DIIS |
| **Manual Flattening** | Yes | No (DPD ops) | Yes |
| **Manual B Matrix** | Yes | No (automatic) | Yes |
| **Manual Solving** | Yes (DGESV/flin) | No (automatic) | Yes (DGESV/flin) |
| **Code Lines (Total)** | ~980 lines | ~260 lines | ~768 lines |
| **RHF Lines** | 258 | 60 | 213 |
| **ROHF Lines** | 357 | 90 | 275 |
| **UHF Lines** | 365 | 110 | 280 |
| **libdiis Support** | Native (DPD) | Yes | Native (DPD) |

**Key Insight**: The algorithmic structure is **IDENTICAL**. Only the amplitude names differ (T → L).

---

## Migration Feasibility

### Technical Feasibility: **EXTREMELY HIGH** ✅

**Why Migration is Straightforward:**

1. **Proven Pattern**: ccenergy migration was successful for all three reference types
2. **DPD Buffers**: cclambda uses dpdfile2/dpdbuf4 just like ccenergy
3. **Identical Algorithm**: DIIS logic is identical, only amplitude names differ
4. **Class Inheritance**: `CCLambdaWavefunction` inherits from `CCEnergyWavefunction`
   - Already has `ccsd_diis_manager_` member variable
   - Can reuse parent class infrastructure

5. **libdiis Compatibility**:
   - Native support for dpdfile2 ✅
   - Native support for dpdbuf4 ✅
   - Handles L_irr symmetry automatically via DPD ✅

### Code Reduction Estimate

Based on ccenergy results (73% reduction):

| Reference | Current Lines | Estimated New Lines | Reduction |
|-----------|---------------|---------------------|-----------|
| RHF       | 213           | ~60                 | 72%       |
| ROHF      | 275           | ~90                 | 67%       |
| UHF       | 280           | ~110                | 61%       |
| **Total** | **~768**      | **~260**            | **66%**   |

**Net Savings**: ~500 lines of complex, duplicate DIIS code eliminated

---

## Key Differences from ccenergy

### 1. Amplitude Names
- ccenergy: `tIA`, `New tIA`, `tIjAb`, `New tIjAb`
- cclambda: `LIA`, `New LIA`, `LIjAb`, `New LIjAb`

### 2. L_irr Parameter
- cclambda DIIS takes `L_irr` parameter (target state irrep)
- This is already handled by DPD operations (`h ^ L_irr` in loops)
- libdiis will handle this transparently via DPD buffers

### 3. Spin-Adaptation
- RHF has legacy spin-adaptation code (lines 269, 284-287)
- Copies `LIjAb` → `LIJAB`, `Lijab`
- Must be preserved in libdiis implementation

### 4. PSIO Files
- ccenergy uses: `PSIF_CC_TAMPS`, `PSIF_CC_OEI`
- cclambda uses: `PSIF_CC_LAMBDA`, `PSIF_CC_TMP0`
- Just need to update file IDs in libdiis implementation

### 5. DIISManager Initialization
- ccenergy: In `CCEnergyWavefunction::CCEnergyWavefunction()`
- cclambda: Should be in `CCLambdaWavefunction::compute_energy()` after `get_params()`
- Can reuse parent's `ccsd_diis_manager_` member

---

## Implementation Plan

### Phase 1: RHF Implementation

**Create**: `psi4/src/psi4/cc/cclambda/diis_RHF_libdiis.cc`

**Structure** (based on ccenergy pattern):
```cpp
void CCLambdaWavefunction::diis_RHF_libdiis(int iter, int L_irr) {
    if (iter < 2) return;

    auto nirreps = moinfo.nirreps;
    dpdfile2 L1_new, L1_old, R1;
    dpdbuf4 L2_new, L2_old, R2;

    // Step 1: Compute error vectors using DPD operations
    global_dpd_->file2_init(&L1_new, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&L1_old, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_copy(&L1_new, PSIF_CC_OEI, "R1_IA");
    global_dpd_->file2_init(&R1, PSIF_CC_OEI, L_irr, 0, 1, "R1_IA");
    global_dpd_->file2_axpy(&L1_old, &R1, -1.0, 0);
    global_dpd_->file2_close(&L1_old);

    global_dpd_->buf4_init(&L2_new, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_init(&L2_old, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_copy(&L2_new, PSIF_CC_LAMBDA, "R2_IjAb");
    global_dpd_->buf4_init(&R2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "R2_IjAb");
    global_dpd_->buf4_axpy(&L2_old, &R2, -1.0);
    global_dpd_->buf4_close(&L2_old);

    // Step 2: Add to DIIS subspace
    bool added = ccsd_diis_manager_->add_entry(&R1, &R2, &L1_new, &L2_new);

    // Step 3: Extrapolate if enough vectors
    int subspace_size = ccsd_diis_manager_->subspace_size();
    if (subspace_size >= 2) {
        ccsd_diis_manager_->extrapolate(&L1_new, &L2_new);
        outfile->Printf("  DIIS: extrapolated with %d vectors\n", subspace_size);
    }

    // Step 4: Spin-adaptation (preserve legacy behavior)
    global_dpd_->file2_copy(&L1_new, PSIF_CC_LAMBDA, "New Lia");

    global_dpd_->buf4_init(&L2_new, PSIF_CC_LAMBDA, L_irr, 2, 7, 0, 5, 1, "New LIjAb");
    global_dpd_->buf4_copy(&L2_new, PSIF_CC_LAMBDA, "New LIJAB");
    global_dpd_->buf4_copy(&L2_new, PSIF_CC_LAMBDA, "New Lijab");
    global_dpd_->buf4_close(&L2_new);

    // Step 5: Cleanup
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->file2_close(&L1_new);
}
```

**Estimated Lines**: ~65-70 (including spin-adaptation)

### Phase 2: ROHF Implementation

**Create**: `psi4/src/psi4/cc/cclambda/diis_ROHF_libdiis.cc`

**Amplitude Components**:
- L1: `LIA`, `Lia`
- L2: `LIJAB`, `Lijab`, `LIjAb`

**Pattern**: Identical to ccenergy ROHF, just with L instead of T

**Estimated Lines**: ~95-100

### Phase 3: UHF Implementation

**Create**: `psi4/src/psi4/cc/cclambda/diis_UHF_libdiis.cc`

**Amplitude Components**: Same as ROHF

**Pattern**: Identical to ccenergy UHF, just with L instead of T

**Estimated Lines**: ~115-120

### Phase 4: Integration

**Modify**: `psi4/src/psi4/cc/cclambda/cclambda.h`
```cpp
// Add function declarations
void diis_RHF_libdiis(int, int);   // libdiis implementation for RHF
void diis_ROHF_libdiis(int, int);  // libdiis implementation for ROHF
void diis_UHF_libdiis(int, int);   // libdiis implementation for UHF
```

**Modify**: `psi4/src/psi4/cc/cclambda/cclambda.cc`
```cpp
#include "psi4/libdiis/diismanager.h"

double CCLambdaWavefunction::compute_energy() {
    // ... existing code ...
    get_params(options_);

    // Initialize DIISManager for lambda extrapolation
    if (params.diis) {
        const char* ref_label = (params.ref == 0) ? "RHF" : (params.ref == 1) ? "ROHF" : "UHF";
        std::string diis_label = std::string("Lambda DIIS ") + ref_label;

        ccsd_diis_manager_ = std::make_shared<DIISManager>(
            8,                                          // max 8 vectors
            diis_label,                                 // label with reference type
            DIISManager::RemovalPolicy::LargestError,  // removal policy
            DIISManager::StoragePolicy::OnDisk         // storage policy
        );
        outfile->Printf("  Using libdiis for Lambda DIIS extrapolation (%s)\n", ref_label);
    }

    // ... rest of compute_energy() ...
}
```

**Modify**: `psi4/src/psi4/cc/cclambda/diis.cc`
```cpp
void CCLambdaWavefunction::diis(int iter, int L_irr) {
    // Use libdiis implementation for all reference types
    if (params.ref == 0)
        diis_RHF_libdiis(iter, L_irr);
    else if (params.ref == 1)
        diis_ROHF_libdiis(iter, L_irr);
    else if (params.ref == 2)
        diis_UHF_libdiis(iter, L_irr);
}
```

**Modify**: `psi4/src/psi4/cc/cclambda/CMakeLists.txt`
```cmake
${CMAKE_CURRENT_SOURCE_DIR}/diis_RHF_libdiis.cc
${CMAKE_CURRENT_SOURCE_DIR}/diis_ROHF_libdiis.cc
${CMAKE_CURRENT_SOURCE_DIR}/diis_UHF_libdiis.cc
```

### Phase 5: Testing

**Leverage Existing Tests**:
- cclambda has extensive test suite
- Check convergence behavior matches original
- Verify energies, gradients, properties

**Test Cases**:
- RHF ground state
- ROHF ground state
- UHF ground state
- Excited states (EOM-CCSD)
- Analytical gradients

### Phase 6: Production Deployment

**Remove Old Code**:
- Delete original `diis.cc` implementation (~768 lines)
- Keep dispatcher with libdiis routing only
- Update documentation

---

## Advantages of Migration

### Immediate Benefits (cclambda)

1. **Code Reduction**: ~500 lines eliminated (66% reduction)
2. **Maintainability**: Centralized DIIS in libdiis
3. **Consistency**: Matches ccenergy, occ, dfocc approach
4. **Reduced Duplication**: Single DIIS algorithm for all CC modules
5. **Proven Approach**: Same methodology as successful ccenergy migration

### Broader Benefits (Combined with ccenergy)

**Total CC DIIS Code Consolidated**: ~1,748 lines → ~520 lines
- ccenergy: ~980 → ~260 lines
- cclambda: ~768 → ~260 lines

**Remaining Opportunities**:
- ccresponse: ~300 lines potential reduction
- **Grand Total Potential**: ~2,200 → ~650 lines (70% reduction)

---

## Risk Assessment

### Technical Risks: **VERY LOW** ✅

| Risk | Impact | Likelihood | Mitigation | Status |
|------|--------|------------|------------|--------|
| Different convergence | High | Very Low | Proven with ccenergy | ✅ No issues expected |
| L_irr handling | Medium | Very Low | DPD handles transparently | ✅ Native support |
| Spin-adaptation | Medium | Very Low | Preserve legacy code | ✅ Easy to maintain |
| Performance | Medium | Very Low | Same as ccenergy | ✅ No regression expected |
| Integration | Low | Very Low | Inherits from CCEnergy | ✅ Simple integration |

**Overall Risk**: **VERY LOW** - Identical to ccenergy migration which was 100% successful

### Performance Considerations

**Expected**: Identical performance to ccenergy
- No overhead from libdiis (proven with ccenergy)
- OnDisk storage policy (same as original)
- DPD operations are native (no copying)

---

## Timeline Estimate

| Phase | Task | Estimated Time |
|-------|------|----------------|
| 1 | RHF implementation | 2 hours |
| 2 | ROHF implementation | 2 hours |
| 3 | UHF implementation | 2 hours |
| 4 | Integration & build | 1 hour |
| 5 | Testing & validation | 2 hours |
| 6 | Documentation | 1 hour |
| **Total** | **Complete migration** | **10 hours** |

**Note**: Much faster than ccenergy (which was exploratory). We now have:
- Proven template from ccenergy
- Understanding of libdiis integration
- Clear testing methodology

---

## Recommendation

### Priority: **HIGH - IMMEDIATE** ✅

**Rationale**:
1. **Proven Methodology**: ccenergy migration was 100% successful
2. **High Code Reduction**: ~500 lines (66%) with minimal effort
3. **Low Risk**: Identical structure to validated ccenergy approach
4. **Natural Continuation**: Builds directly on ccenergy work
5. **Quick Win**: ~10 hours total effort for substantial benefit

### Suggested Sequence

1. ✅ **ccenergy** - COMPLETE (validated, guards removed)
2. 🎯 **cclambda** - **NEXT TARGET** (this analysis)
3. ⏭️ **ccresponse** - Follow-up (~300 lines potential)

### Success Criteria

Same as ccenergy:
- ✅ Energy accuracy < 1e-9 Hartree
- ✅ Iteration count identical to original
- ✅ Convergence pattern matching
- ✅ Performance within 5%
- ✅ All regression tests pass
- ✅ No memory leaks or issues

---

## Technical Notes

### L_irr Handling

The `L_irr` parameter (irrep of target state) is used throughout DPD operations:
```cpp
vector_length += L1.params->rowtot[h] * L1.params->coltot[h ^ L_irr];
```

**libdiis Compatibility**:
- libdiis operates on DPD buffers directly
- DPD buffers already contain `L_irr` information
- No special handling needed - transparent pass-through ✅

### Spin-Adaptation

RHF has legacy spin-adaptation (lines 269, 284-287):
```cpp
global_dpd_->file2_copy(&L1a, PSIF_CC_LAMBDA, "New Lia");

global_dpd_->buf4_init(&L2a, PSIF_CC_LAMBDA, L_irr, 2, 7, 0, 5, 1, "New LIjAb");
global_dpd_->buf4_copy(&L2a, PSIF_CC_LAMBDA, "New LIJAB");
global_dpd_->buf4_copy(&L2a, PSIF_CC_LAMBDA, "New Lijab");
```

**Preservation**: This must be maintained in libdiis implementation to avoid breaking spin-adapted code.

### ROHF vs UHF

Both handle 5 components but with different DPD index conventions:
- **ROHF**: Uses same orbital spaces for alpha/beta
- **UHF**: Uses separate orbital spaces (lines 574-575)

**Implementation**: Follow ccenergy pattern - nearly identical code with proper labels.

---

## Conclusion

The cclambda DIIS implementation is **perfectly suited** for libdiis migration:

✅ **Identical algorithm** to successfully migrated ccenergy
✅ **Same data structures** (dpdfile2, dpdbuf4)
✅ **Proven methodology** from ccenergy applies directly
✅ **Low risk** with high benefit (66% code reduction)
✅ **Fast implementation** (~10 hours) using established pattern

**Recommendation**: **PROCEED IMMEDIATELY** with cclambda migration following ccenergy template.

---

**Analysis Complete**: 2025-11-18
**Next Action**: Implement RHF libdiis version as proof of concept

**Estimated Total Savings (ccenergy + cclambda)**: ~1,220 lines of duplicate DIIS code eliminated
