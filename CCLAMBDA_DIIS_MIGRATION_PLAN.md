# cclambda DIIS Migration to libdiis - Implementation Plan

**Date**: 2025-11-18
**Status**: Ready to implement
**Prerequisite**: ccenergy DIIS migration (✅ COMPLETE)

---

## Overview

This document outlines the implementation plan for migrating cclambda's custom DIIS implementation (~768 lines) to libdiis, following the proven methodology from the successful ccenergy migration.

**Goal**: Reduce cclambda DIIS code by ~66% (~768 → ~260 lines) while maintaining identical numerical behavior.

---

## Prerequisites Checklist

- [x] ccenergy DIIS migration complete and validated
- [x] All ccenergy tests passing (18/18)
- [x] ccenergy guards removed (production deployment)
- [x] libdiis DPD support confirmed (dpdfile2, dpdbuf4)
- [x] cclambda code structure analyzed
- [x] Migration feasibility confirmed (EXTREMELY HIGH)

---

## Implementation Phases

### Phase 1: RHF Implementation

**Objective**: Create libdiis-based DIIS for RHF lambda amplitudes

**File to Create**: `psi4/src/psi4/cc/cclambda/diis_RHF_libdiis.cc`

**Template** (based on ccenergy_RHF_libdiis.cc):

```cpp
/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCLAMBDA
    \brief DIIS extrapolation for RHF Lambda amplitudes using libdiis
*/

#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/psifiles.h"
#include "cclambda.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cclambda {

/*
** diis_RHF_libdiis: DIIS extrapolation for RHF Lambda amplitudes using libdiis
**
** This implementation uses libdiis/DIISManager for DIIS extrapolation,
** leveraging native DPD buffer support to eliminate manual flattening.
**
** Amplitude components (2 total):
**   - L1: LIA (singles)
**   - L2: LIjAb (doubles)
**
** Error vectors: R = L_new - L_old (computed using DPD operations)
**
** Key advantages over original implementation:
**   - No manual vector length calculation
**   - No manual flattening of DPD buffers to 1D arrays
**   - No manual B matrix construction or linear system solving
**   - Automatic storage management by libdiis
**   - Cleaner, more maintainable code
**
** Parameters:
**   iter: Iteration number
**   L_irr: Irrep of target state (for excited states)
*/

void CCLambdaWavefunction::diis_RHF_libdiis(int iter, int L_irr) {
    // Need at least 2 iterations for DIIS extrapolation
    if (iter < 2) {
        return;
    }

    auto nirreps = moinfo.nirreps;
    dpdfile2 L1_new, L1_old, R1;
    dpdbuf4 L2_new, L2_old, R2;

    /*
     * Step 1: Compute error vector for L1 (singles)
     * R1 = L1_new - L1_old using DPD file2_axpy operation
     */
    global_dpd_->file2_init(&L1_new, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&L1_old, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");

    // Create error vector as copy of L1_new, then subtract L1_old
    global_dpd_->file2_copy(&L1_new, PSIF_CC_OEI, "R1_IA");
    global_dpd_->file2_init(&R1, PSIF_CC_OEI, L_irr, 0, 1, "R1_IA");
    global_dpd_->file2_axpy(&L1_old, &R1, -1.0, 0);
    global_dpd_->file2_close(&L1_old);

    /*
     * Step 2: Compute error vector for L2 (doubles)
     * R2 = L2_new - L2_old using DPD buf4_axpy operation
     */
    global_dpd_->buf4_init(&L2_new, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_init(&L2_old, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");

    // Create error vector as copy of L2_new, then subtract L2_old
    global_dpd_->buf4_copy(&L2_new, PSIF_CC_LAMBDA, "R2_IjAb");
    global_dpd_->buf4_init(&R2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "R2_IjAb");
    global_dpd_->buf4_axpy(&L2_old, &R2, -1.0);
    global_dpd_->buf4_close(&L2_old);

    /*
     * Step 3: Add error and amplitude vectors to DIIS subspace
     * libdiis handles:
     *   - DPD → Matrix conversion
     *   - Vector storage (OnDisk policy)
     *   - B matrix construction
     *   - Subspace management (LargestError removal if full)
     */
    bool added = ccsd_diis_manager_->add_entry(
        &R1, &R2,        // Error vectors (2 components)
        &L1_new, &L2_new // Amplitude vectors (2 components)
    );

    /*
     * Step 4: Perform DIIS extrapolation if subspace is large enough
     * Requires at least 2 vectors for extrapolation
     */
    int subspace_size = ccsd_diis_manager_->subspace_size();
    if (subspace_size >= 2) {
        // libdiis computes optimal linear combination and updates amplitudes in place
        ccsd_diis_manager_->extrapolate(&L1_new, &L2_new);
        outfile->Printf("  DIIS: extrapolated with %d vectors\n", subspace_size);
    }

    /*
     * Step 5: Spin-adaptation (legacy code for compatibility)
     * RHF requires copying LIjAb → LIJAB and Lijab for spin-adapted routines
     * This preserves the original behavior
     */
    global_dpd_->file2_copy(&L1_new, PSIF_CC_LAMBDA, "New Lia");

    // Need to reinitialize L2_new with proper symmetry for spin-adaptation
    global_dpd_->buf4_close(&L2_new);
    global_dpd_->buf4_init(&L2_new, PSIF_CC_LAMBDA, L_irr, 2, 7, 0, 5, 1, "New LIjAb");
    global_dpd_->buf4_copy(&L2_new, PSIF_CC_LAMBDA, "New LIJAB");
    global_dpd_->buf4_copy(&L2_new, PSIF_CC_LAMBDA, "New Lijab");
    global_dpd_->buf4_close(&L2_new);

    /*
     * Step 6: Cleanup - close all DPD file handles
     */
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->file2_close(&L1_new);
}

}  // namespace cclambda
}  // namespace psi
```

**Estimated Lines**: ~125 (including comments and spin-adaptation)
**Actual Code**: ~60-65 lines
**Reduction**: 213 → 65 lines (69% reduction)

---

### Phase 2: ROHF Implementation

**Objective**: Create libdiis-based DIIS for ROHF lambda amplitudes

**File to Create**: `psi4/src/psi4/cc/cclambda/diis_ROHF_libdiis.cc`

**Amplitude Components** (5 total):
- L1: `LIA`, `Lia` (2 components)
- L2: `LIJAB`, `Lijab`, `LIjAb` (3 components)

**Template Structure**:
```cpp
void CCLambdaWavefunction::diis_ROHF_libdiis(int iter, int L_irr) {
    if (iter < 2) return;

    auto nirreps = moinfo.nirreps;
    dpdfile2 L1a_new, L1a_old, R1a;
    dpdfile2 L1b_new, L1b_old, R1b;
    dpdbuf4 L2aa_new, L2aa_old, R2aa;
    dpdbuf4 L2bb_new, L2bb_old, R2bb;
    dpdbuf4 L2ab_new, L2ab_old, R2ab;

    // Compute error vectors for all 5 components
    // L1a (LIA)
    global_dpd_->file2_init(&L1a_new, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&L1a_old, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_copy(&L1a_new, PSIF_CC_OEI, "R1a_IA");
    global_dpd_->file2_init(&R1a, PSIF_CC_OEI, L_irr, 0, 1, "R1a_IA");
    global_dpd_->file2_axpy(&L1a_old, &R1a, -1.0, 0);
    global_dpd_->file2_close(&L1a_old);

    // L1b (Lia)
    global_dpd_->file2_init(&L1b_new, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia");
    global_dpd_->file2_init(&L1b_old, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
    global_dpd_->file2_copy(&L1b_new, PSIF_CC_OEI, "R1b_ia");
    global_dpd_->file2_init(&R1b, PSIF_CC_OEI, L_irr, 0, 1, "R1b_ia");
    global_dpd_->file2_axpy(&L1b_old, &R1b, -1.0, 0);
    global_dpd_->file2_close(&L1b_old);

    // L2aa (LIJAB)
    global_dpd_->buf4_init(&L2aa_new, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_init(&L2aa_old, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_copy(&L2aa_new, PSIF_CC_LAMBDA, "R2aa_IJAB");
    global_dpd_->buf4_init(&R2aa, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "R2aa_IJAB");
    global_dpd_->buf4_axpy(&L2aa_old, &R2aa, -1.0);
    global_dpd_->buf4_close(&L2aa_old);

    // L2bb (Lijab)
    global_dpd_->buf4_init(&L2bb_new, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    global_dpd_->buf4_init(&L2bb_old, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_copy(&L2bb_new, PSIF_CC_LAMBDA, "R2bb_ijab");
    global_dpd_->buf4_init(&R2bb, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "R2bb_ijab");
    global_dpd_->buf4_axpy(&L2bb_old, &R2bb, -1.0);
    global_dpd_->buf4_close(&L2bb_old);

    // L2ab (LIjAb)
    global_dpd_->buf4_init(&L2ab_new, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_init(&L2ab_old, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_copy(&L2ab_new, PSIF_CC_LAMBDA, "R2ab_IjAb");
    global_dpd_->buf4_init(&R2ab, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "R2ab_IjAb");
    global_dpd_->buf4_axpy(&L2ab_old, &R2ab, -1.0);
    global_dpd_->buf4_close(&L2ab_old);

    // Add to DIIS subspace (variadic template handles 5 components)
    bool added = ccsd_diis_manager_->add_entry(
        &R1a, &R1b, &R2aa, &R2bb, &R2ab,               // Error vectors
        &L1a_new, &L1b_new, &L2aa_new, &L2bb_new, &L2ab_new  // Amplitude vectors
    );

    // Extrapolate if enough vectors
    int subspace_size = ccsd_diis_manager_->subspace_size();
    if (subspace_size >= 2) {
        ccsd_diis_manager_->extrapolate(&L1a_new, &L1b_new, &L2aa_new, &L2bb_new, &L2ab_new);
        outfile->Printf("  DIIS: extrapolated with %d vectors\n", subspace_size);
    }

    // Cleanup
    global_dpd_->file2_close(&R1a);
    global_dpd_->file2_close(&R1b);
    global_dpd_->buf4_close(&R2aa);
    global_dpd_->buf4_close(&R2bb);
    global_dpd_->buf4_close(&R2ab);
    global_dpd_->file2_close(&L1a_new);
    global_dpd_->file2_close(&L1b_new);
    global_dpd_->buf4_close(&L2aa_new);
    global_dpd_->buf4_close(&L2bb_new);
    global_dpd_->buf4_close(&L2ab_new);
}
```

**Estimated Lines**: ~150 (including comments)
**Actual Code**: ~95-100 lines
**Reduction**: 275 → 100 lines (64% reduction)

---

### Phase 3: UHF Implementation

**Objective**: Create libdiis-based DIIS for UHF lambda amplitudes

**File to Create**: `psi4/src/psi4/cc/cclambda/diis_UHF_libdiis.cc`

**Amplitude Components** (5 total - same as ROHF but different DPD spaces):
- L1: `LIA`, `Lia`
- L2: `LIJAB`, `Lijab`, `LIjAb`

**Key Difference from ROHF**: DPD index spaces
- UHF L1a: `(0, 1)` → Alpha occupied/virtual
- UHF L1b: `(2, 3)` → Beta occupied/virtual
- ROHF uses same spaces for both

**Template**: Nearly identical to ROHF, just update file2_init calls:
```cpp
// UHF uses different orbital space indices
global_dpd_->file2_init(&L1a_new, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");  // Alpha
global_dpd_->file2_init(&L1b_new, PSIF_CC_LAMBDA, L_irr, 2, 3, "New Lia");  // Beta (different!)
```

**Estimated Lines**: ~175 (including comments)
**Actual Code**: ~115-120 lines
**Reduction**: 280 → 120 lines (57% reduction)

---

### Phase 4: Integration

**Objective**: Integrate new libdiis implementations into cclambda build system

#### 4.1 Header File

**Modify**: `psi4/src/psi4/cc/cclambda/cclambda.h`

**Add after line 68** (`void diis(int, int);`):
```cpp
    // libdiis-based DIIS implementations
    void diis_RHF_libdiis(int, int);   // libdiis implementation for RHF
    void diis_ROHF_libdiis(int, int);  // libdiis implementation for ROHF
    void diis_UHF_libdiis(int, int);   // libdiis implementation for UHF
```

#### 4.2 Main Source File

**Modify**: `psi4/src/psi4/cc/cclambda/cclambda.cc`

**Add include** (near top, after existing includes):
```cpp
#include "psi4/libdiis/diismanager.h"
```

**Add DIIS initialization** in `compute_energy()` function (after `get_params(options_);` at line 125):
```cpp
    // Initialize DIISManager for lambda amplitude extrapolation
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
```

#### 4.3 DIIS Dispatcher

**Replace entire file**: `psi4/src/psi4/cc/cclambda/diis.cc`

**New content**:
```cpp
/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCLAMBDA
    \brief DIIS extrapolation dispatcher for Lambda amplitudes
*/

#include "cclambda.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cclambda {

/*
** DIIS: Direct inversion in the iterative subspace routine to
** accelerate convergence of the Lambda amplitude equations.
**
** Implementation uses libdiis/DIISManager for all reference types,
** leveraging centralized DIIS infrastructure and native DPD buffer support.
**
** Original custom implementations replaced 2025-11-18:
** - Eliminates ~768 lines of duplicate DIIS code
** - Uses native DPD operations (file2_axpy, buf4_axpy)
** - Automatic B matrix construction and linear system solving
** - Consistent DIIS behavior with ccenergy, occ, dfocc modules
**
** Migration based on proven ccenergy DIIS consolidation (100% successful).
*/

void CCLambdaWavefunction::diis(int iter, int L_irr) {
    // Use libdiis implementation for all reference types
    if (params.ref == 0)
        diis_RHF_libdiis(iter, L_irr);
    else if (params.ref == 1)
        diis_ROHF_libdiis(iter, L_irr);
    else if (params.ref == 2)
        diis_UHF_libdiis(iter, L_irr);
}

}  // namespace cclambda
}  // namespace psi
```

#### 4.4 Build System

**Modify**: `psi4/src/psi4/cc/cclambda/CMakeLists.txt`

**Find the source files list and add**:
```cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/diis_RHF_libdiis.cc   # libdiis implementation for RHF
    ${CMAKE_CURRENT_SOURCE_DIR}/diis_ROHF_libdiis.cc  # libdiis implementation for ROHF
    ${CMAKE_CURRENT_SOURCE_DIR}/diis_UHF_libdiis.cc   # libdiis implementation for UHF
```

---

### Phase 5: Testing and Validation

**Objective**: Verify migration produces identical results to original

#### 5.1 Basic Functional Test

**Test**: RHF ground state lambda
```bash
pytest psi4/tests/cclambda -k "cc1" -v
```

**Expected**:
- Lambda iterations converge
- Same iteration count as original
- Energies match to < 1e-9 Hartree
- Output shows: "Using libdiis for Lambda DIIS extrapolation (RHF)"

#### 5.2 Comprehensive Test Suite

**Run full cclambda test suite**:
```bash
pytest psi4/tests/cclambda -v
pytest psi4/tests/cc* -k lambda -v
```

**Coverage**:
- RHF ground state
- ROHF ground state
- UHF ground state
- Excited states (EOM-CCSD)
- Analytical gradients
- Properties

#### 5.3 Validation Criteria

| Metric | Target | How to Check |
|--------|--------|--------------|
| Energy accuracy | < 1e-9 Hartree | Compare with reference |
| Iteration count | Identical | Check output logs |
| Convergence pattern | Matching | Compare iteration energies |
| Performance | Within 5% | Time measurements |
| Memory | No increase | Monitor during tests |

#### 5.4 Regression Testing

**Compare side-by-side** (if needed):
1. Build with original code (git stash changes)
2. Run test, save output
3. Build with new code
4. Run test, compare outputs
5. Verify identical energies and iterations

---

### Phase 6: Documentation and Cleanup

**Objective**: Document migration and remove old code

#### 6.1 Update Documentation

**Create**: `CCLAMBDA_DIIS_MIGRATION_SUMMARY.md`
- Migration results
- Test validation summary
- Code reduction metrics
- Performance comparison

**Update**: Main project documentation
- Note libdiis is now used for lambda DIIS
- Reference CCLAMBDA_DIIS_ANALYSIS.md for details

#### 6.2 Code Cleanup

Once testing is complete and validated:

**Delete original implementation**:
- The bulk of `diis.cc` (~768 lines) replaced by dispatcher (~60 lines)
- Net reduction: ~708 lines

**Final state**:
- `diis.cc`: Dispatcher only (~60 lines)
- `diis_RHF_libdiis.cc`: ~65 lines
- `diis_ROHF_libdiis.cc`: ~100 lines
- `diis_UHF_libdiis.cc`: ~120 lines
- **Total**: ~345 lines (was ~768 lines)
- **Reduction**: 55% including new implementations

#### 6.3 Commit Strategy

**Commit 1**: Add RHF implementation
```
Add libdiis-based DIIS for cclambda RHF

- Create diis_RHF_libdiis.cc with DPD-native operations
- Reduces RHF DIIS from ~213 to ~65 lines
- Uses ccsd_diis_manager_ from parent CCEnergyWavefunction
```

**Commit 2**: Add ROHF implementation
```
Add libdiis-based DIIS for cclambda ROHF

- Create diis_ROHF_libdiis.cc handling 5 amplitude components
- Reduces ROHF DIIS from ~275 to ~100 lines
- Follows proven ccenergy ROHF pattern
```

**Commit 3**: Add UHF implementation
```
Add libdiis-based DIIS for cclambda UHF

- Create diis_UHF_libdiis.cc with proper UHF orbital spaces
- Reduces UHF DIIS from ~280 to ~120 lines
- Completes cclambda libdiis migration
```

**Commit 4**: Integration and testing
```
Integrate cclambda libdiis DIIS and validate

- Update cclambda.h, cclambda.cc, CMakeLists.txt
- Initialize DIISManager in compute_energy()
- Replace diis.cc with dispatcher
- All tests passing (energies, gradients, properties)
```

**Commit 5**: Documentation
```
Document cclambda DIIS migration to libdiis

- Create CCLAMBDA_DIIS_MIGRATION_SUMMARY.md
- Total reduction: ~768 → ~345 lines (55%)
- Combined with ccenergy: ~1,748 → ~605 lines (65%)
```

---

## Success Metrics

### Code Quality
- [  ] Total line reduction: ~55% (768 → 345 lines)
- [  ] RHF reduction: ~69% (213 → 65 lines)
- [  ] ROHF reduction: ~64% (275 → 100 lines)
- [  ] UHF reduction: ~57% (280 → 120 lines)
- [  ] Code is cleaner and more maintainable
- [  ] Follows established ccenergy pattern

### Functional Correctness
- [  ] All cclambda tests pass
- [  ] Energies match reference (< 1e-9 Hartree)
- [  ] Iteration counts identical
- [  ] Convergence patterns match
- [  ] Gradients validated (if applicable)
- [  ] Properties validated (if applicable)

### Performance
- [  ] Runtime within 5% of original
- [  ] Memory usage unchanged
- [  ] No performance regressions

### Integration
- [  ] Clean compilation (no warnings)
- [  ] Proper CMake integration
- [  ] DIISManager initialization works
- [  ] Inherits from CCEnergyWavefunction correctly

---

## Risk Mitigation

### Potential Issues and Solutions

**Issue 1**: L_irr parameter not handled correctly
- **Risk**: Low
- **Mitigation**: DPD operations already use L_irr in indexing
- **Test**: Run excited state calculations

**Issue 2**: Spin-adaptation breaks
- **Risk**: Low
- **Mitigation**: Preserved legacy copy operations in RHF
- **Test**: Verify spin-adapted quantities

**Issue 3**: Different convergence behavior
- **Risk**: Very Low (proven with ccenergy)
- **Mitigation**: Same DIIS algorithm, same parameters
- **Test**: Compare iteration-by-iteration convergence

**Issue 4**: Performance regression
- **Risk**: Very Low (no overhead in ccenergy)
- **Mitigation**: OnDisk storage, no copying overhead
- **Test**: Benchmark before/after

---

## Timeline

| Phase | Task | Estimated Time | Cumulative |
|-------|------|----------------|------------|
| 1 | RHF implementation | 1.5 hours | 1.5 hours |
| 2 | ROHF implementation | 1.5 hours | 3 hours |
| 3 | UHF implementation | 1.5 hours | 4.5 hours |
| 4 | Integration | 1 hour | 5.5 hours |
| 5 | Testing & validation | 2 hours | 7.5 hours |
| 6 | Documentation | 1 hour | 8.5 hours |
| **Total** | **Complete migration** | **8.5 hours** | **8.5 hours** |

**Note**: Faster than ccenergy (~12 hours) due to:
- Proven template from ccenergy
- Established testing methodology
- No exploratory work needed

---

## Next Steps

1. **Create RHF implementation** (Phase 1)
   - Copy ccenergy_RHF_libdiis.cc template
   - Update T → L amplitude names
   - Add L_irr parameter
   - Add spin-adaptation code

2. **Test RHF** immediately
   - Run cc1 test
   - Verify convergence
   - Check energy accuracy

3. **Proceed to ROHF** (Phase 2)
   - Follow ccenergy ROHF pattern
   - Handle 5 components
   - Test with ROHF cases

4. **Complete UHF** (Phase 3)
   - Adjust orbital space indices
   - Test with UHF cases

5. **Integrate and validate** (Phases 4-5)
   - Full test suite
   - Performance benchmarking

6. **Document and deploy** (Phase 6)
   - Final summary
   - Commit to repository

---

## Conclusion

This migration follows the proven, successful methodology from ccenergy DIIS consolidation. The implementation is straightforward, low-risk, and delivers substantial code reduction (~55%) with improved maintainability.

**Combined Impact** (ccenergy + cclambda):
- **Before**: ~1,748 lines of duplicate DIIS code
- **After**: ~605 lines of clean, maintainable code
- **Reduction**: 65% (~1,143 lines eliminated)

**Recommendation**: **PROCEED** with implementation following this plan.

---

**Plan Author**: Claude (Anthropic AI)
**Plan Date**: 2025-11-18
**Status**: Ready to implement
**Expected Completion**: ~8.5 hours total effort
