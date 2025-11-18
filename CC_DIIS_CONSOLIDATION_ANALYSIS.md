# CC Module DIIS Consolidation Analysis

## Executive Summary

This analysis investigates the feasibility of consolidating DIIS implementations in the coupled cluster (`cc/`) modules to use the centralized `libdiis/DIISManager`. Unlike the occ/dfocc modules which already use libdiis, **the cc modules have substantial custom DIIS implementations totaling ~2,238 lines of code**.

### Key Findings:

1. ✅ **libdiis DOES support dpdbuf4 and dpdfile2** - Technical compatibility confirmed
2. ⚠️ **~2,238 lines of custom DIIS code** across 6 files in cc modules
3. 💡 **Consolidation is technically feasible** but requires careful refactoring
4. **Estimated savings: 1,500-1,800 lines** after accounting for necessary wrapper code
5. **Risk: Medium** - cc modules are critical and widely used

---

## Code Volume Analysis

### Current DIIS Implementations

| Module | File | Lines | Purpose |
|--------|------|-------|---------|
| ccenergy | diis.cc | 118 | Dispatcher + conditioning utilities |
| ccenergy | diis_RHF.cc | 258 | RHF-specific DIIS |
| ccenergy | diis_ROHF.cc | 357 | ROHF-specific DIIS |
| ccenergy | diis_UHF.cc | 365 | UHF-specific DIIS |
| cclambda | diis.cc | 857 | Lambda amplitude DIIS (all references) |
| ccresponse | diis.cc | 283 | Response amplitude DIIS |
| **TOTAL** | | **2,238** | |

### Breakdown by Module

**ccenergy**: 1,098 lines
- Handles T1 and T2 amplitudes
- Separate implementations for RHF, ROHF, UHF
- Custom conditioning in `diis_invert_B()`

**cclambda**: 857 lines
- Handles L1 and L2 lambda amplitudes
- Includes RHF, ROHF, and UHF in single file
- Uses irrep parameter for excited states

**ccresponse**: 283 lines
- Handles response amplitudes
- Works with perturbation-dependent vectors
- Uses omega (frequency) parameter

---

## Technical Analysis

### Current Implementation Pattern

All cc DIIS implementations follow the same pattern:

```cpp
void diis(int iter) {
    // 1. Calculate vector length by summing DPD buffer sizes
    dpdfile2 T1;
    dpdbuf4 T2;
    int vector_length = 0;
    for (h = 0; h < nirreps; h++) {
        vector_length += T1.params->rowtot[h] * T1.params->coltot[h];
        vector_length += T2.params->rowtot[h] * T2.params->coltot[h];
    }

    // 2. Manually flatten DPD buffers into 1D error vector
    double **error = dpd_block_matrix(1, vector_length);
    int word = 0;
    for (h = 0; h < nirreps; h++)
        for (row = 0; row < T1.params->rowtot[h]; row++)
            for (col = 0; col < T1.params->coltot[h]; col++)
                error[0][word++] = T1_new[h][row][col] - T1_old[h][row][col];
    // ... same for T2

    // 3. Store error vector to PSIO file
    psio_write(PSIF_CC_DIIS_ERR, "DIIS Error Vectors", ...);

    // 4. Store amplitude vector to PSIO file
    psio_write(PSIF_CC_DIIS_AMP, "DIIS Amplitude Vectors", ...);

    // 5. Build B matrix by reading back error vectors
    B = block_matrix(nvector+1, nvector+1);
    for (p = 0; p < nvector; p++) {
        psio_read(PSIF_CC_DIIS_ERR, ...);
        for (q = 0; q <= p; q++) {
            psio_read(PSIF_CC_DIIS_ERR, ...);
            B[p][q] = dot_product(error_p, error_q);
        }
    }

    // 6. Solve linear system with custom conditioning
    diis_invert_B(B, C, nvector+1, tolerance);

    // 7. Extrapolate new amplitudes
    for (p = 0; p < nvector; p++) {
        psio_read(PSIF_CC_DIIS_AMP, ...);
        new_amps += C[p] * old_amps[p];
    }

    // 8. Unpack 1D vector back to DPD buffers
    word = 0;
    for (h = 0; h < nirreps; h++)
        for (row = 0; row < T1.params->rowtot[h]; row++)
            for (col = 0; col < T1.params->coltot[h]; col++)
                T1[h][row][col] = new_amps[word++];
}
```

### libdiis Capabilities

**Key discovery**: libdiis DOES support DPD buffers!

From `psi4/driver/procrouting/diis.py`:
```python
# Line 24-25: axpy operation
elif isinstance(y, (core.dpdbuf4, core.dpdfile2)):
    y.axpy_matrix(x, alpha)

# Line 47-48: template helper
elif isinstance(arg, (core.Matrix, core.dpdfile2, core.dpdbuf4)):
    template.append([arg.rowdim(), arg.coldim()])

# Line 66: documentation
# Currently supported storage types: ..., Psi.dpdfile2, Psi.dpdbuf4, ...

# Line 124-125: copier
elif isinstance(x, (core.dpdbuf4, core.dpdfile2)):
    copy = core.Matrix(x)
```

**How libdiis handles DPD buffers:**
1. Converts dpdbuf4/dpdfile2 to Matrix internally
2. Stores as Matrix in PSIO
3. Performs all DIIS operations on Matrix representation
4. Can extrapolate back to DPD format

### Proposed Migration Pattern

Using libdiis, the code would simplify to:

```cpp
void diis_RHF(int iter) {
    // 1. Open DPD buffers (no manual flattening needed!)
    dpdfile2 T1_new, T1_old;
    dpdbuf4 T2_new, T2_old;
    dpdfile2 R1;  // error vectors
    dpdbuf4 R2;

    // 2. Compute error vectors in DPD format
    global_dpd_->file2_init(&R1, ...);
    global_dpd_->file2_init(&T1_new, ...);
    global_dpd_->file2_init(&T1_old, ...);
    // R1 = T1_new - T1_old (using DPD operations)

    global_dpd_->buf4_init(&R2, ...);
    global_dpd_->buf4_init(&T2_new, ...);
    global_dpd_->buf4_init(&T2_old, ...);
    // R2 = T2_new - T2_old (using DPD operations)

    // 3. Add to DIIS (libdiis handles everything else!)
    diis_manager->add_entry(&R1, &R2, &T1_new, &T2_new);

    // 4. Extrapolate (libdiis handles everything!)
    if (diis_manager->subspace_size() >= min_diis_vectors) {
        diis_manager->extrapolate(&T1_new, &T2_new);
    }

    // 5. Close buffers
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->file2_close(&T1_new);
    global_dpd_->buf4_close(&T2_new);

    // That's it! ~20 lines vs ~250 lines
}
```

---

## Detailed Module Analysis

### 1. ccenergy Module

**Current Structure:**
- `diis.cc` (118 lines) - Dispatcher and utilities
  - `diis()` - Calls reference-specific function
  - `diis_invert_B()` - Custom conditioning for B matrix
- `diis_RHF.cc` (258 lines) - RHF amplitudes
- `diis_ROHF.cc` (357 lines) - ROHF amplitudes
- `diis_UHF.cc` (365 lines) - UHF amplitudes

**Commonalities:**
- All use fixed subspace size of 8 vectors
- All use PSIF_CC_DIIS_ERR and PSIF_CC_DIIS_AMP
- All use same B matrix construction
- All use same custom conditioning

**Differences:**
- Number of error/amplitude components:
  - RHF: T1(ia), T2(ijab)
  - ROHF: T1(IA), T1(ia), T2(IJAB), T2(ijab), T2(IjAb)
  - UHF: T1(IA), T1(ia), T2(IJAB), T2(ijab), T2(IjAb)

**Migration Strategy:**
```cpp
// Create DIISManager instance (once per calculation)
std::shared_ptr<DIISManager> ccsd_diis;

void CCEnergyWavefunction::diis_init() {
    ccsd_diis = std::make_shared<DIISManager>(
        8, "CCSD DIIS",
        DIISManager::RemovalPolicy::LargestError,
        DIISManager::StoragePolicy::OnDisk
    );
    // No need to set sizes - libdiis auto-detects from first add_entry()
}

void CCEnergyWavefunction::diis_RHF(int iter) {
    dpdfile2 R1, T1_new;
    dpdbuf4 R2, T2_new;

    // Compute errors in-place (using existing DPD operations)
    compute_t1_error(&R1);  // R1 = T1_new - T1_old
    compute_t2_error(&R2);  // R2 = T2_new - T2_old

    // DIIS magic happens here
    ccsd_diis->add_entry(&R1, &R2, &T1_new, &T2_new);
    if (ccsd_diis->subspace_size() >= 2) {
        ccsd_diis->extrapolate(&T1_new, &T2_new);
    }
}

// Similar for ROHF and UHF, just with more amplitude components
```

**Code Reduction:**
- Before: 1,098 lines
- After: ~150-200 lines (setup + thin wrappers for each reference)
- **Savings: ~900 lines**

---

### 2. cclambda Module

**Current Structure:**
- `diis.cc` (857 lines) - Single file handling all references
  - Includes RHF, ROHF, and UHF cases
  - Takes `L_irr` parameter for excited state symmetry
  - Much more complex due to handling all references

**Unique Features:**
- Irrep parameter: `L_irr` for excited states
- Longer iteration checks and error handling
- More complex UHF case (~400 lines just for UHF)

**Migration Complexity:**
The cclambda DIIS is more complex because it handles excited states with different irreps. However, libdiis supports this through irrep-aware DPD buffers.

**Migration Strategy:**
```cpp
void CCLambdaWavefunction::diis(int iter, int L_irr) {
    // libdiis automatically handles irrep through DPD buffer structure
    dpdfile2 R1, L1_new;
    dpdbuf4 R2, L2_new;

    if (params.ref == 0) {  // RHF
        // Initialize DPD buffers with L_irr
        global_dpd_->file2_init(&R1, PSIF_CC_LAMBDA, L_irr, 0, 1, "R_IA");
        global_dpd_->file2_init(&L1_new, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
        global_dpd_->buf4_init(&R2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "R_IjAb");
        global_dpd_->buf4_init(&L2_new, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

        // Compute errors (existing logic)
        compute_lambda_errors(&R1, &R2, L_irr);

        // DIIS
        lambda_diis->add_entry(&R1, &R2, &L1_new, &L2_new);
        if (lambda_diis->subspace_size() >= 2) {
            lambda_diis->extrapolate(&L1_new, &L2_new);
        }
    }
    else if (params.ref == 1) {  // ROHF
        // Similar pattern with more components
    }
    else if (params.ref == 2) {  // UHF
        // Similar pattern with more components
    }
}
```

**Code Reduction:**
- Before: 857 lines
- After: ~200-250 lines (more complex due to irrep handling)
- **Savings: ~600-650 lines**

---

### 3. ccresponse Module

**Current Structure:**
- `diis.cc` (283 lines) - Response amplitudes
  - Takes perturbation label and frequency
  - Builds labels dynamically based on perturbation

**Unique Features:**
- Perturbation-specific: `const char *pert`
- Frequency-dependent: `double omega`
- Dynamic label construction: `sprintf(lbl, "X_%s_%d (%%d|%%d)", pert, irrep)`

**Migration Strategy:**
```cpp
void diis(int iter, const char *pert, int irrep, double omega) {
    dpdfile2 R1, X1;
    dpdbuf4 R2, X2;

    // Build labels (keep existing logic)
    char lbl[32];
    sprintf(lbl, "New X_%s_%d", pert, irrep);

    // Initialize with dynamic labels
    global_dpd_->file2_init(&X1, PSIF_CC_MISC, irrep, 0, 1, lbl);
    // ...

    // Compute errors
    compute_response_errors(&R1, &R2, pert, irrep, omega);

    // DIIS
    response_diis->add_entry(&R1, &R2, &X1, &X2);
    if (response_diis->subspace_size() >= 2) {
        response_diis->extrapolate(&X1, &X2);
    }
}
```

**Code Reduction:**
- Before: 283 lines
- After: ~80-100 lines
- **Savings: ~180-200 lines**

---

## Consolidation Feasibility Assessment

### ✅ Technical Compatibility: CONFIRMED

| Requirement | cc modules | libdiis | Compatible? |
|-------------|-----------|---------|-------------|
| DPD buffers (dpdbuf4) | ✅ Used extensively | ✅ Supported | ✅ YES |
| DPD files (dpdfile2) | ✅ Used extensively | ✅ Supported | ✅ YES |
| OnDisk storage | ✅ PSIF_CC_DIIS_* | ✅ PSIF_LIBDIIS | ✅ YES |
| Subspace management | ✅ Manual cycling | ✅ Automatic | ✅ YES |
| Error vector construction | ✅ Manual | ✅ Automatic | ✅ YES |
| B matrix construction | ✅ Manual | ✅ Automatic | ✅ YES |
| Linear system solving | ✅ Custom conditioning | ✅ NumPy with conditioning | ⚠️ DIFFERENT |

### ⚠️ Key Differences to Address

#### 1. Conditioning Methodology

**Current (cc modules):**
```cpp
void CCEnergyWavefunction::diis_invert_B(double** B, double* C, int dimension, double tolerance) {
    // Balanced diagonal preconditioning
    for (int i = 0; i < dimension-1; i++) {
        Sp[i] = pow(Bp[i][i], -1.0/2.0);  // Diagonal scaling
    }
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            Bp[i][j] *= Sp[i] * Sp[j];  // Apply conditioning
        }
    }
    B2->power(-1.0, tolerance);  // Pseudoinverse with tolerance
}
```

**libdiis approach:**
```python
def diis_coefficients(self):
    # Similar diagonal conditioning in Python
    diagonals = B.diagonal().copy()
    diagonals[-1] = 1
    if np.all(diagonals > 0):
        diagonals = diagonals ** (-0.5)
        B = np.einsum("i,ij,j -> ij", diagonals, B, diagonals)
        return np.linalg.lstsq(B, rhs, rcond=None)[0][:-1] * diagonals[:-1]
    else:
        return np.linalg.lstsq(B, rhs, rcond=None)[0][:-1]
```

**Assessment:** The approaches are **essentially equivalent**. Both use:
- Diagonal preconditioning
- Least-squares solver (LAPACK's DGESV vs NumPy's lstsq, which also uses LAPACK)
- Tolerance for ill-conditioning

**Risk:** Low - algorithms are mathematically equivalent

#### 2. Storage Files

**Current:** Uses `PSIF_CC_DIIS_ERR` and `PSIF_CC_DIIS_AMP`
**libdiis:** Uses `PSIF_LIBDIIS`

**Impact:** Cosmetic only - both use PSIO subsystem

#### 3. Fixed Subspace Size

**Current:** Hard-coded `nvector = 8`
**libdiis:** Configurable via constructor

**Migration:** Simply pass `8` to DIISManager constructor

#### 4. DPD Buffer Conversion

**Key Technical Detail:**

libdiis converts dpdbuf4/dpdfile2 to Matrix internally:
```python
# psi4/driver/procrouting/diis.py, line 124-125
elif isinstance(x, (core.dpdbuf4, core.dpdfile2)):
    copy = core.Matrix(x)  # <-- Conversion happens here
```

This means libdiis:
1. Converts DPD → Matrix for storage
2. Stores Matrix in PSIO
3. Can extrapolate back to DPD format via `axpy_matrix()`

The conversion is handled by Psi4's core Matrix constructor which accepts DPD buffers.

**Performance Impact:**
- Extra conversion overhead per DIIS iteration
- Likely negligible compared to CC iteration cost
- Could benchmark to verify

---

## Implementation Roadmap

### Phase 1: Proof of Concept (Low Risk)

**Target:** ccenergy RHF only
- Smallest scope (258 lines → ~30-40 lines)
- Most commonly used reference
- Easy to validate

**Steps:**
1. Create `ccenergy_diis_libdiis.cc` (new file, don't modify existing)
2. Implement `diis_RHF_libdiis()` using DIISManager
3. Add compile flag to switch between implementations
4. Run comprehensive test suite
5. Compare energies to machine precision
6. Benchmark performance

**Success Criteria:**
- All tests pass with identical energies
- Performance within 5% of original
- Code reduction achieved

**Estimated Effort:** 4-6 hours
**Risk:** Very Low (can easily revert)

### Phase 2: Expand to All ccenergy (Medium Risk)

**Target:** All ccenergy DIIS implementations
- diis_ROHF, diis_UHF

**Steps:**
1. Implement ROHF and UHF versions
2. Run full test suite for all references
3. Validate against regression tests
4. Monitor for any numerical differences

**Estimated Effort:** 6-8 hours
**Risk:** Low

### Phase 3: cclambda Migration (Higher Risk)

**Target:** cclambda DIIS (857 lines)

**Challenge:** More complex due to:
- Excited state irreps
- Larger codebase to refactor
- Less frequently used (harder to validate)

**Steps:**
1. Careful analysis of irrep handling
2. Create libdiis version alongside existing
3. Extensive testing with excited states
4. Compare results for multiple irreps
5. Performance validation

**Estimated Effort:** 12-16 hours
**Risk:** Medium

### Phase 4: ccresponse Migration (Lower Priority)

**Target:** ccresponse DIIS (283 lines)

**Challenge:**
- Perturbation-dependent
- Response properties are sensitive
- Requires careful validation

**Steps:**
1. Implement libdiis version
2. Validate against analytical derivatives
3. Test with multiple perturbations
4. Performance comparison

**Estimated Effort:** 6-10 hours
**Risk:** Medium

### Phase 5: Remove Legacy Code (After Validation)

**Only after extensive testing:**
1. Remove old DIIS implementations
2. Update documentation
3. Update developer guides
4. Add migration notes

---

## Risk Assessment

### Technical Risks

| Risk | Severity | Mitigation |
|------|----------|------------|
| Numerical differences | Medium | Extensive regression testing, bit-for-bit comparison |
| Performance regression | Low | Benchmark suite, profiling |
| DPD conversion overhead | Low | Profile-guided optimization if needed |
| Excited state handling | Medium | Thorough testing with multiple irreps |
| Response property errors | Medium | Analytical derivative validation |

### Implementation Risks

| Risk | Severity | Mitigation |
|------|----------|------------|
| Breaking existing code | High | Phased approach, compile flags for fallback |
| Incomplete testing | Medium | Comprehensive test suite expansion |
| User impact | Medium | Transparent migration, backward compatibility |
| Maintenance burden | Low | Better long-term maintainability |

### Mitigation Strategy

1. **Parallel Implementation:** Keep both versions during transition
2. **Compile-Time Switching:** `#ifdef USE_LIBDIIS_CC` for easy A/B testing
3. **Extensive Testing:**Run full regression suite for each module
4. **Performance Monitoring:** Benchmark critical calculations
5. **Gradual Rollout:** Phase 1 → 2 → 3 → 4, validate each phase

---

## Code Savings Estimate

| Module | Current Lines | Estimated After | Savings | Percentage |
|--------|---------------|-----------------|---------|------------|
| ccenergy | 1,098 | 150-200 | ~900 | 82% |
| cclambda | 857 | 200-250 | ~600 | 70% |
| ccresponse | 283 | 80-100 | ~180 | 64% |
| **TOTAL** | **2,238** | **430-550** | **~1,680** | **75%** |

**Conservative Estimate:** 1,500-1,700 lines of code reduction

---

## Recommendations

### Recommended Actions

#### 1. **High Priority: Implement Proof of Concept** ✅
   - Start with ccenergy RHF
   - Low risk, high reward
   - Validates feasibility
   - **Timeline:** 1 week
   - **Effort:** 4-6 hours

#### 2. **Medium Priority: Expand to All ccenergy** ✅
   - If POC successful
   - Consolidate all reference types
   - **Timeline:** 2-3 weeks after POC
   - **Effort:** 6-8 hours

#### 3. **Lower Priority: Migrate cclambda** ⚠️
   - More complex, requires careful testing
   - Significant code reduction benefit
   - **Timeline:** 1-2 months after ccenergy
   - **Effort:** 12-16 hours

#### 4. **Lower Priority: Migrate ccresponse** ⚠️
   - Smaller benefit, specialized use case
   - **Timeline:** After cclambda success
   - **Effort:** 6-10 hours

### Not Recommended (Yet)

❌ **Immediate full migration** - Too risky without validation
❌ **Big-bang approach** - Prefer phased implementation
❌ **Removing conditioning** - Keep similar numerical behavior

---

## Testing Strategy

### Required Tests

1. **Regression Tests:**
   - All existing cc module tests must pass
   - Energies must match to machine precision (< 1e-10 Eh)
   - Gradients must match (< 1e-8 Eh/bohr)
   - Properties must match (< 1e-6)

2. **New Tests:**
   - DIIS with 2, 4, 6, 8 vectors
   - Near-singular DIIS matrices
   - Excited states (multiple irreps)
   - Response properties (multiple perturbations)

3. **Performance Tests:**
   - Benchmark CC iterations before/after
   - Measure DIIS overhead separately
   - Profile memory usage
   - Compare I/O patterns

4. **Edge Cases:**
   - Calculations that previously failed DIIS
   - Very small molecules (DIIS dominates)
   - Very large molecules (I/O critical)
   - Difficult convergence cases

---

## Success Metrics

### Must Have
- ✅ All regression tests pass
- ✅ Bit-for-bit identical energies (single-point)
- ✅ Gradients within 1e-8 Eh/bohr
- ✅ No performance regression > 5%
- ✅ Code reduction > 1,000 lines

### Nice to Have
- Performance improvement from better I/O
- Easier maintenance
- Unified DIIS behavior across all modules
- Better documentation

---

## Conclusion

**The migration of cc modules to libdiis is TECHNICALLY FEASIBLE and RECOMMENDED** with proper phased implementation and testing.

### Summary

**Pros:**
- ✅ Eliminates ~1,680 lines of duplicate code
- ✅ Technical compatibility confirmed
- ✅ libdiis supports all required data types
- ✅ Improved maintainability
- ✅ Consistent DIIS behavior across Psi4
- ✅ Proven approach (occ/dfocc/dct already use it)

**Cons:**
- ⚠️ Requires careful refactoring
- ⚠️ Medium implementation risk
- ⚠️ Extensive testing needed
- ⚠️ User-facing changes (though transparent)

**Recommendation:** **PROCEED with phased approach**

1. Start with ccenergy RHF proof of concept
2. Validate thoroughly
3. Expand incrementally if successful
4. Maintain backward compatibility during transition

**Overall Assessment:** This is a worthwhile refactoring effort that will significantly reduce code duplication and improve long-term maintainability, with manageable risk through careful phased implementation.

---

## Next Steps

1. ✅ Present this analysis for review
2. ⏳ Get stakeholder approval for POC
3. ⏳ Implement ccenergy RHF POC
4. ⏳ Run validation tests
5. ⏳ Decision point: proceed or abort
6. ⏳ If successful, continue with phases 2-4

---

## Files Referenced

### cc/ccenergy
- `psi4/src/psi4/cc/ccenergy/diis.cc`
- `psi4/src/psi4/cc/ccenergy/diis_RHF.cc`
- `psi4/src/psi4/cc/ccenergy/diis_ROHF.cc`
- `psi4/src/psi4/cc/ccenergy/diis_UHF.cc`

### cc/cclambda
- `psi4/src/psi4/cc/cclambda/diis.cc`

### cc/ccresponse
- `psi4/src/psi4/cc/ccresponse/diis.cc`

### libdiis
- `psi4/src/psi4/libdiis/diismanager.h`
- `psi4/src/psi4/libdiis/diismanager.cc`
- `psi4/driver/procrouting/diis.py`

---

*Analysis completed: 2025-11-18*
*Analyst: Claude (Sonnet 4.5)*
*Status: Consolidation feasible with phased approach*
