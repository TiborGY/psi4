# OCC Module Refactoring Plan - Integral Permutation Utilities

## Executive Summary

**Objective**: Refactor the OCC (Orbital-Optimized Coupled Cluster) module to use the integral permutation utilities from libtrans, improving code readability and maintainability.

**Scope**: 244 buf4_sort calls across 54 source files
**Addressable**: 207 calls (85%) use common permutations supported by utilities
**Timeline**: 3-4 working days
**Risk Level**: Low (proven pattern from DCT refactoring)

---

## Module Overview

The OCC module implements orbital-optimized methods including:
- OMP2, OMP3, OCEPA (orbital-optimized MP2/MP3/CEPA)
- Gradients and response properties
- Ionization potential calculations

**API Compatibility**: ✅ 100% - Uses ID() macro style identical to DCT

---

## Impact Analysis

### File-Level Statistics

**Top 20 Files by buf4_sort Count:**

| File | Total | prqs | Priority | Category |
|------|-------|------|----------|----------|
| gfock_ea.cc | 38 | 8 | High | Green's function |
| trans_ints_uhf.cc | 24 | 21 | **⭐ Critical** | Integral transform |
| t2_2nd_general.cc | 23 | 9 | High | Amplitude |
| t2_2nd_sc.cc | 17 | 4 | Medium | Amplitude |
| t2_amps_remp.cc | 12 | 4 | Medium | Amplitude |
| t2_amps.cc | 12 | 4 | Medium | Amplitude |
| corr_tpdm.cc | 12 | 3 | Medium | Density |
| kappa_orb_resp_iter.cc | 11 | 0 | Medium | Response |
| z_vector.cc | 10 | 0 | Medium | Response |
| set_t2_amplitudes_mp2.cc | 10 | 5 | Medium | Amplitude |
| kappa_orb_resp.cc | 10 | 0 | Medium | Response |
| cepa_iterations.cc | 10 | 5 | Medium | CEPA |
| trans_ints_rhf.cc | 6 | **6** | **⭐ Critical** | Integral transform |
| trans_ints_rmp2.cc | 4 | 1 | Low | Integral transform |
| trans_ints_ump2.cc | 4 | 3 | Low | Integral transform |

**Key Observations:**
- `trans_ints_rhf.cc`: 100% prqs (perfect starting point)
- `trans_ints_uhf.cc`: 88% prqs (high impact)
- Integral transformation files are cleanest targets

### Permutation Type Distribution

| Permutation | Count | Utility Function | Complexity |
|-------------|-------|------------------|------------|
| **prqs** | 74 | `chemist_to_physicist()` | Low |
| **rspq** | 62 | `transpose_bra_ket()` | Low |
| **psrq** | 26 | `apply_permutation()` | Medium |
| **qpsr** | 22 | `apply_permutation()` | Medium |
| **pqsr** | 11 | `apply_permutation()` | Medium |
| **qprs** | 9 | `apply_permutation()` | Medium |
| **prsq** | 3 | `apply_permutation()` | Medium |
| **Others** | 37 | Various patterns | Variable |

**Coverage**: 207 of 244 calls (85%) are addressable with current utilities

---

## Implementation Strategy

### Phase 1: Integral Transformation Files ⭐ **START HERE**

**Priority**: Critical
**Effort**: 0.5-1 day
**Risk**: Minimal

**Files**:
1. `trans_ints_rhf.cc` (6 calls, all prqs)
2. `trans_ints_uhf.cc` (24 calls, 21 prqs)
3. `trans_ints_rmp2.cc` (4 calls, 1 prqs)
4. `trans_ints_ump2.cc` (4 calls, 3 prqs)

**Rationale**:
- Cleanest, most straightforward code
- Pure chemist-to-physicist transformations
- Well-commented with timer blocks
- Minimal conditional logic
- Serves as template for other files

**Example Refactoring** (trans_ints_rhf.cc:102-103):

```cpp
// BEFORE
global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                       ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[O,O]"),
                       "MO Ints <OO|OO>");
global_dpd_->buf4_close(&K);

// AFTER
global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                       ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
libtrans::IntegralPermutations::chemist_to_physicist(&K, PSIF_LIBTRANS_DPD,
                       ID("[O,O]"), ID("[O,O]"), "MO Ints <OO|OO>");
global_dpd_->buf4_close(&K);
```

**Testing**:
- OMP2 energy calculation
- OMP3 energy calculation
- Verify energies match to 10^-10 Eh

---

### Phase 2: Amplitude Manipulation Files

**Priority**: High
**Effort**: 1-1.5 days
**Risk**: Low

**Files**:
1. `t2_amps.cc` (12 calls, 4 prqs)
2. `t2_amps_remp.cc` (12 calls, 4 prqs)
3. `set_t2_amplitudes_mp2.cc` (10 calls, 5 prqs)
4. `cepa_iterations.cc` (10 calls, 5 prqs)
5. `t2_2nd_general.cc` (23 calls, 9 prqs)
6. `t2_2nd_sc.cc` (17 calls, 4 prqs)

**Characteristics**:
- Mix of prqs, rspq, and other permutations
- Amplitude sorting and manipulation
- Some conditional logic based on method type

**Common Pattern**:
```cpp
// Amplitude sorting (cepa_iterations.cc:309)
global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,V]"), ID("[O,V]"),
                       "T2 (OV|OV)");

// With utility
libtrans::IntegralPermutations::chemist_to_physicist(&T, PSIF_OCC_DPD,
                       ID("[O,V]"), ID("[O,V]"), "T2 (OV|OV)");
```

**Testing**:
- CEPA(0) energy calculation
- OMP2 T2 amplitudes verification
- Gradient calculations

---

### Phase 3: Density and Response Files

**Priority**: Medium
**Effort**: 1-1.5 days
**Risk**: Low-Medium

**Files**:
1. `corr_tpdm.cc` (12 calls, 3 prqs)
2. `coord_grad.cc` (moderate complexity)
3. `kappa_orb_resp.cc` (10 calls, 0 prqs)
4. `kappa_orb_resp_iter.cc` (11 calls, 0 prqs)
5. `z_vector.cc` (10 calls, 0 prqs)

**Characteristics**:
- Two-particle density matrices
- Orbital response calculations
- Mix of permutation types
- More complex logic

**Special Considerations**:
- TPDM files use various permutations (rspq, psrq, qpsr)
- Need to use `transpose_bra_ket()` and `apply_permutation()`
- Response files may have fewer prqs but other permutations

**Example** (corr_tpdm.cc:736):
```cpp
// BEFORE
global_dpd_->buf4_sort(&G, PSIF_OCC_DENSITY, rspq, ID("[V,o]"), ID("[O,v]"),
                       "TPDM <Vo|Ov>");

// AFTER
libtrans::IntegralPermutations::transpose_bra_ket(&G, PSIF_OCC_DENSITY,
                       ID("[V,o]"), ID("[O,v]"), "TPDM <Vo|Ov>");
```

**Testing**:
- Analytic gradient calculations
- TPDM trace verification
- Orbital response convergence

---

### Phase 4: Green's Function and Specialized Files

**Priority**: Medium-Low
**Effort**: 0.5-1 day
**Risk**: Medium

**Files**:
1. `gfock_ea.cc` (38 calls, 8 prqs) - Largest file
2. `gfock_diag.cc` (4 calls)
3. `ep2_ip.cc` (5 calls)
4. `omp2_ip_poles.cc` (5 calls)
5. `omp3_ip_poles.cc` (5 calls)

**Characteristics**:
- Green's function calculations
- Ionization potential methods
- More specialized code
- Mix of permutations

**Special Considerations**:
- `gfock_ea.cc` has 38 calls but only 8 prqs
- Many use other permutation types
- Need careful testing with EKT calculations

**Testing**:
- EKT-OMP2 calculations
- IP-OMP2/OMP3 calculations
- Verify poles and residues

---

## Detailed Implementation Checklist

### Pre-Implementation

- [x] Analyze module structure
- [x] Identify file priorities
- [x] Create implementation plan
- [ ] Review DCT refactoring patterns
- [ ] Set up test suite baseline

### Phase 1: Integral Transforms (Days 1-2)

#### trans_ints_rhf.cc
- [ ] Add include: `#include "psi4/libtrans/integral_permutations.h"`
- [ ] Refactor line 102: (OO|OO) → <OO|OO>
- [ ] Refactor line 110: (OO|OV) → <OO|OV>
- [ ] Refactor line 118: (OV|OV) → <OO|VV>
- [ ] Refactor line 126: (OO|VV) → <OV|OV>
- [ ] Refactor line 139: (OV|VV) → <OV|VV> (conditional)
- [ ] Refactor line 149: (VV|VV) → <VV|VV> (conditional)
- [ ] Test: OMP2/He energy
- [ ] Test: OMP3/H2O energy
- [ ] Commit: "Refactor trans_ints_rhf.cc to use integral utilities"

#### trans_ints_uhf.cc
- [ ] Add include
- [ ] Refactor 21 prqs calls (lines TBD after file read)
- [ ] Test: UMP2 energy
- [ ] Test: Verify spin-separated integrals
- [ ] Commit: "Refactor trans_ints_uhf.cc to use integral utilities"

#### trans_ints_rmp2.cc & trans_ints_ump2.cc
- [ ] Add includes
- [ ] Refactor remaining prqs calls
- [ ] Test: RMP2 and UMP2 energies
- [ ] Commit: "Refactor MP2 integral transforms to use utilities"

### Phase 2: Amplitudes (Days 2-3)

#### set_t2_amplitudes_mp2.cc
- [ ] Add include
- [ ] Refactor 5 prqs calls
- [ ] Test: T2 amplitude initialization
- [ ] Commit: "Refactor T2 amplitude initialization"

#### cepa_iterations.cc
- [ ] Add include
- [ ] Refactor 5 prqs calls
- [ ] Handle rspq permutations with transpose_bra_ket()
- [ ] Test: CEPA(0) convergence
- [ ] Commit: "Refactor CEPA iterations to use utilities"

#### t2_amps.cc & t2_amps_remp.cc
- [ ] Add includes
- [ ] Refactor prqs calls
- [ ] Handle mixed permutations
- [ ] Test: Amplitude updates
- [ ] Commit: "Refactor T2 amplitude updates"

#### t2_2nd_general.cc & t2_2nd_sc.cc
- [ ] Add includes
- [ ] Refactor 9 prqs in general, 4 in sc
- [ ] Handle complex permutation patterns
- [ ] Test: Second-order corrections
- [ ] Commit: "Refactor second-order T2 updates"

### Phase 3: Density/Response (Days 3-4)

#### corr_tpdm.cc
- [ ] Add include
- [ ] Refactor 3 prqs calls
- [ ] Handle rspq with transpose_bra_ket()
- [ ] Handle other permutations with apply_permutation()
- [ ] Test: TPDM trace
- [ ] Commit: "Refactor TPDM construction"

#### Response files (kappa_orb_resp*.cc, z_vector.cc)
- [ ] Add includes
- [ ] Refactor non-prqs permutations
- [ ] Use apply_permutation() for custom permutations
- [ ] Test: Orbital response convergence
- [ ] Commit: "Refactor orbital response to use utilities"

### Phase 4: Specialized (Day 4)

#### gfock_ea.cc
- [ ] Add include
- [ ] Refactor 8 prqs calls
- [ ] Handle 30 other permutations selectively
- [ ] Test: EKT calculations
- [ ] Commit: "Refactor Green's function calculations"

#### IP/EA methods
- [ ] Add includes to ep2_ip.cc, omp*_ip_poles.cc
- [ ] Refactor remaining calls
- [ ] Test: IP-OMP2/OMP3
- [ ] Commit: "Refactor IP/EA methods"

### Post-Implementation

- [ ] Run full OCC test suite
- [ ] Verify all energies match baseline
- [ ] Check gradients match baseline
- [ ] Update module documentation
- [ ] Create summary commit
- [ ] Update IMPLEMENTATION_SUMMARY.md

---

## Testing Strategy

### Baseline Establishment

**Before ANY refactoring**, capture baseline results:

```bash
# Run comprehensive OCC test suite
pytest tests/pytests/test_occ.py -v --capture=no

# Capture specific test outputs
psi4 tests/omp2-1/test.in > baseline_omp2.out
psi4 tests/omp3-1/test.in > baseline_omp3.out
psi4 tests/ocepa1/test.in > baseline_ocepa.out
psi4 tests/omp2-grad1/test.in > baseline_grad.out
```

### Progressive Testing

**After each file refactoring:**

1. **Unit-level**: Re-run tests specific to modified functionality
2. **Integration-level**: Run full method calculations
3. **Regression**: Compare energies/gradients to baseline

### Test Hierarchy

**Level 1: Quick Smoke Tests (Run after each file)**
- He atom OMP2 energy (fast, ~1 second)
- H2O OMP2 energy (medium, ~5 seconds)
- Verify energies match to 10^-10 Eh

**Level 2: Method-Specific Tests (Run after each phase)**
- OMP2, OMP3, OCEPA energies
- RHF and UHF references
- Various basis sets (STO-3G, cc-pVDZ)

**Level 3: Comprehensive Tests (Run at end)**
- Analytic gradients
- TPDM properties
- Response properties
- IP-OMP2/OMP3 calculations
- EKT calculations

**Level 4: Full Regression (Final validation)**
- Entire OCC test suite
- All combinations of methods/references
- Timing benchmarks (ensure no performance regression)

### Validation Criteria

**Energy Accuracy**:
- Must match baseline to 10^-10 Eh (within numerical precision)
- No exceptions for any test

**Gradient Accuracy**:
- Must match baseline to 10^-8 Eh/Bohr
- RMS and maximum differences both checked

**Integral Accuracy**:
- Spot-check integral values in debug mode
- Verify antisymmetry properties maintained

**Performance**:
- No more than 2% slowdown (utilities are inline-eligible)
- Actually expect negligible difference

### Test Automation

Create test script: `tests/occ_refactor_validation.sh`

```bash
#!/bin/bash
# OCC Refactoring Validation Script

echo "Running OCC refactoring validation..."

# Phase 1: Quick tests
echo "Phase 1: Quick validation"
psi4 tests/omp2-1/test.in
if [ $? -ne 0 ]; then
    echo "FAIL: OMP2 basic test failed"
    exit 1
fi

# Phase 2: Method coverage
echo "Phase 2: Method coverage"
for test in omp2-1 omp2-2 omp3-1 omp3-2 ocepa1; do
    psi4 tests/$test/test.in > /tmp/$test.out
    # Compare energies
    grep "Total Energy" /tmp/$test.out > /tmp/$test.energy
    diff -q /tmp/$test.energy baseline/$test.energy
    if [ $? -ne 0 ]; then
        echo "FAIL: $test energy mismatch"
        exit 1
    fi
done

# Phase 3: Gradients
echo "Phase 3: Gradient tests"
psi4 tests/omp2-grad1/test.in
psi4 tests/omp3-grad1/test.in

echo "All validations passed!"
```

---

## Risk Mitigation

### Risk 1: Numerical Precision Issues
**Likelihood**: Very Low
**Impact**: Medium
**Mitigation**:
- Utilities are direct wrappers (no algorithmic changes)
- DPD library handles all numerics
- Test thoroughly with tight thresholds

### Risk 2: Conditional Logic Errors
**Likelihood**: Low
**Impact**: Medium
**Mitigation**:
- Carefully handle conditional buf4_sort calls
- Test all code paths (RHF/UHF, different methods)
- Preserve all conditionals exactly

### Risk 3: Label/Index Mismatches
**Likelihood**: Low
**Impact**: Medium
**Mitigation**:
- Use exact same labels as original code
- Use exact same indices as original code
- Verify with debug prints if uncertain

### Risk 4: Performance Regression
**Likelihood**: Very Low
**Impact**: Low
**Mitigation**:
- Utilities are inline-eligible
- No additional function call overhead expected
- Benchmark before/after

### Risk 5: Incomplete Testing
**Likelihood**: Medium
**Impact**: High
**Mitigation**:
- Follow phased testing approach
- Don't move to next file until current validated
- Maintain baseline results throughout

---

## Code Review Guidelines

### Before Committing Each File:

1. **Visual Inspection**:
   - [ ] All buf4_sort calls properly replaced
   - [ ] Correct utility function chosen
   - [ ] Indices match exactly
   - [ ] Labels match exactly
   - [ ] Include statement added
   - [ ] No unintended changes

2. **Pattern Verification**:
   ```bash
   # Ensure no naked prqs remain
   grep "buf4_sort.*prqs" <file>
   # Should return 0 matches for fully refactored file
   ```

3. **Testing**:
   - [ ] Quick smoke test passes
   - [ ] Relevant method tests pass
   - [ ] Energy/gradient accuracy verified

4. **Documentation**:
   - [ ] Commit message describes changes
   - [ ] Notes any special cases or decisions
   - [ ] References issue/plan document

---

## Success Metrics

### Quantitative Goals

| Metric | Target | Measurement |
|--------|--------|-------------|
| buf4_sort calls refactored | 200+ | Count remaining buf4_sort |
| prqs permutations eliminated | 70+ | grep count |
| Files updated | 15-20 | git diff --numstat |
| Test pass rate | 100% | pytest results |
| Energy accuracy | 10^-10 Eh | Baseline comparison |
| Gradient accuracy | 10^-8 Eh/Bohr | Baseline comparison |
| Performance change | < 2% | Timing benchmarks |

### Qualitative Goals

- [ ] Code is more readable (semantic function names vs cryptic codes)
- [ ] Pattern consistency between DCT and OCC modules
- [ ] Easier maintenance for future developers
- [ ] Documentation updated appropriately

---

## Timeline Estimate

### Optimistic (3 days)
- Day 1: Phase 1 complete, Phase 2 started
- Day 2: Phase 2 complete, Phase 3 started
- Day 3: Phase 3 & 4 complete, final testing

### Realistic (4 days)
- Day 1: Phase 1 complete
- Day 2: Phase 2 complete
- Day 3: Phase 3 complete
- Day 4: Phase 4 complete, comprehensive testing

### Conservative (5 days)
- Days 1-2: Phase 1 & 2 with extra testing
- Day 3: Phase 3
- Day 4: Phase 4
- Day 5: Full regression testing and cleanup

**Recommended**: Plan for 4 days with buffer

---

## Future Enhancements

After successful OCC refactoring, consider:

1. **Additional Utility Functions**: If OCC reveals common patterns not yet in utilities
2. **cctransort Module**: Natural next target (88 buf4_sort calls)
3. **CC Family Extension**: Requires numeric index support (265+ calls)
4. **Performance Profiling**: Detailed timing comparison
5. **Documentation**: Tutorial for refactoring other modules

---

## Conclusion

The OCC module refactoring represents a significant expansion of the integral permutation utilities' adoption. With 244 buf4_sort calls and 100% API compatibility, it offers:

✅ **High Impact**: Nearly doubles utility usage (from 55 to ~300 calls)
✅ **Low Risk**: Proven pattern from DCT refactoring
✅ **Clear Value**: Improved readability and maintainability
✅ **Good Timeline**: 3-4 days of focused work

The phased approach, comprehensive testing strategy, and risk mitigation plan ensure a successful refactoring with minimal disruption.

**Recommendation**: Proceed with Phase 1 implementation starting with `trans_ints_rhf.cc`.
