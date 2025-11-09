# CC2 Logic Bug Regression Test

## Purpose

This test serves as a regression test and documentation for a boolean logic bug in `psi4/src/psi4/cc/ccenergy/t2.cc` at line 48.

## The Bug

**Location:** `psi4/src/psi4/cc/ccenergy/t2.cc:48`

**Buggy Code:**
```cpp
if (params_.wfn != "CC2" || params_.wfn != "EOM_CC2") { /* skip all this is wfn=CC2 */
    FaetT2();
    FmitT2();
    WmnijT2();
    // ... more CCSD-specific terms
}
```

**Problem:**

The condition uses `||` (OR) instead of `&&` (AND), making it **always evaluate to true**:

| wfn value | wfn != "CC2" | wfn != "EOM_CC2" | Result with OR | Expected Result |
|-----------|--------------|------------------|----------------|-----------------|
| "CC2"     | false        | true             | **true** ❌     | false ✓         |
| "EOM_CC2" | true         | false            | **true** ❌     | false ✓         |
| "CCSD"    | true         | true             | true ✓         | true ✓          |

**Correct Code:**
```cpp
if (params_.wfn != "CC2" && params_.wfn != "EOM_CC2") { /* skip all this if wfn=CC2 */
```

## Impact

### Theoretical Impact (if triggered)

If `t2_build()` were called for CC2 or EOM_CC2, the buggy logic would cause:

1. **Incorrect algorithm:** CC2 would execute CCSD-specific intermediate terms (FaetT2, FmitT2, WmnijT2, etc.) instead of the proper CC2 approximations (cc2_faeT2, cc2_fmiT2, etc.)

2. **Wrong physics:**
   - CC2 uses **T1-only intermediates** for O(N⁵) scaling
   - CCSD uses **T1+T2 intermediates** for O(N⁶) scaling
   - The bug would create an inconsistent hybrid

3. **Incorrect energies:** Results would be neither proper CC2 nor proper CCSD

### Actual Impact (current code)

**Good news:** Normal CC2 calculations are **NOT affected** by this bug because:

- The main iteration loop (`ccenergy.cc:241`) correctly calls `cc2_t2_build()` for CC2
- `cc2_t2_build()` does not contain this bug
- `t2_build()` is only called for CCSD, CC3, and other non-CC2 methods

**Potential issue:** The bug could manifest in edge cases:

- `one_step()` function (`ccenergy.cc:564`) calls `t2_build()` without checking the method type
- If used with `just_residuals` option for CC2, incorrect results could occur
- Any future code that inadvertently calls `t2_build()` for CC2/EOM_CC2

## Test Description

This test:

1. Runs CC2 calculation on H2O with cc-pVDZ basis
2. Runs CCSD calculation for comparison
3. Verifies that CC2 gives higher energy than CCSD (as expected from theory)
4. Uses reference values from existing test (cc36)

### Expected Behavior

- **CC2 energy:** -76.22960579193078 Ha
- **CCSD energy:** -76.24046955364111 Ha
- **Difference:** CC2 > CCSD (CC2 recovers less correlation)

### Test Success Criteria

1. CC2 energy matches reference
2. CCSD energy matches reference
3. CC2 energy > CCSD energy (assertion at line 73)

If the assertion fails (CC2 ≤ CCSD), it may indicate:
- The bug is affecting results
- The test references need updating
- An implementation problem

## Files

- `input.dat` - Test input with detailed bug documentation
- `test_input.py` - pytest integration
- `CMakeLists.txt` - Build system integration
- `README.md` - This file

## References

- Original CC2 test: `tests/cc36/`
- CC2 implementation: `psi4/src/psi4/cc/ccenergy/cc2_t2.cc`
- CCSD implementation: `psi4/src/psi4/cc/ccenergy/t2.cc`
- CC2 theory: Christiansen et al., Chem. Phys. Lett. 243, 409 (1995)

## Related Code

### Correct CC2 Implementation
- `cc2_t2_build()` in `cc2_t2.cc:43` ✓
- Uses `cc2_faeT2()`, `cc2_fmiT2()` - T1-only intermediates
- O(N⁵) scaling

### Buggy Code (not normally executed for CC2)
- `t2_build()` in `t2.cc:44` with bug at line 48 ❌
- Would incorrectly use `FaetT2()`, `FmitT2()`, etc. for CC2
- Would break CC2 approximation

## Fix Required

Change line 48 in `psi4/src/psi4/cc/ccenergy/t2.cc` from:
```cpp
if (params_.wfn != "CC2" || params_.wfn != "EOM_CC2") {
```

To:
```cpp
if (params_.wfn != "CC2" && params_.wfn != "EOM_CC2") {
```

This test will serve as a regression test after the fix is applied.
