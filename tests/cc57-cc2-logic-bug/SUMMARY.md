# Test Case for Boolean Logic Bug in t2.cc:48

## Executive Summary

I've created a comprehensive test case and demonstration for the critical boolean logic bug found in `psi4/src/psi4/cc/ccenergy/t2.cc` at line 48.

**Location:** `/home/user/psi4/tests/cc57-cc2-logic-bug/`

## The Bug

**File:** `psi4/src/psi4/cc/ccenergy/t2.cc:48`

```cpp
// BUGGY:
if (params_.wfn != "CC2" || params_.wfn != "EOM_CC2") { /* skip all this is wfn=CC2 */

// CORRECT:
if (params_.wfn != "CC2" && params_.wfn != "EOM_CC2") { /* skip all this if wfn=CC2 */
```

The bug uses `||` (OR) instead of `&&` (AND), making the condition **always evaluate to true**, even for CC2/EOM_CC2 when it should be false.

## Why It's a Bug

Using `||` (OR):
- If `wfn == "CC2"`: `wfn != "EOM_CC2"` is TRUE → condition is TRUE ❌
- If `wfn == "EOM_CC2"`: `wfn != "CC2"` is TRUE → condition is TRUE ❌
- If `wfn == "CCSD"`: both are TRUE → condition is TRUE ✓

Using `&&` (AND):
- If `wfn == "CC2"`: condition is FALSE ✓ (correctly skips CCSD terms)
- If `wfn == "EOM_CC2"`: condition is FALSE ✓ (correctly skips CCSD terms)
- If `wfn == "CCSD"`: condition is TRUE ✓ (correctly executes CCSD terms)

## Impact

### What the Bug Would Do (if triggered)

If `t2_build()` were called for CC2/EOM_CC2, the buggy logic would cause:

1. **Wrong algorithm:** Execute CCSD-specific intermediates (FaetT2, FmitT2, WmnijT2, etc.) instead of CC2 approximations
2. **Broken physics:** CC2 uses T1-only intermediates (O(N⁵)), but buggy code would mix in T2-based CCSD terms (O(N⁶))
3. **Incorrect results:** Energy would be neither proper CC2 nor proper CCSD

### Why Normal CC2 Calculations Are Safe

**Good news:** Most CC2 calculations are NOT affected because:

- Main iteration loop (`ccenergy.cc:241`) calls `cc2_t2_build()` for CC2 ✓
- `cc2_t2_build()` is the correct implementation without the bug ✓
- `t2_build()` is only called for CCSD/CC3/other methods ✓

**Potential issue:** Edge cases where `t2_build()` might be called for CC2:
- `one_step()` function with `just_residuals` option (`ccenergy.cc:564`)
- Any future code that inadvertently calls `t2_build()` for CC2

## Test Files Created

### 1. `input.dat`
- Main test file for pytest/ctest integration
- Tests CC2 and CCSD on H2O with cc-pVDZ
- Verifies CC2 > CCSD (CC2 recovers less correlation)
- Uses reference values from existing test (cc36)
- **This test currently PASSES** because normal CC2 uses `cc2_t2_build()`

### 2. `demonstrate_bug.py`
- Python script that demonstrates the boolean logic error
- Shows exactly why the condition is always true
- Tests multiple wavefunction types
- **Run this to see the bug clearly**

### 3. `README.md`
- Comprehensive documentation of the bug
- Theoretical background on CC2 vs CCSD
- Code references and related files
- Fix instructions

### 4. `SUMMARY.md`
- This file - overview and usage guide

### 5. Supporting files
- `CMakeLists.txt` - Build system integration
- `test_input.py` - pytest wrapper

## How to Use This Test

### Option 1: See the Bug Logic (Recommended First)

```bash
cd /home/user/psi4/tests/cc57-cc2-logic-bug
python3 demonstrate_bug.py
```

This clearly shows why the OR logic is always true.

### Option 2: Run the Regression Test

```bash
cd /home/user/psi4/tests/cc57-cc2-logic-bug
psi4 input.dat
```

**Note:** This test will likely PASS even with the bug, because normal CC2 calculations don't trigger the buggy code path. The test serves as:
- Documentation of correct behavior
- Regression test after the fix
- Verification that CC2 < CCSD (correct physics)

### Option 3: Truly Expose the Bug (Advanced)

To truly expose the bug, you would need to:

1. Modify `ccenergy.cc:241` to call `t2_build()` instead of `cc2_t2_build()` for CC2
2. Run the test - it should fail or give wrong energies
3. Fix the bug in `t2.cc:48` (change `||` to `&&`)
4. Run again - should now pass

## Demonstration Output

```
Method       | Expected | Buggy OR  | Fixed AND  | Status
----------------------------------------------------------------------
CC2          | False    | True      | False      | ❌ BUG → ✓ FIXED
EOM_CC2      | False    | True      | False      | ❌ BUG → ✓ FIXED
CCSD         | True     | True      | True       | ✓ OK
CC3          | True     | True      | True       | ✓ OK
```

**The buggy OR logic is ALWAYS TRUE!**

## Recommendations

1. **Fix the bug** in `t2.cc:48` by changing `||` to `&&`
2. **Keep this test** as a regression test and documentation
3. **Review `one_step()` function** to ensure it handles CC2 correctly
4. **Consider adding explicit checks** in `t2_build()` to error if called for CC2

## Fix Required

**File:** `psi4/src/psi4/cc/ccenergy/t2.cc`

**Line:** 48

**Change:**
```cpp
// Before:
if (params_.wfn != "CC2" || params_.wfn != "EOM_CC2") { /* skip all this is wfn=CC2 */

// After:
if (params_.wfn != "CC2" && params_.wfn != "EOM_CC2") { /* skip all this if wfn=CC2 */
```

**Alternative (more explicit):**
```cpp
// More readable version:
if (params_.wfn != "CC2" && params_.wfn != "EOM_CC2") {
    // Execute CCSD-specific terms for all methods except CC2/EOM_CC2
    FaetT2();
    FmitT2();
    // ...
}
```

**Or using positive logic:**
```cpp
// Most readable version:
bool is_cc2_type = (params_.wfn == "CC2" || params_.wfn == "EOM_CC2");
if (!is_cc2_type) {
    // Execute CCSD-specific terms
    FaetT2();
    FmitT2();
    // ...
}
```

## References

- Original issue: Code review of ccenergy files
- Related test: `tests/cc36/` (CC2 reference)
- CC2 theory: Christiansen et al., Chem. Phys. Lett. 243, 409 (1995)
- CCSD theory: Purvis & Bartlett, J. Chem. Phys. 76, 1910 (1982)

## Theory Background

**CC2 (Second-order Coupled Cluster):**
- Approximation to CCSD
- T2 amplitudes computed to first order in fluctuation potential
- Intermediates use **only T1 amplitudes**
- Scaling: O(N⁵) - same as MP2
- Less accurate but faster than CCSD

**CCSD (Coupled Cluster Singles and Doubles):**
- Full treatment of singles and doubles
- Intermediates use **both T1 and T2 amplitudes**
- Scaling: O(N⁶)
- More accurate, more expensive

**Expected:** CC2 energy > CCSD energy (CC2 recovers less correlation)

**With bug:** If triggered, CC2 would incorrectly mix CCSD terms and give wrong results
