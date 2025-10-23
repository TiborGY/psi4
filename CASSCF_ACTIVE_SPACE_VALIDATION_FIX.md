# CASSCF Active Space Validation Fix

## Issue Summary

CASSCF calculations would fail with a cryptic "DSYEV diagonalizer failed" error when users specified active orbitals in symmetry irreps that had no (or insufficient) available molecular orbitals.

## Symptom

### Error Message (Before Fix)
```
RuntimeError:
Fatal Error: DSYEV diagonalizer failed in DETCI SEM!
Error occurred in file: /path/to/psi4/src/psi4/detci/sem.cc on line: 561
```

### Example Failing Input
```python
molecule {
H
H 1 0.60
}

set basis cc-pVDZ
set docc [ 1, 0, 0, 0, 0, 0, 0, 0 ]
set active [1, 0, 0, 0, 0, 1, 0, 0 ]  # Invalid: requests orbital in irrep 5
energy('casscf')
```

In this case, the H₂ molecule in D2h point group has no molecular orbitals in irrep 5, but the user requested 1 active orbital there (`active[5] = 1`).

## Root Cause Analysis

### Primary Cause

**Missing validation in `ras_set3()` function** ([psi4/src/psi4/libqt/ras_set.cc](psi4/src/psi4/libqt/ras_set.cc))

The code validated that:
- ✅ Array sizes matched the number of irreps
- ✅ Total orbitals didn't exceed available orbitals (aggregate check)

But it did NOT validate that:
- ❌ Individual `ACTIVE[irrep]` values were ≤ `orbspi[irrep]` (available orbitals per irrep)
- ❌ Individual `RAS1[irrep]`, `RAS2[irrep]`, `RAS3[irrep]`, `RAS4[irrep]` values were valid

### Cascading Failure Mechanism

1. **Invalid Input Accepted**
   ```cpp
   options.fill_int_array("ACTIVE", ras_opi[1]);  // Line 297
   // Sets ras_opi[1][5] = 1, but orbspi[5] = 0
   // No validation occurs!
   ```

2. **Negative Orbital Count Calculated** (Lines 314-315)
   ```cpp
   restruocc[irrep] = orbspi[irrep] - frdocc[irrep] - restrdocc[irrep] -
                      ras_opi[0][irrep] - ras_opi[1][irrep] - fruocc[irrep];
   // For irrep 5: restruocc[5] = 0 - 0 - 0 - 0 - 1 - 0 = -1
   ```

3. **Orbital Reordering Corruption** (Lines 399-413)
   ```cpp
   while (tras[i][irrep]) {  // -1 is non-zero, enters loop!
       point = used[irrep] + offset[irrep];
       order[point] = cnt++;
       used[irrep]++;
       tras[i][irrep]--;  // Becomes -2, -3, -4... potentially infinite!
   }
   ```

   The negative value caused either:
   - An "Invalid point value" exception (if bounds checking caught it), OR
   - Corrupted orbital configuration data

4. **CI Hamiltonian Matrix Corruption**
   - Invalid orbital configuration led to incorrect determinant counting
   - G matrix (subspace Hamiltonian) became ill-conditioned or contained NaN/Inf values
   - LAPACK's DSYEV eigensolver failed on the corrupted matrix

### Contributing Factors

1. **Delayed Error Detection**: The error manifested far from the root cause (in the diagonalizer rather than at input validation)

2. **Uninformative Error Message**: "DSYEV diagonalizer failed" gave no indication that the problem was with the user's active space specification

3. **Silent Integer Underflow**: Negative values from subtraction were not checked, allowing corruption to propagate

4. **Complex Call Chain**: Input validation → orbital setup → CI matrix construction → diagonalization spans multiple modules, making debugging difficult

## Solution

### Changes Made

Modified [psi4/src/psi4/libqt/ras_set.cc](psi4/src/psi4/libqt/ras_set.cc):

#### 1. Added `#include <string>` (Line 44)
Required for `std::to_string()` in error messages.

#### 2. Added Validation for RAS1/RAS2/RAS3/RAS4 (Lines 246-310)
After each `options.fill_int_array()` call:
```cpp
// Validate that RAS1 values don't exceed available orbitals in each irrep
for (irrep = 0; irrep < nirreps; irrep++) {
    if (ras_opi[0][irrep] > orbspi[irrep]) {
        throw InputException("ras_set3(): RAS1[" + std::to_string(irrep) + "] = " +
                            std::to_string(ras_opi[0][irrep]) + " exceeds available orbitals (" +
                            std::to_string(orbspi[irrep]) + ") in irrep " + std::to_string(irrep),
                            "RAS1", __FILE__, __LINE__);
    }
}
```
(Similar validation added for RAS2, RAS3, RAS4)

#### 3. Added Validation for ACTIVE Keyword (Lines 344-354)
```cpp
// Validate that ACTIVE values don't exceed available orbitals in each irrep
for (irrep = 0; irrep < nirreps; irrep++) {
    if (ras_opi[1][irrep] > orbspi[irrep]) {
        throw InputException("ras_set3(): ACTIVE[" + std::to_string(irrep) + "] = " +
                            std::to_string(ras_opi[1][irrep]) + " exceeds available orbitals (" +
                            std::to_string(orbspi[irrep]) + ") in irrep " + std::to_string(irrep),
                            "ACTIVE", __FILE__, __LINE__);
    }
}
```

### New Behavior

#### Error Message (After Fix)
```
InputException: ras_set3(): ACTIVE[5] = 1 exceeds available orbitals (0) in irrep 5
```

This error:
- ✅ Appears immediately at input validation time
- ✅ Identifies the exact problematic parameter (ACTIVE[5])
- ✅ Shows both the requested value (1) and available value (0)
- ✅ Clearly indicates which irrep is problematic (irrep 5)

## Impact

### Affected Calculations
- CASSCF (`energy('casscf')`)
- RASSCF (`energy('rasscf')`)
- Any DETCI-based multireference calculations with explicit active space specifications

### User Perspective

**Before**: Confusing error deep in the calculation, difficult to debug

**After**: Clear, actionable error message at the start of the calculation

### Developer Perspective

**Prevents**:
- Integer underflow/overflow bugs
- Corrupted internal data structures
- Wasted computational resources on doomed calculations

## Testing

### Test Case 1: Invalid ACTIVE Array
```python
molecule {
H
H 1 0.60
}

set basis cc-pVDZ
set docc [ 1, 0, 0, 0, 0, 0, 0, 0 ]
set active [1, 0, 0, 0, 0, 1, 0, 0 ]  # Should raise InputException
energy('casscf')
```

**Expected**: `InputException` with message about ACTIVE[5] exceeding available orbitals

### Test Case 2: Valid ACTIVE Array
```python
molecule {
H
H 1 0.60
}

set basis cc-pVDZ
set docc [ 1, 0, 0, 0, 0, 0, 0, 0 ]
set active [1, 0, 0, 0, 0, 0, 0, 0 ]  # Valid
energy('casscf')
```

**Expected**: Calculation proceeds normally

### Test Case 3: Invalid RAS2 Specification
```python
molecule {
N
N 1 1.1
}

set basis cc-pVDZ
set reference rohf
set ras2 [10, 0, 0, 0, 0, 10, 0, 0]  # Likely exceeds available orbitals
energy('rasscf')
```

**Expected**: `InputException` if RAS2 values exceed available orbitals in any irrep

## Related Code Locations

### Modified File
- [psi4/src/psi4/libqt/ras_set.cc](psi4/src/psi4/libqt/ras_set.cc): Input validation and orbital space setup

### Related (Unmodified) Files
- [psi4/src/psi4/detci/get_mo_info.cc:114-118](psi4/src/psi4/detci/get_mo_info.cc): Calls `ras_set3()`
- [psi4/src/psi4/detci/params.cc:930-968](psi4/src/psi4/detci/params.cc): Additional RAS parameter validation
- [psi4/src/psi4/detci/sem.cc:560](psi4/src/psi4/detci/sem.cc): DSYEV diagonalization (where error manifested)
- [psi4/driver/procrouting/proc.py:5512-5580](psi4/driver/procrouting/proc.py): Python driver for DETCAS methods

## Design Improvements

### Fail-Fast Principle
The fix implements **fail-fast** error handling: catch invalid input as early as possible, before it can corrupt internal state.

### User-Friendly Error Messages
Error messages now follow best practices:
- Identify the problematic parameter by name
- Show both expected and actual values
- Indicate which irrep is problematic
- Appear at input validation time, not deep in computation

### Defensive Programming
The validation serves as a **precondition check** that prevents:
- Invalid assumptions from propagating through the codebase
- Arithmetic underflow/overflow
- Undefined behavior from corrupted data structures

## Future Enhancements

### Additional Validation Opportunities

1. **Total Active Space Check**: Validate that total active electrons can fit in total active orbitals
   ```cpp
   int total_active_orbs = 0;
   int total_active_elec = 0;
   for (irrep = 0; irrep < nirreps; irrep++) {
       total_active_orbs += ras_opi[1][irrep];
       total_active_elec += /* count from docc/socc */;
   }
   if (total_active_elec > 2 * total_active_orbs) {
       // Error: too many electrons for active space
   }
   ```

2. **Cross-Validation**: Check that `DOCC + ACTIVE + VIRTUAL ≤ total_orbitals` per irrep

3. **Informative Suggestions**: When validation fails, suggest valid alternatives
   ```
   InputException: ACTIVE[5] = 1 exceeds available orbitals (0) in irrep 5
   Suggestion: Irreps with available orbitals: 0 (5 orbs), 1 (2 orbs), ...
   ```

4. **Unit Tests**: Add regression tests to ensure this specific failure mode is caught

## References

### Relevant Psi4 Documentation
- [DETCI Module Documentation](https://psicode.org/psi4manual/master/detci.html)
- [MCSCF Documentation](https://psicode.org/psi4manual/master/mcscf.html)
- [Symmetry in Psi4](https://psicode.org/psi4manual/master/psithonmol.html#symmetry)

### Point Group Theory
The issue often arises from misunderstanding which irreps contain molecular orbitals for a given molecule/basis combination. Users should:
1. Run an initial SCF calculation with `print 3` to see orbital counts per irrep
2. Ensure ACTIVE/RAS arrays only specify orbitals in irreps that exist
3. Use `symmetry c1` to disable symmetry if unsure

## Authors & Timeline

- **Issue Discovered**: User-reported CASSCF failure with "DSYEV diagonalizer failed" error
- **Root Cause Identified**: Missing validation in `ras_set3()` function
- **Fix Implemented**: 2025-10-23
- **Files Changed**: 1 file ([psi4/src/psi4/libqt/ras_set.cc](psi4/src/psi4/libqt/ras_set.cc))
- **Lines Added**: ~55 lines (validation code + comments)
- **Lines Deleted**: 0
- **Status**: Ready for testing and integration into main branch
