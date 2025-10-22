# Fix for Invalid ACTIVE Space Specification Error

## Issue Description

CASSCF calculations with invalid ACTIVE orbital specifications were failing deep in the DETCI code with a cryptic error:

```
RuntimeError:
Fatal Error: DSYEV diagonalizer failed in DETCI SEM!
Error occurred in file: /home/work/psi4/psi4/src/psi4/detci/sem.cc on line: 561
```

This occurred when users specified ACTIVE orbital spaces that were incompatible with the actual molecular orbitals available from the basis set.

### Example Failing Input

```python
molecule {
H
H 1 0.60
}

set basis cc-pVDZ
set docc [ 1, 0, 0, 0, 0, 0, 0, 0 ]
# STO-3G active space (wrong for cc-pVDZ!)
set active [1, 0, 0, 0, 0, 1, 0, 0 ]
energy('casscf')
```

This input requests active orbitals in irreps 0 and 5, but the actual H2/cc-pVDZ orbitals may not be distributed that way.

## Root Cause

1. **Insufficient Validation**: The `ras_set3()` function in `psi4/src/psi4/libqt/ras_set.cc` had only a warning (not an error) when orbital specifications were invalid.

2. **Late Failure**: Invalid orbital configurations were allowed to proceed through setup and into the DETCI CI diagonalizer, where they caused numerical failures in LAPACK's DSYEV routine.

3. **Unclear Error Messages**: The DSYEV failure provided no context about why the diagonalization failed.

## Solution

### Changes Made

**File: `psi4/src/psi4/libqt/ras_set.cc`**

1. **Added comprehensive validation** (lines 303-348) to check:
   - Requested orbital spaces don't exceed available orbitals per irrep
   - Sufficient active orbitals exist for all occupied electrons
   - No negative orbital counts after allocation

2. **Replaced warnings with exceptions**: Changed from `outfile->Printf()` warnings to `throw InputException()` with detailed error messages

3. **Added informative error messages** that show:
   - Which irrep has the problem
   - How many orbitals are available
   - How many were requested in each category (FROZEN_DOCC, RESTRICTED_DOCC, ACTIVE, etc.)
   - What the user should adjust

4. **Added header**: Included `<string>` for `std::string` and `std::to_string`

### Example New Error Message

```
InputException: Orbital specification exceeds available orbitals for irrep 5.
  Available orbitals in this irrep: 0
  Frozen DOCC:      0
  Restricted DOCC:  0
  RAS 1:            0
  ACTIVE (RAS 2):   1
  Frozen UOCC:      0
  Total requested:  1
Please adjust your ACTIVE, DOCC, FROZEN_DOCC, RESTRICTED_DOCC, FROZEN_UOCC, or RESTRICTED_UOCC specification.
```

This clearly shows that irrep 5 has 0 orbitals available, but the ACTIVE specification requests 1.

## Building and Testing

### Build Instructions

```bash
# From the psi4 repository root
cd /c/Users/Szuperhab/Documents/psi4

# Configure (adjust paths as needed)
cmake -S. -Bobjdir -DCMAKE_INSTALL_PREFIX=/path/to/install

# Build with parallel compilation
cmake --build objdir -j$(getconf _NPROCESSORS_ONLN)

# Install
cmake --install objdir
```

### Testing the Fix

**Quick Manual Test:**

```bash
cd /c/Users/Szuperhab/Documents/psi4
python test_invalid_active.py
```

This should now print:
```
SUCCESS: Caught expected validation error
Error message:
[Informative error about orbital mismatch]
```

**Run Pytest Suite:**

```bash
cd objdir
pytest ../tests/pytests/test_invalid_active_space.py -v
```

All three tests should pass:
- `test_invalid_active_exceeds_available_orbitals` - Catches over-specified ACTIVE
- `test_invalid_active_insufficient_for_electrons` - Catches under-specified ACTIVE
- `test_valid_active_space` - Ensures valid inputs still work

**Run Full Test Suite:**

```bash
cd objdir
ctest -j$(getconf _NPROCESSORS_ONLN)
```

## Benefits

1. **Early Detection**: Invalid configurations caught during setup, not deep in computation
2. **Clear Diagnostics**: Users get actionable error messages explaining exactly what's wrong
3. **Time Savings**: No need to wait for SCF and integral transformation before discovering the error
4. **Better User Experience**: Guides users toward correct input specification

## Backward Compatibility

- Valid ACTIVE space specifications continue to work exactly as before
- Only previously invalid (and failing) inputs now raise exceptions earlier with better messages
- No changes to API or valid input syntax

## Related Files

- `psi4/src/psi4/libqt/ras_set.cc` - Main fix location
- `tests/pytests/test_invalid_active_space.py` - New test cases
- `test_invalid_active.py` - Quick manual test script

## Future Improvements

Potential enhancements for even better validation:

1. Add check in CIWavefunction constructor for zero-dimensional CI spaces
2. Pre-compute and display expected orbital distributions for common molecules
3. Add validation for RAS1/RAS3/RAS4 specifications similar to ACTIVE
4. Provide suggestions for correct ACTIVE space based on SCF orbitals
