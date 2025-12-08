# CASSCF/RASSCF Active Space Input Validation

## Summary

This change adds early input validation to `ras_set3()` to catch invalid active space specifications before they cause confusing downstream errors.

## What This Fix Does

Adds per-irrep validation for RAS1, RAS2, RAS3, RAS4, and ACTIVE keywords to ensure that requested orbital counts do not exceed available orbitals in each irreducible representation.

### Before This Fix

Invalid active space specifications (e.g., requesting orbitals in an irrep with no available MOs) could lead to:
- Negative orbital counts in internal calculations
- Corrupted orbital reordering arrays
- Cryptic errors deep in the CI code (e.g., "Invalid point value" or potentially diagonalization failures)

### After This Fix

Invalid specifications are caught immediately with a clear error message:
```
InputException: ras_set3(): ACTIVE[4] = 1 exceeds available orbitals (0) in irrep 4
```

## Example: Invalid Input That Is Now Caught

For H₂ in D2h symmetry with cc-pVDZ basis, the orbital distribution is:
- Ag (0): 3 orbitals
- B1g (1): 0 orbitals
- B2g (2): 1 orbital
- B3g (3): 1 orbital
- Au (4): 0 orbitals
- B1u (5): 3 orbitals
- B2u (6): 1 orbital
- B3u (7): 1 orbital

The following input is **invalid** because Au (index 4) has no orbitals:
```python
set active [1, 0, 0, 0, 1, 0, 0, 0]  # Requests 1 orbital in Au - invalid!
```

This fix catches this error immediately rather than allowing it to propagate.

## Changes Made

### Modified File
- [psi4/src/psi4/libqt/ras_set.cc](psi4/src/psi4/libqt/ras_set.cc)

### Validation Added
After each `options.fill_int_array()` call for RAS1, RAS2, RAS3, RAS4, and ACTIVE:
```cpp
for (irrep = 0; irrep < nirreps; irrep++) {
    if (ras_opi[X][irrep] > orbspi[irrep]) {
        throw InputException("ras_set3(): KEYWORD[" + std::to_string(irrep) + "] = " +
                            std::to_string(ras_opi[X][irrep]) + " exceeds available orbitals (" +
                            std::to_string(orbspi[irrep]) + ") in irrep " + std::to_string(irrep),
                            "KEYWORD", __FILE__, __LINE__);
    }
}
```

## Affected Calculations

- CASSCF (`energy('casscf')`)
- RASSCF (`energy('rasscf')`)
- Any DETCI-based multireference calculations with explicit RAS/ACTIVE specifications

## Relationship to Issue #3096

**Important clarification:** This fix does **not** resolve GitHub issue #3096.

Issue #3096 reported a CASSCF failure with a **valid** active space specification:
```python
set active [1, 0, 0, 0, 0, 1, 0, 0]  # Valid (2e,2o) in Ag and B1u
```

This input requests orbitals in Ag (index 0) and B1u (index 5), both of which have available orbitals for H₂/cc-pVDZ. The validation added here would not catch this case because the input is valid.

The actual cause of the DSYEV failure in issue #3096 remains undiagnosed. As noted by TiborGY in the issue: "No idea what is causing DSYEV to fail."

## Design Rationale

### Fail-Fast Principle
Invalid input is caught at validation time rather than causing mysterious failures deep in the computation.

### User-Friendly Error Messages
Error messages identify:
- The problematic keyword (ACTIVE, RAS1, etc.)
- The specific irrep that's problematic
- Both the requested and available orbital counts

### Defensive Programming
The validation serves as a precondition check that prevents:
- Integer underflow in subsequent calculations
- Corrupted orbital reordering arrays
- Wasted computational resources on doomed calculations

## Testing

Test cases added to `tests/pytests/test_raises.py`:

1. `test_casscf_invalid_active_irrep()` - Verifies that invalid ACTIVE specifications raise clear errors
2. `test_casscf_valid_active_succeeds()` - Verifies that valid CASSCF calculations still work

## Future Work

The underlying cause of the DSYEV failure in issue #3096 (with valid input) remains to be investigated. Possible areas to explore:
- Numerical issues in CI matrix construction for small active spaces
- Edge cases in the Davidson/SEM algorithm
- Symmetry-specific orbital handling
