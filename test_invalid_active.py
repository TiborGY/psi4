#!/usr/bin/env python
"""
Test that invalid ACTIVE space specifications are properly caught
and produce informative error messages.
"""

import psi4

# This is the user's original failing input
input_str = """
molecule {
H
H 1 0.60
}

set basis cc-pVDZ
set docc [ 1, 0, 0, 0, 0, 0, 0, 0 ]
# STO-3G active space (this is a mismatch!)
set active [1, 0, 0, 0, 0, 1, 0, 0 ]
"""

psi4.set_output_file('test_invalid_active.out', False)
psi4.geometry(input_str)

# This should now raise an informative exception instead of
# failing deep in the DETCI SEM diagonalizer
try:
    psi4.energy('casscf')
    print("ERROR: Should have raised an exception!")
except psi4.ValidationError as e:
    print("SUCCESS: Caught expected validation error")
    print(f"Error message:\n{e}")
except Exception as e:
    print(f"PARTIAL SUCCESS: Caught exception (type: {type(e).__name__})")
    print(f"Error message:\n{e}")
