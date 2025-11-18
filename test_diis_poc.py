#!/usr/bin/env python3
"""
Basic functional test for DIIS POC implementation

This script runs a simple CCSD calculation to verify that the
libdiis-based DIIS implementation works correctly.

Usage:
    python test_diis_poc.py [--verbose]

Expected behavior:
    - Should complete without errors
    - Should show "POC: Using libdiis for DIIS extrapolation" message
    - Should show "DIIS: extrapolated with N vectors" messages
    - Final energy should match reference value within tolerance
"""

import psi4
import sys

def test_basic_ccsd():
    """Run basic CCSD test with DIIS POC"""

    # Set up molecule (H2O from cc1 test)
    psi4.core.clean()
    psi4.set_output_file("test_diis_poc_output.dat", False)

    h2o = psi4.geometry("""
        O
        H 1 0.97
        H 1 0.97 2 103.0
    """)

    psi4.set_options({
        'basis': '6-31G**',
        'scf_type': 'pk',
        'reference': 'rhf',
        'e_convergence': 10,
        'd_convergence': 8,
        'r_convergence': 7,
        'diis': True,
    })

    print("="*70)
    print("DIIS POC Functional Test")
    print("="*70)
    print("Testing: RHF-CCSD/6-31G** on H2O")
    print("Expected: Calculation should complete successfully")
    print("="*70)

    try:
        # Run CCSD energy calculation
        energy = psi4.energy('ccsd')

        # Reference values from cc1 test
        ref_scf = -76.0229427274435
        ref_ccsd_corr = -0.20823570806196
        ref_total = -76.2311784355056

        # Get computed values
        scf_energy = psi4.variable("SCF TOTAL ENERGY")
        ccsd_corr = psi4.variable("CCSD CORRELATION ENERGY")

        # Check results
        print("\n" + "="*70)
        print("RESULTS")
        print("="*70)
        print(f"SCF Energy:        {scf_energy:16.10f}")
        print(f"Reference:         {ref_scf:16.10f}")
        print(f"Difference:        {abs(scf_energy - ref_scf):16.10e}")
        print()
        print(f"CCSD Correlation:  {ccsd_corr:16.10f}")
        print(f"Reference:         {ref_ccsd_corr:16.10f}")
        print(f"Difference:        {abs(ccsd_corr - ref_ccsd_corr):16.10e}")
        print()
        print(f"Total Energy:      {energy:16.10f}")
        print(f"Reference:         {ref_total:16.10f}")
        print(f"Difference:        {abs(energy - ref_total):16.10e}")
        print("="*70)

        # Validate
        tol = 1.0e-9
        if abs(energy - ref_total) < tol:
            print("\n✓ TEST PASSED: Energy matches reference within tolerance")
            return True
        else:
            print(f"\n✗ TEST FAILED: Energy difference {abs(energy - ref_total):.2e} exceeds tolerance {tol:.2e}")
            return False

    except Exception as e:
        print(f"\n✗ TEST FAILED: Exception occurred: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    verbose = "--verbose" in sys.argv

    success = test_basic_ccsd()

    if success:
        print("\n" + "="*70)
        print("DIIS POC FUNCTIONAL TEST: SUCCESS")
        print("="*70)
        sys.exit(0)
    else:
        print("\n" + "="*70)
        print("DIIS POC FUNCTIONAL TEST: FAILURE")
        print("="*70)
        sys.exit(1)
