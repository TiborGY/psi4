#!/usr/bin/env python3
"""
Compare original DIIS vs libdiis POC implementation

This script runs the same CCSD calculation twice:
1. With USE_LIBDIIS_POC disabled (original implementation)
2. With USE_LIBDIIS_POC enabled (POC implementation)

Then compares:
- Final energies
- Iteration counts
- Convergence behavior
- Performance (optional)

Usage:
    python compare_diis_implementations.py [--molecule MOLECULE] [--basis BASIS]

Note: This script assumes you have TWO builds of Psi4:
    - Default build (original DIIS)
    - Build with -DCMAKE_CXX_FLAGS="-DUSE_LIBDIIS_POC" (POC DIIS)

For now, this is a template - actual comparison requires build infrastructure.
"""

import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Compare DIIS implementations")
    parser.add_argument("--molecule", default="h2o",
                        help="Molecule to test (h2o, nh3, ch4)")
    parser.add_argument("--basis", default="6-31G**",
                        help="Basis set to use")
    parser.add_argument("--method", default="ccsd",
                        help="CC method to test")
    parser.add_argument("--verbose", action="store_true",
                        help="Print detailed output")
    return parser.parse_args()

def get_molecule(name):
    """Return molecule geometry"""
    molecules = {
        'h2o': """
            O
            H 1 0.97
            H 1 0.97 2 103.0
        """,
        'nh3': """
            N
            H 1 1.008
            H 1 1.008 2 106.7
            H 1 1.008 2 106.7 3 120.0
        """,
        'ch4': """
            C
            H 1 1.09
            H 1 1.09 2 109.5
            H 1 1.09 2 109.5 3 120.0
            H 1 1.09 2 109.5 3 -120.0
        """
    }
    return molecules.get(name.lower(), molecules['h2o'])

def run_comparison(molecule_name, basis, method, verbose=False):
    """
    Run comparison between implementations

    NOTE: This is a template for the comparison workflow.
    Actual implementation requires:
    1. Two separate Psi4 builds (with/without POC flag)
    2. Subprocess management to run each build
    3. Output parsing to extract energies and iteration counts
    """

    print("="*70)
    print("DIIS Implementation Comparison")
    print("="*70)
    print(f"Molecule: {molecule_name}")
    print(f"Basis:    {basis}")
    print(f"Method:   {method}")
    print("="*70)

    print("\nNOTE: Full comparison requires two separate Psi4 builds:")
    print("  1. Default build (original DIIS)")
    print("  2. Build with -DCMAKE_CXX_FLAGS=\"-DUSE_LIBDIIS_POC\"")
    print("\nCurrently showing analysis framework only.")
    print("="*70)

    # Template for what the comparison would check
    comparison_metrics = {
        'Energy Accuracy': {
            'description': 'Final energy difference',
            'tolerance': 1.0e-9,
            'status': 'PENDING'
        },
        'Iteration Count': {
            'description': 'Number of CCSD iterations',
            'tolerance': 0,  # Should be identical
            'status': 'PENDING'
        },
        'Convergence Rate': {
            'description': 'Energy convergence per iteration',
            'tolerance': 1.0e-10,
            'status': 'PENDING'
        },
        'DIIS Subspace Size': {
            'description': 'Number of DIIS vectors used',
            'tolerance': 0,  # Should be identical
            'status': 'PENDING'
        }
    }

    print("\nComparison Metrics:")
    print("-"*70)
    for metric, details in comparison_metrics.items():
        print(f"  {metric:20s}: {details['description']}")
        print(f"  {'':20s}  Tolerance: {details['tolerance']}")
        print(f"  {'':20s}  Status: {details['status']}")
        print()

    return True

def main():
    args = parse_args()

    success = run_comparison(
        args.molecule,
        args.basis,
        args.method,
        args.verbose
    )

    if success:
        print("\nComparison framework ready.")
        print("To run actual comparison:")
        print("  1. Build Psi4 without POC flag")
        print("  2. Run test and save output")
        print("  3. Build Psi4 with POC flag")
        print("  4. Run test and save output")
        print("  5. Compare outputs")
        return 0
    else:
        return 1

if __name__ == "__main__":
    sys.exit(main())
