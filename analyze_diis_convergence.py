#!/usr/bin/env python3
"""
Analyze DIIS convergence behavior from Psi4 output

This script parses Psi4 CCSD output files and extracts:
- Iteration-by-iteration energy values
- RMS and MAX amplitude changes
- DIIS subspace sizes
- Convergence patterns

Usage:
    python analyze_diis_convergence.py output.dat

Output:
    - Prints convergence table
    - Optionally plots convergence graph (if matplotlib available)
"""

import re
import sys
import argparse
from pathlib import Path

def parse_ccsd_output(filename):
    """
    Parse CCSD iteration data from output file

    Returns dict with:
        - iterations: list of iteration numbers
        - energies: list of energies
        - delta_E: list of energy changes
        - rms: list of RMS values
        - diis_info: list of DIIS messages
    """

    with open(filename, 'r') as f:
        content = f.read()

    # Find CCSD iteration block
    # Typical line:
    # @RHF CCSD Iteration   1: Energy = -76.228... dE = -1.000... t1 = 0.01 t2 = 0.05
    pattern = r'@\w+\s+CCSD\s+Iteration\s+(\d+):\s+Energy\s+=\s+([-\d.]+)'

    iterations = []
    energies = []

    for match in re.finditer(pattern, content):
        iteration = int(match.group(1))
        energy = float(match.group(2))
        iterations.append(iteration)
        energies.append(energy)

    # Find DIIS messages
    diis_pattern = r'DIIS:\s+(.+)'
    diis_messages = re.findall(diis_pattern, content)

    # Calculate energy changes
    delta_E = [0.0]  # First iteration has no previous
    for i in range(1, len(energies)):
        delta_E.append(energies[i] - energies[i-1])

    return {
        'iterations': iterations,
        'energies': energies,
        'delta_E': delta_E,
        'diis_messages': diis_messages
    }

def print_convergence_table(data):
    """Print nicely formatted convergence table"""

    print("\n" + "="*80)
    print("CCSD CONVERGENCE ANALYSIS")
    print("="*80)

    if not data['iterations']:
        print("No CCSD iteration data found in output file")
        return

    print(f"\n{'Iter':>5} {'Energy':>18} {'ΔE':>15} {'DIIS Info':>30}")
    print("-"*80)

    for i, (iter_num, energy, delta) in enumerate(zip(
        data['iterations'],
        data['energies'],
        data['delta_E']
    )):
        diis_info = data['diis_messages'][i] if i < len(data['diis_messages']) else ""
        diis_info = diis_info[:30]  # Truncate if too long

        print(f"{iter_num:>5d} {energy:>18.10f} {delta:>15.8e} {diis_info:>30}")

    print("-"*80)
    print(f"Total iterations: {len(data['iterations'])}")
    print(f"Final energy:     {data['energies'][-1]:.10f}")
    print(f"Final ΔE:         {data['delta_E'][-1]:.2e}")
    print("="*80)

def plot_convergence(data):
    """Plot convergence graph (if matplotlib available)"""

    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("\nNote: matplotlib not available, skipping plot")
        return

    if not data['iterations']:
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Plot 1: Energy vs iteration
    ax1.plot(data['iterations'], data['energies'], 'b-o', label='Energy')
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Energy (Hartree)')
    ax1.set_title('CCSD Energy Convergence')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Plot 2: Energy change (log scale)
    abs_delta = [abs(d) for d in data['delta_E'][1:]]  # Skip first (0.0)
    iters = data['iterations'][1:]

    ax2.semilogy(iters, abs_delta, 'r-o', label='|ΔE|')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('|ΔE| (Hartree, log scale)')
    ax2.set_title('CCSD Energy Change per Iteration')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.legend()

    plt.tight_layout()
    plt.savefig('diis_convergence.png', dpi=150)
    print("\nConvergence plot saved to: diis_convergence.png")

def main():
    parser = argparse.ArgumentParser(description='Analyze DIIS convergence')
    parser.add_argument('output_file', help='Psi4 output file to analyze')
    parser.add_argument('--plot', action='store_true', help='Generate plot')
    args = parser.parse_args()

    if not Path(args.output_file).exists():
        print(f"Error: File not found: {args.output_file}")
        return 1

    print(f"Analyzing: {args.output_file}")

    data = parse_ccsd_output(args.output_file)
    print_convergence_table(data)

    if args.plot:
        plot_convergence(data)

    return 0

if __name__ == "__main__":
    sys.exit(main())
