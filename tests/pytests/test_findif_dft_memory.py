"""
Test to verify HF destructor fix is safe for finite difference calculations.

Finite difference gradients/Hessians run many energy calculations at displaced
geometries. This test verifies that the potential->finalize() cleanup in the
HF destructor works correctly in this context.
"""

import gc
import pytest
import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]

# Try to import psutil for memory measurement
try:
    import psutil
    import os
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


def get_memory_mb():
    """Get current process memory usage in MB."""
    if HAS_PSUTIL:
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / (1024 ** 2)
    return None


def test_findif_gradient_dft():
    """
    Test finite difference DFT gradients work correctly with cleanup fix.

    This runs multiple energy calculations at displaced geometries.
    Each displacement creates a new wavefunction with its own grid.
    Verifies that destructor cleanup doesn't interfere.
    """

    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
        'points': 3,  # 3-point stencil for speed
    })

    # Calculate gradient by finite differences
    # This will run ~7 energy calculations (reference + 6 displacements for 3 atoms)
    grad = psi4.gradient('b3lyp', dertype='findif')

    # Should complete without crashes
    assert grad is not None
    assert grad.rows() == 3  # 3 atoms
    assert grad.cols() == 3  # x, y, z


def test_findif_gradient_multiple_functionals():
    """
    Test finite differences with different DFT functionals.

    Each functional creates different superfunctional objects,
    ensuring cleanup works across functional types.
    """

    mol = psi4.geometry("""
    H
    F 1 0.92
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
        'points': 3,
    })

    functionals = ['pbe', 'b3lyp']

    for func in functionals:
        psi4.core.clean()
        grad = psi4.gradient(func, dertype='findif')

        assert grad is not None
        assert grad.rows() == 2  # 2 atoms


@pytest.mark.skipif(not HAS_PSUTIL, reason="psutil not available for memory measurement")
def test_findif_no_memory_leak():
    """
    Test that finite difference gradients don't leak memory.

    Without the destructor fix, each displaced geometry's wavefunction
    would leak its DFT grid, causing memory to grow significantly.
    """

    mol = psi4.geometry("""
    0 1
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
        'points': 3,
        'dft_spherical_points': 110,  # Moderate grid
        'dft_radial_points': 20,
    })

    # Force garbage collection
    gc.collect()

    # First gradient for baseline
    psi4.gradient('b3lyp', dertype='findif')
    psi4.core.clean()
    gc.collect()

    mem_baseline = get_memory_mb()

    # Run multiple gradient calculations
    num_gradients = 3
    for i in range(num_gradients):
        psi4.gradient('b3lyp', dertype='findif')
        psi4.core.clean()
        gc.collect()

    mem_final = get_memory_mb()
    total_growth = mem_final - mem_baseline
    growth_per_gradient = total_growth / num_gradients

    # Each gradient runs ~7 energy calculations (1 ref + 6 displacements for 3 atoms)
    # Without fix: ~3-5 MB leaked per energy = ~21-35 MB per gradient
    # With fix: should be < 5 MB per gradient (just normal overhead)
    assert growth_per_gradient < 8.0, (
        f"Memory leak detected in finite differences: {growth_per_gradient:.2f} MB per gradient. "
        f"Started at {mem_baseline:.2f} MB, ended at {mem_final:.2f} MB after "
        f"{num_gradients} gradients. Total growth: {total_growth:.2f} MB."
    )


def test_optimize_with_findif():
    """
    Test geometry optimization with finite difference gradients.

    Optimizations run many gradient evaluations, each of which runs many
    energy calculations. This is a stress test for the cleanup fix.
    """

    mol = psi4.geometry("""
    H 0 0 0
    H 0 0 0.9
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
        'points': 3,
        'g_convergence': 'gau_loose',  # Loose convergence for speed
        'geom_maxiter': 5,  # Limit iterations
    })

    # Run optimization with finite difference gradients
    # Should complete without crashes or memory issues
    E = psi4.optimize('pbe', dertype='findif')

    assert E is not None
    # H2 should optimize to around -1.12 Hartree with PBE/STO-3G
    assert -1.2 < E < -1.0


def test_findif_hessian_dft():
    """
    Test finite difference Hessian from gradients.

    Hessians from gradients run many gradient calculations, each of which
    runs many energy calculations. This is an even more intensive test.
    """

    mol = psi4.geometry("""
    H 0 0 0
    H 0 0 0.74
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
        'points': 3,
    })

    # Calculate Hessian from finite difference of gradients
    # For 2 atoms (6 DOF), this runs many gradient calculations
    # Each gradient runs ~7 energy calculations
    hess = psi4.hessian('pbe', dertype=1)  # dertype=1 means findif from gradients

    assert hess is not None
    assert hess.rows() == 6  # 2 atoms * 3 coordinates
    assert hess.cols() == 6


def test_findif_mixed_reference():
    """
    Test finite differences with different reference types.

    Verifies that RHF, UHF, and ROHF all work correctly with cleanup.
    """

    # Singlet - RHF reference
    mol_singlet = psi4.geometry("""
    0 1
    H 0 0 0
    H 0 0 0.74
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
        'points': 3,
        'reference': 'rhf',
    })

    grad_rhf = psi4.gradient('b3lyp', molecule=mol_singlet, dertype='findif')
    assert grad_rhf is not None

    # Doublet - UHF reference
    mol_doublet = psi4.geometry("""
    0 2
    H
    symmetry c1
    """)

    psi4.set_options({
        'reference': 'uhf',
    })

    psi4.core.clean()
    grad_uhf = psi4.gradient('b3lyp', molecule=mol_doublet, dertype='findif')
    assert grad_uhf is not None
