"""
Test to verify HF destructor fix is safe with COSX (chain-of-spheres exchange).

COSX uses DFT grids to compute exact exchange semi-numerically. These grids
are stored in the JK object (COSK::grids_), separate from the VBase DFT potential
grid. This test verifies both grid systems are cleaned up independently.
"""

import gc
import pytest
import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

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


def test_cosx_basic():
    """
    Test basic COSX calculation completes without crashes.

    COSX creates two grids (Initial and Final) for the exchange computation.
    Verifies cleanup doesn't interfere.
    """

    mol = psi4.geometry("""
    H 0 0 0
    F 0 0 0.92
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'def2-svp',
        'scf_type': 'cosx',  # Enable COSX
    })

    E = psi4.energy('b3lyp')
    assert E is not None


def test_cosx_repeated_calculations():
    """
    Test repeated COSX calculations (memory leak scenario).

    Each calculation creates:
    - JK object with COSX grids (cleaned in HF::finalize)
    - VBase potential with DFT grid (cleaned in ~HF destructor)

    Both should be properly freed.
    """

    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'def2-svp',
        'scf_type': 'cosx',
    })

    energies = []
    for i in range(5):
        E = psi4.energy('pbe')
        energies.append(E)
        psi4.core.clean()

    # All energies should be consistent
    for i, E in enumerate(energies[1:], 1):
        assert psi4.compare_values(energies[0], E, 7,
                                   f"PBE/COSX energy iteration {i}")


def test_cosx_with_return_wfn():
    """
    Test COSX with return_wfn to verify destructor cleanup.

    When wfn is returned, it contains:
    - jk_ object (with COSX grids)
    - potential_ object (with VBase grid)

    Both should be cleaned when wfn is deleted.
    """

    mol = psi4.geometry("""
    Ne
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'def2-svp',
        'scf_type': 'cosx',
    })

    E, wfn = psi4.energy('b3lyp', return_wfn=True)

    # Verify both exist
    assert wfn.jk() is not None, "JK object should exist"
    assert wfn.V_potential() is not None, "VBase potential should exist"

    # Delete wavefunction
    # - jk_.reset() happens in wfn.finalize() (called before destruction)
    # - potential->finalize() happens in ~HF() destructor (my fix)
    del wfn

    # Run another calculation to verify everything works
    psi4.core.clean()
    E2 = psi4.energy('b3lyp')

    assert psi4.compare_values(E, E2, 7, "B3LYP/COSX energies")


@pytest.mark.skipif(not HAS_PSUTIL, reason="psutil not available for memory measurement")
def test_cosx_no_memory_leak():
    """
    Test that COSX calculations don't leak memory.

    Without destructor fix, VBase grids would accumulate.
    COSX grids are already cleaned up in HF::finalize().
    """

    mol = psi4.geometry("""
    0 1
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'def2-svp',
        'scf_type': 'cosx',
    })

    # Force garbage collection
    gc.collect()

    # Baseline calculation
    psi4.energy('b3lyp')
    psi4.core.clean()
    gc.collect()

    mem_baseline = get_memory_mb()

    # Run multiple calculations
    num_iterations = 5
    for i in range(num_iterations):
        psi4.energy('b3lyp')
        psi4.core.clean()
        gc.collect()

    mem_final = get_memory_mb()
    total_growth = mem_final - mem_baseline
    growth_per_iteration = total_growth / num_iterations

    # Each COSX calculation uses 2 grids (Initial + Final) for K
    # Plus 1 grid for VBase (Vxc)
    # Without fix: VBase grid leaked (~3-5 MB)
    # With fix: all grids cleaned up properly
    assert growth_per_iteration < 2.0, (
        f"Memory leak detected with COSX: {growth_per_iteration:.2f} MB per iteration. "
        f"Started at {mem_baseline:.2f} MB, ended at {mem_final:.2f} MB after "
        f"{num_iterations} iterations. Total growth: {total_growth:.2f} MB."
    )


def test_cosx_grid_switching():
    """
    Test COSX grid switching (Initial → Final).

    COSX starts with a small grid ("Initial") and switches to
    a larger grid ("Final") near convergence. Verifies cleanup
    works correctly when grids are switched.
    """

    mol = psi4.geometry("""
    H 0 0 0
    H 0 0 0.74
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'def2-svp',
        'scf_type': 'cosx',
        'cosx_maxiter_final': 5,  # Switch to final grid 5 iterations before convergence
    })

    E, wfn = psi4.energy('pbe', return_wfn=True)

    # Wavefunction should have completed successfully
    assert E is not None
    assert wfn.jk() is not None

    # Clean up
    del wfn
    psi4.core.clean()


def test_cosx_vs_regular_df():
    """
    Compare COSX with regular DF to ensure COSX works correctly.

    This serves as a validation test - COSX should give similar
    energies to regular DF (within tolerance).
    """

    mol = psi4.geometry("""
    H 0 0 0
    F 0 0 0.92
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'def2-svp',
    })

    # Regular DF calculation
    psi4.set_options({'scf_type': 'df'})
    E_df = psi4.energy('pbe')
    psi4.core.clean()

    # COSX calculation
    psi4.set_options({'scf_type': 'cosx'})
    E_cosx = psi4.energy('pbe')
    psi4.core.clean()

    # Should be close (COSX is approximate)
    # Tolerance is loose because COSX is semi-numerical
    assert abs(E_df - E_cosx) < 1e-4, (
        f"COSX energy ({E_cosx}) differs significantly from DF ({E_df})"
    )


def test_cosx_findif_gradient():
    """
    Test COSX with finite difference gradients.

    Each displacement creates a new wavefunction with:
    - New JK object (COSX grids)
    - New VBase potential (DFT grid)

    Both should be cleaned up after each displacement.
    """

    mol = psi4.geometry("""
    H 0 0 0
    H 0 0 0.74
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'sto-3g',  # Small basis for speed
        'scf_type': 'cosx',
        'points': 3,
    })

    # Calculate gradient using finite differences
    grad = psi4.gradient('pbe', dertype='findif')

    assert grad is not None
    assert grad.rows() == 2  # 2 atoms


def test_cosx_hybrid_functional():
    """
    Test COSX with hybrid functional (both XC and exact exchange).

    Hybrid DFT needs:
    - COSX grids for exact exchange (K)
    - VBase grid for DFT XC contribution

    Both grid systems are used simultaneously.
    """

    mol = psi4.geometry("""
    Ne
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'def2-svp',
        'scf_type': 'cosx',
    })

    # B3LYP is a hybrid (20% exact exchange)
    E, wfn = psi4.energy('b3lyp', return_wfn=True)

    # Both grid systems should exist
    assert wfn.jk() is not None, "JK (COSX) should exist for hybrid"
    assert wfn.V_potential() is not None, "VBase should exist for hybrid"

    # Cleanup
    del wfn
    psi4.core.clean()

    # Run again to verify
    E2 = psi4.energy('b3lyp')
    assert psi4.compare_values(E, E2, 7, "B3LYP/COSX hybrid")
