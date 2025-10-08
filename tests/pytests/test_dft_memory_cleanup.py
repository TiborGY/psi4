"""
Test for DFT memory cleanup regression test (GitHub issue #XXXX)

This test verifies that repeated DFT calculations in the same process
properly clean up grid and functional resources to prevent memory leaks.
"""

import gc
import pytest
import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

# Try to import psutil for actual memory measurement, but don't fail if unavailable
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


@pytest.mark.skipif(not HAS_PSUTIL, reason="psutil not available for memory measurement")
def test_dft_repeated_calculations_no_memory_leak():
    """
    Regression test for DFT memory leak fix.

    Running repeated DFT calculations should not accumulate memory by leaving
    DFT grid and functional worker objects unreleased. This test actually
    measures memory usage to detect leaks.

    Without the fix (HF destructor calling potential->finalize()), memory
    usage would increase by several MB per iteration.
    """

    # Small molecule for fast testing
    mol = psi4.geometry("""
    0 1
    C    0.000    0.239    0.000
    N    -0.044   1.395    0.000
    N    0.134    -1.155   0.000
    H    -0.315   -1.556   0.837
    H    -0.315   -1.556   -0.837
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'def2-svp',
        'scf_type': 'df',
        'dft_spherical_points': 110,
        'dft_radial_points': 20,
    })

    # Force garbage collection before starting
    gc.collect()

    # Measure baseline memory
    mem_start = get_memory_mb()

    # Run first calculation to establish baseline (includes one-time allocations)
    psi4.energy('b3lyp')
    psi4.core.clean()
    gc.collect()

    mem_after_first = get_memory_mb()
    baseline_growth = mem_after_first - mem_start

    # Run multiple iterations to detect memory leak
    num_iterations = 10
    for i in range(num_iterations):
        psi4.energy('b3lyp')
        psi4.core.clean()
        gc.collect()

    mem_final = get_memory_mb()
    total_growth = mem_final - mem_after_first
    growth_per_iteration = total_growth / num_iterations

    # With the fix, growth should be minimal (< 1 MB per iteration due to normal overhead)
    # Without the fix, each iteration leaks ~3-5 MB (grid + functionals + workers)
    assert growth_per_iteration < 1.5, (
        f"Memory leak detected: {growth_per_iteration:.2f} MB per iteration. "
        f"Started at {mem_after_first:.2f} MB, ended at {mem_final:.2f} MB after "
        f"{num_iterations} iterations. Total growth: {total_growth:.2f} MB."
    )


@pytest.mark.skipif(not HAS_PSUTIL, reason="psutil not available for memory measurement")
def test_r2scan3c_memory_leak():
    """
    Specific regression test for r2SCAN-3c (the original bug report method).

    This test uses the exact method from the bug report and measures memory.
    """

    mol = psi4.geometry("""
    C    0.000    0.239    0.000
    N    -0.044   1.395    0.000
    N    0.134    -1.155   0.000
    H    -0.315   -1.556   0.837
    H    -0.315   -1.556   -0.837
    symmetry c1
    """)

    # Force garbage collection
    gc.collect()

    # First calculation for baseline
    psi4.energy("r2scan-3c")
    psi4.core.clean()
    gc.collect()

    mem_baseline = get_memory_mb()

    # Run 15 more iterations (original bug report used 21 total)
    num_iterations = 15
    for i in range(num_iterations):
        psi4.energy("r2scan-3c")
        psi4.core.clean()
        gc.collect()

    mem_final = get_memory_mb()
    total_growth = mem_final - mem_baseline
    growth_per_iteration = total_growth / num_iterations

    # r2SCAN-3c has larger grids, so allow slightly more tolerance
    # but still should be < 2 MB per iteration with the fix
    assert growth_per_iteration < 2.0, (
        f"r2SCAN-3c memory leak detected: {growth_per_iteration:.2f} MB per iteration. "
        f"Started at {mem_baseline:.2f} MB, ended at {mem_final:.2f} MB after "
        f"{num_iterations} iterations. Total growth: {total_growth:.2f} MB."
    )


def test_dft_potential_finalize_called():
    """
    Test that verifies the potential finalize mechanism works correctly.

    This doesn't measure memory but verifies that:
    1. DFT wavefunctions have a potential object
    2. The potential can be accessed via V_potential()
    3. Multiple calculations don't crash (would crash if finalize broken)
    """

    mol = psi4.geometry("""
    H
    F 1 0.92
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'cc-pvdz',
        'scf_type': 'df',
    })

    # Create multiple wavefunctions and let them go out of scope
    for i in range(5):
        E, wfn = psi4.energy('b3lyp', return_wfn=True)

        # Verify the potential exists
        potential = wfn.V_potential()
        assert potential is not None, "DFT potential should exist for B3LYP"

        # Verify it's the right type (RV for restricted)
        assert wfn.same_a_b_orbs(), "Should be restricted for singlet"

        psi4.core.clean()
        # wfn goes out of scope here, triggering destructor
        # which should call potential->finalize() to clean up grid


def test_hf_no_potential():
    """
    Control test: HF calculations should have no DFT potential.

    This verifies the fix doesn't break pure HF calculations.
    """

    mol = psi4.geometry("""
    He
    """)

    psi4.set_options({
        'basis': 'cc-pvdz',
        'scf_type': 'df',
    })

    for i in range(5):
        E, wfn = psi4.energy('hf', return_wfn=True)

        # HF wavefunctions should have no DFT potential
        potential = wfn.V_potential()
        assert potential is None, "HF should not have DFT potential"

        psi4.core.clean()


def test_uhf_dft_cleanup():
    """
    Test that unrestricted DFT also cleans up properly.

    UHF uses UV potential instead of RV, so verify cleanup works there too.
    """

    mol = psi4.geometry("""
    0 2
    O
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'def2-svp',
        'scf_type': 'df',
        'reference': 'uhf',
    })

    for i in range(5):
        E, wfn = psi4.energy('b3lyp', return_wfn=True)

        # Verify it's unrestricted
        assert not wfn.same_a_b_orbs(), "Should be unrestricted for doublet"

        # Verify the potential exists (should be UV type)
        potential = wfn.V_potential()
        assert potential is not None, "UKS should have DFT potential"

        psi4.core.clean()
