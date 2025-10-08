"""
Test to verify HF destructor fix doesn't cause issues.

This test verifies that calling potential->finalize() in the HF destructor
is safe even when HF::finalize() was already called during SCF.
"""

import pytest
import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def test_double_finalize_safe():
    """
    Test that calling finalize() twice doesn't cause crashes.

    During normal SCF flow:
    1. SCF completes and calls wfn.finalize() (which cleans up JK, etc.)
    2. Wavefunction destructor calls potential->finalize() (which resets grid_)

    This should be safe since grid_.reset() is idempotent.
    """

    mol = psi4.geometry("""
    H
    F 1 0.92
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
    })

    # Run a DFT calculation - this will call wfn.finalize() at the end
    E = psi4.energy('b3lyp')

    # The wavefunction goes out of scope here, destructor calls potential->finalize() again
    # This should not crash

    # Run another calculation to verify everything still works
    psi4.core.clean()
    E2 = psi4.energy('b3lyp')

    assert psi4.compare_values(E, E2, 9, "B3LYP energies should match")


def test_wfn_manual_deletion():
    """
    Test that explicitly deleting wavefunction works correctly.
    """

    mol = psi4.geometry("""
    Ne
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
    })

    E, wfn = psi4.energy('pbe', return_wfn=True)

    # Verify potential exists
    assert wfn.V_potential() is not None

    # Explicitly delete the wavefunction
    del wfn

    # Should not crash, run another calculation
    psi4.core.clean()
    E2 = psi4.energy('pbe')

    assert psi4.compare_values(E, E2, 9, "PBE energies should match")


def test_multiple_wfns_same_scope():
    """
    Test having multiple wavefunctions in the same scope.

    This ensures destructors are called correctly when multiple
    DFT wavefunctions exist simultaneously.
    """

    mol1 = psi4.geometry("""
    H
    symmetry c1
    """)

    mol2 = psi4.geometry("""
    He
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
    })

    # Create two wavefunctions
    E1, wfn1 = psi4.energy('b3lyp', molecule=mol1, return_wfn=True)
    psi4.core.clean()

    E2, wfn2 = psi4.energy('b3lyp', molecule=mol2, return_wfn=True)
    psi4.core.clean()

    # Both should have potentials
    assert wfn1.V_potential() is not None
    assert wfn2.V_potential() is not None

    # Delete in reverse order
    del wfn2
    del wfn1

    # Should not crash


def test_hf_destructor_with_nullptr():
    """
    Test that HF destructor works correctly when potential is nullptr.

    HF (not DFT) wavefunctions have V_potential() == None, so the
    destructor should handle this gracefully.
    """

    mol = psi4.geometry("""
    H 0 0 0
    H 0 0 0.74
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
    })

    E, wfn = psi4.energy('hf', return_wfn=True)

    # HF should have no potential
    assert wfn.V_potential() is None

    # Delete should not crash (destructor checks if potential != nullptr)
    del wfn

    # Run another calculation
    psi4.core.clean()
    E2 = psi4.energy('hf')

    assert psi4.compare_values(E, E2, 10, "HF energies should match")
