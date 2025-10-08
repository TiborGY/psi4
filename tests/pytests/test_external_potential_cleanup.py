"""
Test to verify HF destructor fix is safe with external potentials.

External potentials (point charges, electric fields, PCM, etc.) are stored
separately from the VBase DFT potential. This test verifies that cleaning up
the DFT grid in the destructor doesn't interfere with external potentials.
"""

import pytest
import numpy as np
import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def test_dft_with_external_charges():
    """
    Test DFT calculation with external point charges.

    External charges are stored in external_pot_ (separate from VBase potential).
    Verifies destructor cleanup doesn't interfere.
    """

    mol = psi4.geometry("""
    0 1
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)

    # Define external point charges (e.g., from MM environment)
    # Positive charge at (3, 0, 0) to perturb the wavefunction
    external_charges = [(1.0, 3.0, 0.0, 0.0)]  # (charge, x, y, z) in Bohr

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
    })

    # Run DFT with external potential
    E, wfn = psi4.energy('b3lyp', external_potentials=external_charges, return_wfn=True)

    # Verify both potentials exist
    assert wfn.V_potential() is not None, "DFT potential should exist"
    assert wfn.external_pot() is not None, "External potential should exist"

    # Get the external potential charges
    ext_pot = wfn.external_pot()
    charges = ext_pot.getCharges()
    assert len(charges) == 1
    assert charges[0][0] == 1.0  # Check charge value

    # Delete wavefunction - should clean up DFT grid but preserve external pot until last reference
    energy_ref = E
    del wfn

    # Run another calculation to verify everything still works
    psi4.core.clean()
    E2 = psi4.energy('b3lyp', external_potentials=external_charges)

    # Energy should be consistent
    assert psi4.compare_values(energy_ref, E2, 7, "B3LYP energy with external charges")


def test_dft_repeated_with_external_charges():
    """
    Test repeated DFT calculations with external charges (memory leak scenario).

    Each calculation creates a new wavefunction with:
    - Its own VBase potential (grid)
    - Its own external potential object
    Both should be cleaned up properly.
    """

    mol = psi4.geometry("""
    Ne
    symmetry c1
    """)

    external_charges = [(0.5, 2.0, 0.0, 0.0), (-0.5, -2.0, 0.0, 0.0)]

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
    })

    energies = []
    for i in range(5):
        E = psi4.energy('pbe', external_potentials=external_charges)
        energies.append(E)
        psi4.core.clean()

    # All energies should be identical
    for i, E in enumerate(energies[1:], 1):
        assert psi4.compare_values(energies[0], E, 9,
                                   f"PBE energy with external charges, iteration {i}")


def test_hf_with_external_charges():
    """
    Test pure HF with external charges (no DFT potential).

    HF has no VBase potential, only external potential.
    Verifies destructor handles nullptr case correctly.
    """

    mol = psi4.geometry("""
    H 0 0 0
    H 0 0 0.74
    symmetry c1
    """)

    external_charges = [(0.1, 0.0, 0.0, 3.0)]

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
    })

    E, wfn = psi4.energy('hf', external_potentials=external_charges, return_wfn=True)

    # HF should have no DFT potential
    assert wfn.V_potential() is None, "HF should not have DFT potential"
    # But should have external potential
    assert wfn.external_pot() is not None, "Should have external potential"

    # Delete wavefunction - destructor should skip VBase cleanup (nullptr)
    del wfn

    # Run another calculation
    psi4.core.clean()
    E2 = psi4.energy('hf', external_potentials=external_charges)

    assert psi4.compare_values(E, E2, 10, "HF energy with external charges")


def test_dft_with_perturb_h():
    """
    Test DFT with PERTURB_H option (electric field perturbation).

    PERTURB_H adds a dipole perturbation to the Hamiltonian.
    This is stored in external_potentials_ (separate from VBase).
    """

    mol = psi4.geometry("""
    H 0 0 0
    F 0 0 0.92
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
        'perturb_h': True,
        'perturb_with': 'dipole',
        'perturb_dipole': [0.0, 0.0, 0.1],  # Z-direction electric field
    })

    E, wfn = psi4.energy('b3lyp', return_wfn=True)

    # Should have DFT potential
    assert wfn.V_potential() is not None

    # Delete and re-run
    del wfn
    psi4.core.clean()

    E2 = psi4.energy('b3lyp')
    assert psi4.compare_values(E, E2, 7, "B3LYP energy with PERTURB_H")

    # Reset perturb_h for other tests
    psi4.set_options({'perturb_h': False})


def test_findif_gradient_with_external_field():
    """
    Test finite difference gradients with external potential.

    This combines:
    - Multiple energy calculations (finite differences)
    - External potentials
    - DFT grid cleanup

    Each displacement should have independent DFT grid and external potential.
    """

    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)

    external_charges = [(0.2, 3.0, 0.0, 0.0)]

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
        'points': 3,
    })

    # Calculate gradient with external field using finite differences
    # Each displacement will create its own wavefunction with DFT grid + external pot
    grad = psi4.gradient('pbe', dertype='findif', external_potentials=external_charges)

    assert grad is not None
    assert grad.rows() == 3  # 3 atoms


def test_external_potential_independence():
    """
    Verify that VBase potential and external potential are truly independent.

    Create multiple wavefunctions with different combinations and verify
    cleanup works correctly in all cases.
    """

    mol = psi4.geometry("""
    He
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'df',
    })

    # Case 1: DFT only (has VBase, no external)
    E1, wfn1 = psi4.energy('b3lyp', return_wfn=True)
    assert wfn1.V_potential() is not None
    assert wfn1.external_pot() is None
    del wfn1

    psi4.core.clean()

    # Case 2: HF only (no VBase, no external)
    E2, wfn2 = psi4.energy('hf', return_wfn=True)
    assert wfn2.V_potential() is None
    assert wfn2.external_pot() is None
    del wfn2

    psi4.core.clean()

    # Case 3: DFT with external (has VBase, has external)
    external = [(0.1, 2.0, 0.0, 0.0)]
    E3, wfn3 = psi4.energy('b3lyp', external_potentials=external, return_wfn=True)
    assert wfn3.V_potential() is not None
    assert wfn3.external_pot() is not None
    del wfn3

    psi4.core.clean()

    # Case 4: HF with external (no VBase, has external)
    E4, wfn4 = psi4.energy('hf', external_potentials=external, return_wfn=True)
    assert wfn4.V_potential() is None
    assert wfn4.external_pot() is not None
    del wfn4

    # All should complete without crashes
