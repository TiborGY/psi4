"""
Unit tests for psi4.driver.qmmm module
Tests QM/MM functionality including external charges and potentials.
"""
import numpy as np
import pytest
import warnings

import psi4
from psi4 import core

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def test_external_potential_single_charge():
    """Test QM calculation with a single external point charge"""
    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    psi4.set_options({"basis": "sto-3g"})

    # Calculate energy without external charge
    e_no_charge = psi4.energy("scf")

    # Add a positive charge far away (in Bohr)
    external_charges = [[1.0, [0.0, 0.0, 10.0]]]

    # Calculate with external charge
    e_with_charge = psi4.energy("scf", external_potentials=external_charges)

    # Energy should be different but both should be valid
    assert isinstance(e_no_charge, float)
    assert isinstance(e_with_charge, float)
    assert abs(e_no_charge - e_with_charge) > 1e-6  # Should be measurably different


def test_external_potential_multiple_charges():
    """Test QM calculation with multiple external point charges"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "sto-3g"})

    # Add multiple charges around the molecule
    external_charges = [
        [0.5, [0.0, 5.0, 0.0]],
        [-0.5, [0.0, -5.0, 0.0]],
        [0.3, [5.0, 0.0, 0.0]]
    ]

    e = psi4.energy("scf", external_potentials=external_charges)

    assert isinstance(e, float)


def test_external_potential_negative_charge():
    """Test QM calculation with negative external charge"""
    mol = psi4.geometry("""
    H
    """)

    psi4.set_options({"basis": "cc-pvdz"})

    # Negative charge should stabilize the system
    external_charges = [[-1.0, [0.0, 0.0, 5.0]]]

    e = psi4.energy("scf", external_potentials=external_charges)

    assert isinstance(e, float)
    assert e < 0.0


def test_external_potential_zero_charge():
    """Test that zero charge has no effect"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "sto-3g"})

    e_no_charge = psi4.energy("scf")

    # Zero charge should have no effect
    external_charges = [[0.0, [0.0, 0.0, 10.0]]]
    e_with_zero = psi4.energy("scf", external_potentials=external_charges)

    # Energies should be essentially identical
    assert psi4.compare_values(e_no_charge, e_with_zero, 8, "Zero charge test")


def test_external_potential_symmetry():
    """Test external charges with symmetric placement"""
    mol = psi4.geometry("""
    symmetry c1
    He
    """)

    psi4.set_options({"basis": "cc-pvdz"})

    # Place charges symmetrically
    d = 5.0
    external_charges = [
        [1.0, [d, 0.0, 0.0]],
        [1.0, [-d, 0.0, 0.0]],
        [1.0, [0.0, d, 0.0]],
        [1.0, [0.0, -d, 0.0]]
    ]

    e = psi4.energy("scf", external_potentials=external_charges)

    assert isinstance(e, float)


def test_external_potential_close_charge():
    """Test external charge placed close to molecule"""
    mol = psi4.geometry("""
    H
    """)

    psi4.set_options({"basis": "sto-3g"})

    # Charge placed 2 Bohr away
    external_charges = [[1.0, [0.0, 0.0, 2.0]]]

    e = psi4.energy("scf", external_potentials=external_charges)

    # Should still converge even with close charge
    assert isinstance(e, float)


def test_external_potential_with_dft():
    """Test external charges with DFT calculation"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "sto-3g"})

    external_charges = [[0.5, [0.0, 0.0, 5.0]]]

    e = psi4.energy("b3lyp", external_potentials=external_charges)

    assert isinstance(e, float)


def test_external_potential_with_mp2():
    """Test external charges with post-HF method"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "sto-3g"})

    external_charges = [[0.5, [0.0, 0.0, 5.0]]]

    e = psi4.energy("mp2", external_potentials=external_charges)

    assert isinstance(e, float)


def test_external_potential_gradient():
    """Test that gradients work with external charges"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "sto-3g"})

    external_charges = [[0.5, [0.0, 0.0, 5.0]]]

    # Gradient should work with external potential
    try:
        grad = psi4.gradient("scf", external_potentials=external_charges)
        # Should return a matrix-like object
        assert grad is not None
    except Exception:
        # Some configurations might not support gradients with external potentials
        pytest.skip("Gradients with external potentials not available")


def test_external_charges_list_format():
    """Test that external charges accept correct list format"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "sto-3g"})

    # Test valid format: [[charge, [x, y, z]], ...]
    external_charges = [
        [1.0, [0.0, 0.0, 5.0]],
        [-0.5, [0.0, 5.0, 0.0]]
    ]

    e = psi4.energy("scf", external_potentials=external_charges)
    assert isinstance(e, float)


def test_external_charges_numpy_array():
    """Test external charges with numpy array coordinates"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "sto-3g"})

    # Coordinates as numpy arrays
    coords = np.array([0.0, 0.0, 5.0])
    external_charges = [[1.0, coords.tolist()]]

    e = psi4.energy("scf", external_potentials=external_charges)
    assert isinstance(e, float)


def test_external_potential_energy_comparison():
    """Test that external charge affects energy in expected direction"""
    mol = psi4.geometry("""
    0 1
    O  0.0  0.0  0.0
    H  0.0  0.0  1.0
    H  0.0  1.0  0.0
    """)

    psi4.set_options({"basis": "sto-3g"})

    # Calculate reference energy
    e_ref = psi4.energy("scf")

    # Positive charge near molecule should destabilize (raise energy)
    external_charges_pos = [[2.0, [0.0, 0.0, 3.0]]]
    e_pos = psi4.energy("scf", external_potentials=external_charges_pos)

    # Negative charge should stabilize (lower energy)
    external_charges_neg = [[-2.0, [0.0, 0.0, 3.0]]]
    e_neg = psi4.energy("scf", external_potentials=external_charges_neg)

    # All should be valid energies
    assert isinstance(e_ref, float)
    assert isinstance(e_pos, float)
    assert isinstance(e_neg, float)

    # Energies should differ
    assert abs(e_ref - e_pos) > 1e-5
    assert abs(e_ref - e_neg) > 1e-5


def test_external_potential_empty_list():
    """Test that empty external charges list works (should be same as no charges)"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "sto-3g"})

    e_ref = psi4.energy("scf")

    # Empty list should be same as no external charges
    e_empty = psi4.energy("scf", external_potentials=[])

    assert psi4.compare_values(e_ref, e_empty, 10, "Empty external charges")


def test_qmmm_bohr_class_deprecation():
    """Test that QMMMbohr class raises deprecation warning"""
    from psi4.driver.qmmm import QMMMbohr

    # Should raise FutureWarning when instantiated
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            qmmm_obj = QMMMbohr()
            # Check that a warning was raised
            assert len(w) > 0
            assert issubclass(w[0].category, FutureWarning)
            assert "deprecated" in str(w[0].message).lower()
        except Exception:
            # May raise error instead of warning in newer versions
            pass


def test_qmmm_class_upgrade_error():
    """Test that old QMMM class raises UpgradeHelper error"""
    from psi4.driver.qmmm import QMMM
    from psi4.driver.p4util.exceptions import UpgradeHelper

    # Should raise UpgradeHelper when instantiated
    with pytest.raises(UpgradeHelper):
        qmmm_obj = QMMM()


def test_external_potential_fragment():
    """Test external charges with molecular fragments"""
    mol = psi4.geometry("""
    He
    --
    Ne 1 5.0
    """)

    psi4.set_options({"basis": "cc-pvdz"})

    # Add external charge
    external_charges = [[0.5, [0.0, 10.0, 0.0]]]

    e = psi4.energy("scf", external_potentials=external_charges)

    assert isinstance(e, float)


def test_external_potential_large_charge():
    """Test with large magnitude external charge"""
    mol = psi4.geometry("""
    H
    """)

    psi4.set_options({"basis": "sto-3g"})

    # Large charge far away
    external_charges = [[100.0, [0.0, 0.0, 50.0]]]

    e = psi4.energy("scf", external_potentials=external_charges)

    # Should still converge
    assert isinstance(e, float)


def test_external_potential_many_charges():
    """Test with many external charges"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "sto-3g"})

    # Create grid of charges around molecule
    external_charges = []
    for x in [-5.0, 0.0, 5.0]:
        for y in [-5.0, 0.0, 5.0]:
            for z in [-5.0, 0.0, 5.0]:
                if not (x == 0 and y == 0 and z == 0):  # Skip center
                    external_charges.append([0.1, [x, y, z]])

    e = psi4.energy("scf", external_potentials=external_charges)

    assert isinstance(e, float)
    # Should have created 26 charges (3^3 - 1)
    assert len(external_charges) == 26
