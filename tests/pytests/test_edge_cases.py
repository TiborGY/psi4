"""
Unit tests for numerical edge cases and boundary conditions
Tests unusual but valid inputs, extreme values, and numerical stability.
"""
import numpy as np
import pytest

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]


def test_single_atom_helium():
    """Test single atom calculation - simplest possible system"""
    mol = psi4.geometry("""
    He
    """)

    psi4.set_options({"basis": "cc-pvdz"})
    e = psi4.energy("scf")

    assert isinstance(e, float)
    assert e < 0.0  # Energy should be negative
    assert psi4.compare_values(-2.855, e, 2, "He SCF energy")


def test_single_atom_neon():
    """Test single atom with more electrons"""
    mol = psi4.geometry("""
    Ne
    """)

    psi4.set_options({"basis": "cc-pvdz"})
    e = psi4.energy("scf")

    assert isinstance(e, float)
    assert e < -100.0  # Neon energy should be significantly negative


def test_h2_very_short_bond():
    """Test H2 with very short bond distance (repulsive region)"""
    mol = psi4.geometry("""
    H
    H 1 0.3
    """)

    psi4.set_options({"basis": "cc-pvdz"})
    e = psi4.energy("scf")

    # Very short bond should have high energy due to repulsion
    assert isinstance(e, float)
    # Energy will be higher than equilibrium but should still complete


def test_h2_very_long_bond():
    """Test H2 with very long bond distance (dissociation limit)"""
    mol = psi4.geometry("""
    H
    H 1 10.0
    """)

    psi4.set_options({"basis": "cc-pvdz"})
    e = psi4.energy("scf")

    # At large separation, should approach 2x H atom energy
    assert isinstance(e, float)
    # Should be close to -1.0 hartree (2 separate H atoms)
    assert -1.1 < e < -0.9


def test_tight_convergence():
    """Test SCF with very tight convergence thresholds"""
    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    psi4.set_options({
        "basis": "cc-pvdz",
        "e_convergence": 1e-10,
        "d_convergence": 1e-10
    })

    e = psi4.energy("scf")

    assert isinstance(e, float)
    assert psi4.compare_values(-76.026, e, 3, "Tight convergence SCF")


def test_loose_convergence():
    """Test SCF with loose convergence thresholds"""
    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    psi4.set_options({
        "basis": "cc-pvdz",
        "e_convergence": 1e-4,
        "d_convergence": 1e-4
    })

    e = psi4.energy("scf")

    assert isinstance(e, float)
    # Should converge quickly with loose thresholds
    assert psi4.compare_values(-76.026, e, 2, "Loose convergence SCF")


def test_high_spin_multiplicity():
    """Test high-spin open shell system"""
    mol = psi4.geometry("""
    0 5
    C
    """)

    psi4.set_options({
        "basis": "cc-pvdz",
        "reference": "uhf"
    })

    e = psi4.energy("scf")

    assert isinstance(e, float)
    # Quintet C atom
    assert e < 0.0


def test_minimal_basis_large_atom():
    """Test large atom with minimal basis set"""
    mol = psi4.geometry("""
    Ar
    """)

    psi4.set_options({
        "basis": "sto-3g",
        "reference": "rhf"
    })

    e = psi4.energy("scf")

    assert isinstance(e, float)
    assert e < -500.0  # Argon has 18 electrons


@pytest.mark.long
def test_very_diffuse_basis():
    """Test molecule with highly diffuse basis functions"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "aug-cc-pv5z"})
    e = psi4.energy("scf")

    assert isinstance(e, float)
    assert psi4.compare_values(-1.133, e, 2, "H2 with very diffuse basis")


def test_near_linear_molecule():
    """Test molecule with nearly linear geometry"""
    mol = psi4.geometry("""
    O
    H 1 1.0
    H 1 1.0 2 179.0
    """)

    psi4.set_options({"basis": "sto-3g"})
    e = psi4.energy("scf")

    # Near-linear water should still work
    assert isinstance(e, float)


def test_linear_molecule_symmetry():
    """Test linear molecule symmetry detection"""
    mol = psi4.geometry("""
    H 0.0 0.0 0.0
    C 0.0 0.0 1.0
    C 0.0 0.0 2.0
    H 0.0 0.0 3.0
    """)

    psi4.set_options({"basis": "sto-3g"})

    # Linear molecule - symmetry should be detected
    pg = mol.schoenflies_symbol()
    assert isinstance(pg, str)

    e = psi4.energy("scf")
    assert isinstance(e, float)


def test_cation_radical():
    """Test cation radical (odd electron, positive charge)"""
    mol = psi4.geometry("""
    1 2
    H
    H 1 1.0
    """)

    psi4.set_options({
        "basis": "cc-pvdz",
        "reference": "uhf"
    })

    e = psi4.energy("scf")

    assert isinstance(e, float)
    # H2+ should have energy around -0.6 hartree at this geometry


def test_anion():
    """Test anion calculation"""
    mol = psi4.geometry("""
    -1 1
    H
    """)

    psi4.set_options({
        "basis": "aug-cc-pvdz",  # Need diffuse functions for anions
        "reference": "uhf"
    })

    e = psi4.energy("scf")

    assert isinstance(e, float)
    # H- should have negative energy


def test_dication():
    """Test doubly charged cation"""
    mol = psi4.geometry("""
    2 1
    O
    """)

    psi4.set_options({
        "basis": "cc-pvdz",
        "reference": "rhf"
    })

    e = psi4.energy("scf")

    assert isinstance(e, float)
    # O2+ has 6 electrons


def test_high_angular_momentum_basis():
    """Test basis set with high angular momentum functions"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    # cc-pV5Z includes up to h functions
    psi4.set_options({"basis": "cc-pv5z"})
    e = psi4.energy("scf")

    assert isinstance(e, float)


def test_mixed_charge_fragments():
    """Test system with mixed charge fragments"""
    mol = psi4.geometry("""
    1 1
    H
    --
    -1 1
    H
    """)

    psi4.set_options({
        "basis": "aug-cc-pvdz",
        "reference": "rhf"
    })

    # Net charge should be 0
    assert mol.molecular_charge() == 0
    e = psi4.energy("scf")
    assert isinstance(e, float)


def test_heavy_atom_hydrogen():
    """Test molecule with heavy atom and hydrogen"""
    mol = psi4.geometry("""
    Cl
    H 1 1.27
    """)

    psi4.set_options({"basis": "cc-pvdz"})
    e = psi4.energy("scf")

    assert isinstance(e, float)
    assert e < -450.0  # Chlorine has many electrons


def test_symmetry_c1():
    """Test calculation in C1 symmetry (no symmetry)"""
    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)

    psi4.set_options({"basis": "cc-pvdz"})
    e = psi4.energy("scf")

    assert mol.schoenflies_symbol() == "c1"
    assert isinstance(e, float)


def test_uhf_closed_shell():
    """Test UHF on closed shell system (should work but inefficient)"""
    mol = psi4.geometry("""
    0 1
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    psi4.set_options({
        "basis": "sto-3g",
        "reference": "uhf"
    })

    e = psi4.energy("scf")

    # UHF on closed shell should give same energy as RHF
    assert isinstance(e, float)


def test_small_basis_correlation():
    """Test correlation method with minimal basis"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "sto-3g"})
    e = psi4.energy("mp2")

    assert isinstance(e, float)
    # Should be close to STO-3G SCF energy with small correlation correction


def test_dft_minimal_system():
    """Test DFT on minimal system"""
    mol = psi4.geometry("""
    He
    """)

    psi4.set_options({"basis": "cc-pvdz"})
    e = psi4.energy("b3lyp")

    assert isinstance(e, float)
    assert e < 0.0


def test_single_point_at_equilibrium():
    """Test that energy at known equilibrium geometry is reasonable"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "cc-pvdz"})
    e = psi4.energy("scf")

    # H2 equilibrium SCF energy should be around -1.13 hartree
    assert psi4.compare_values(-1.133, e, 2, "H2 equilibrium energy")


def test_fragment_with_one_atom():
    """Test fragment calculation with single atom fragment"""
    mol = psi4.geometry("""
    He
    --
    Ne 1 5.0
    """)

    psi4.set_options({"basis": "cc-pvdz"})

    assert mol.nfragments() == 2
    e = psi4.energy("scf")
    assert isinstance(e, float)


def test_three_fragments():
    """Test system with three molecular fragments"""
    mol = psi4.geometry("""
    He
    --
    He 1 5.0
    --
    He 1 10.0 2 90.0
    """)

    psi4.set_options({"basis": "cc-pvdz"})

    assert mol.nfragments() == 3
    e = psi4.energy("scf")
    assert isinstance(e, float)


def test_numerical_stability_small_molecule():
    """Test numerical stability with small well-behaved molecule"""
    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    psi4.set_options({
        "basis": "cc-pvdz",
        "scf_type": "pk"
    })

    e1 = psi4.energy("scf")

    # Run again - should get same result
    e2 = psi4.energy("scf")

    # Results should be identical to machine precision
    assert abs(e1 - e2) < 1e-10


def test_different_scf_types_agreement():
    """Test that different SCF algorithms give same result"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "cc-pvdz"})

    psi4.set_options({"scf_type": "pk"})
    e_pk = psi4.energy("scf")

    psi4.set_options({"scf_type": "direct"})
    e_direct = psi4.energy("scf")

    psi4.set_options({"scf_type": "df"})
    e_df = psi4.energy("scf")

    # PK and DIRECT should be identical
    assert psi4.compare_values(e_pk, e_direct, 10, "PK vs DIRECT")

    # DF should be very close
    assert psi4.compare_values(e_pk, e_df, 5, "PK vs DF")


def test_zero_nuclear_repulsion():
    """Test single atom (zero nuclear repulsion)"""
    mol = psi4.geometry("""
    H
    """)

    nre = mol.nuclear_repulsion_energy()
    assert nre == 0.0

    psi4.set_options({"basis": "cc-pvdz"})
    e = psi4.energy("scf")

    assert isinstance(e, float)
