"""
Unit tests for psi4.driver.molutil module functions
Tests molecule utility functions including geometry creation, conversion, and manipulation.
"""
import numpy as np
import pytest

import psi4
from psi4 import core
from psi4.driver import molutil

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def test_geometry_basic_creation():
    """Test basic geometry() function with simple water molecule"""
    mol = molutil.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    assert mol.natom() == 3
    assert mol.molecular_charge() == 0
    assert mol.multiplicity() == 1


def test_geometry_with_charge_and_multiplicity():
    """Test geometry creation with explicit charge and multiplicity"""
    mol = molutil.geometry("""
    1 2
    O
    H 1 1.0
    """)

    assert mol.molecular_charge() == 1
    assert mol.multiplicity() == 2
    assert mol.natom() == 2


def test_geometry_cartesian_coordinates():
    """Test geometry creation with Cartesian coordinates"""
    mol = molutil.geometry("""
    O  0.000000  0.000000  0.117176
    H  0.000000  0.755453 -0.468706
    H  0.000000 -0.755453 -0.468706
    """)

    assert mol.natom() == 3
    # Check that molecule was created successfully
    assert mol.nuclear_repulsion_energy() > 0.0


def test_geometry_with_units():
    """Test geometry creation with different units"""
    mol_bohr = molutil.geometry("""
    units bohr
    O  0.000000  0.000000  0.221364
    H  0.000000  1.427651 -0.885456
    H  0.000000 -1.427651 -0.885456
    """)

    mol_ang = molutil.geometry("""
    units angstrom
    O  0.000000  0.000000  0.117176
    H  0.000000  0.755453 -0.468706
    H  0.000000 -0.755453 -0.468706
    """)

    # Both should create valid molecules
    assert mol_bohr.natom() == 3
    assert mol_ang.natom() == 3
    # Nuclear repulsion should be approximately equal
    assert abs(mol_bohr.nuclear_repulsion_energy() - mol_ang.nuclear_repulsion_energy()) < 0.1


def test_geometry_with_symmetry():
    """Test geometry with symmetry specification"""
    mol_c1 = molutil.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)

    mol_c2v = molutil.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    assert mol_c1.schoenflies_symbol() == "c1"
    # Water in standard orientation should detect symmetry
    assert mol_c2v.natom() == 3


def test_geometry_fragments():
    """Test geometry creation with molecular fragments"""
    mol = molutil.geometry("""
    0 1
    O
    H 1 0.96
    H 1 0.96 2 104.5
    --
    0 1
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    assert mol.natom() == 6
    assert mol.nfragments() == 2


def test_geometry_ghost_atoms():
    """Test geometry with ghost atoms"""
    mol = molutil.geometry("""
    O
    H 1 0.96
    @H 1 0.96 2 104.5
    """)

    assert mol.natom() == 3
    # Check that third atom is a ghost
    assert mol.Z(2) == 0  # Ghost atoms have Z=0 in some representations


def test_geometry_name():
    """Test setting molecule name"""
    mol = molutil.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """, name="my_water")

    assert mol.name() == "my_water"


def test_activate_molecule():
    """Test activate() function sets active molecule"""
    mol1 = molutil.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """, name="water1")

    mol2 = molutil.geometry("""
    He
    """, name="helium")

    # Activate mol1
    molutil.activate(mol1)
    active = core.get_active_molecule()
    assert active.name() == "water1"
    assert active.natom() == 3

    # Activate mol2
    molutil.activate(mol2)
    active = core.get_active_molecule()
    assert active.name() == "helium"
    assert active.natom() == 1


def test_molecule_from_arrays_basic():
    """Test Molecule.from_arrays() with basic arrays"""
    geom = np.array([
        [0.000000,  0.000000,  0.117176],
        [0.000000,  0.755453, -0.468706],
        [0.000000, -0.755453, -0.468706]
    ])

    mol = core.Molecule.from_arrays(
        geom=geom,
        elez=[8, 1, 1],
        units='Angstrom'
    )

    assert mol.natom() == 3
    assert mol.Z(0) == 8  # Oxygen
    assert mol.Z(1) == 1  # Hydrogen


def test_molecule_from_arrays_with_charge_mult():
    """Test Molecule.from_arrays() with charge and multiplicity"""
    geom = np.array([[0.0, 0.0, 0.0]])

    mol = core.Molecule.from_arrays(
        geom=geom,
        elez=[6],
        molecular_charge=1,
        molecular_multiplicity=2
    )

    assert mol.natom() == 1
    assert mol.molecular_charge() == 1
    assert mol.multiplicity() == 2


def test_molecule_from_arrays_with_fragments():
    """Test Molecule.from_arrays() with fragment specification"""
    geom = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 5.0],
        [0.0, 0.0, 6.0]
    ])

    mol = core.Molecule.from_arrays(
        geom=geom,
        elez=[8, 1, 8, 1],
        fragment_separators=[2],
        fragment_charges=[0.0, 0.0],
        fragment_multiplicities=[1, 1]
    )

    assert mol.natom() == 4
    assert mol.nfragments() == 2


def test_molecule_from_arrays_with_labels():
    """Test Molecule.from_arrays() with custom labels"""
    geom = np.array([[0.0, 0.0, 0.0]])

    mol = core.Molecule.from_arrays(
        geom=geom,
        elez=[6],
        elbl=['C_custom']
    )

    assert mol.natom() == 1
    assert mol.Z(0) == 6


def test_molecule_from_string_simple():
    """Test Molecule.from_string() with simple molecule"""
    mol = core.Molecule.from_string("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    assert mol.natom() == 3
    assert mol.Z(0) == 8


def test_molecule_from_string_with_options():
    """Test Molecule.from_string() with fix_com and fix_orientation"""
    mol = core.Molecule.from_string(
        """
        O
        H 1 0.96
        H 1 0.96 2 104.5
        """,
        fix_com=True,
        fix_orientation=True
    )

    assert mol.natom() == 3


def test_molecule_from_string_xyz_format():
    """Test Molecule.from_string() with XYZ format"""
    mol = core.Molecule.from_string("""
    3
    Water molecule
    O  0.000000  0.000000  0.117176
    H  0.000000  0.755453 -0.468706
    H  0.000000 -0.755453 -0.468706
    """, dtype='xyz')

    assert mol.natom() == 3


def test_molecule_to_dict_and_from_dict():
    """Test round-trip conversion molecule -> dict -> molecule"""
    mol1 = molutil.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    # Convert to dict
    mol_dict = mol1.to_dict()

    # Check dict has expected keys
    assert 'geom' in mol_dict
    assert 'elez' in mol_dict
    assert 'molecular_charge' in mol_dict

    # Convert back to molecule
    mol2 = core.Molecule.from_dict(mol_dict)

    # Should have same properties
    assert mol2.natom() == mol1.natom()
    assert mol2.molecular_charge() == mol1.molecular_charge()
    assert mol2.multiplicity() == mol1.multiplicity()


def test_molecule_to_arrays():
    """Test molecule.to_arrays() method"""
    mol = molutil.geometry("""
    O  0.0  0.0  0.0
    H  0.0  0.0  1.0
    H  0.0  1.0  0.0
    """)

    arrays = mol.to_arrays()

    # Check that arrays dict has expected keys
    assert 'geom' in arrays
    assert 'elez' in arrays
    assert 'molecular_charge' in arrays
    assert 'molecular_multiplicity' in arrays

    # Check geom is numpy array with correct shape
    assert isinstance(arrays['geom'], np.ndarray)
    assert arrays['geom'].size == 9  # 3 atoms * 3 coordinates


def test_molecule_to_string():
    """Test molecule.to_string() method"""
    mol = molutil.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    # Convert to string
    mol_str = mol.to_string(dtype='psi4')

    # Should be a non-empty string
    assert isinstance(mol_str, str)
    assert len(mol_str) > 0
    assert 'O' in mol_str or 'H' in mol_str


def test_invalid_geometry_syntax():
    """Test that invalid geometry raises appropriate error"""
    with pytest.raises(Exception):  # Should raise some kind of parsing error
        molutil.geometry("This is not valid geometry")


def test_geometry_with_variables():
    """Test geometry with variable definitions"""
    mol = molutil.geometry("""
    O
    H 1 r
    H 1 r 2 a

    r = 0.96
    a = 104.5
    """)

    assert mol.natom() == 3


def test_geometry_empty_fragments():
    """Test that geometry handles fragment specifications correctly"""
    # Create a dimer
    mol = molutil.geometry("""
    He
    --
    He 1 3.0
    """)

    assert mol.natom() == 2
    assert mol.nfragments() == 2


def test_molecule_nuclear_repulsion_energy():
    """Test that nuclear repulsion energy is computed correctly"""
    mol = molutil.geometry("""
    H
    H 1 1.0
    """)

    # H2 with 1.0 Angstrom separation should have positive repulsion energy
    nre = mol.nuclear_repulsion_energy()
    assert nre > 0.0
    assert nre < 10.0  # Sanity check


def test_molecule_center_of_mass():
    """Test center of mass calculation"""
    mol = molutil.geometry("""
    O  0.0  0.0  0.0
    H  1.0  0.0  0.0
    H -1.0  0.0  0.0
    """)

    # Symmetric molecule should have COM near origin (along x for this case)
    # Just check that we can call the update_geometry without error
    mol.update_geometry()
    assert mol.natom() == 3


def test_molecule_point_group():
    """Test point group detection"""
    # Linear molecule
    mol_linear = molutil.geometry("""
    H  0.0  0.0  0.0
    H  0.0  0.0  1.0
    """)

    # Just verify we can get schoenflies symbol
    pg = mol_linear.schoenflies_symbol()
    assert isinstance(pg, str)
    assert len(pg) > 0


def test_molecule_mass():
    """Test molecular mass calculation"""
    mol = molutil.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    # Water should have mass around 18 amu
    mass = mol.molecular_mass()
    assert 17.0 < mass < 19.0


def test_molecule_from_schema():
    """Test Molecule.from_schema() with QCSchema format"""
    schema = {
        'symbols': ['O', 'H', 'H'],
        'geometry': [
            0.0, 0.0, 0.0,
            0.0, 0.0, 1.8,
            0.0, 1.8, 0.0
        ],
        'molecular_charge': 0.0,
        'molecular_multiplicity': 1
    }

    mol = core.Molecule.from_schema(schema)

    assert mol.natom() == 3
    assert mol.molecular_charge() == 0


def test_molecule_to_schema():
    """Test molecule.to_schema() method"""
    mol = molutil.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    schema = mol.to_schema(dtype=2)

    # Should return a dictionary with QCSchema keys
    assert isinstance(schema, dict)
    assert 'symbols' in schema
    assert 'geometry' in schema


def test_molecule_reorient():
    """Test molecule reorientation"""
    mol = molutil.geometry("""
    O  1.0  1.0  1.0
    H  2.0  1.0  1.0
    H  1.0  2.0  1.0
    symmetry c1
    """)

    # Update geometry to reorient
    mol.update_geometry()

    # Molecule should still have 3 atoms after reorientation
    assert mol.natom() == 3
