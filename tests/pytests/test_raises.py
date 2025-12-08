import pytest
import psi4


pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

def test_dft_grid_threaded_raise():
    dimer = psi4.geometry("""
      1 1
      K -4.067042 -1.894214 0.002270
    """)
    
    psi4.set_options({
        "dft_grid_name": "SG1",
        "dft_vv10_radial_points": 50,
        "dft_vv10_spherical_points": 194,
        "dft_nuclear_scheme": "treutler",
        "dft_radial_scheme": "EM",
        "basis": "def2-TZVPPD",
    })
    
    with pytest.raises(RuntimeError) as e:
        ene = psi4.energy("wB97M-V")

    assert "There is no SG-1 grid defined for the requested atomic number" in str(e.value)

def test_cc_uhf_raise1():
    psi4.geometry("""
    0 1
    O
    H 1 1.0
    H 1 1.0 2 104.5
    """)

    psi4.set_options({"reference": "uhf"})
    wfn = psi4.energy("scf/sto-3g", return_wfn=True)[1]
    psi4.set_options({"reference": "rhf"})
    with pytest.raises(psi4.ValidationError) as e:
        psi4.properties("ccsd/sto-3g", properties=['polarizability'], ref_wfn=wfn)

    assert "Non-RHF CC response properties are not implemented." in str(e.value)

def test_cc_uhf_raise2():
    psi4.geometry("""
    0 1
    O
    H 1 1.0
    H 1 1.0 2 104.5
    """)

    psi4.set_options({"reference": "uhf"})
    with pytest.raises(psi4.ValidationError) as e:
        psi4.properties("ccsd/sto-3g", properties=['polarizability'])

    assert "Non-RHF CC response properties are not implemented." in str(e.value)


def test_casscf_invalid_active_irrep():
    """Test that CASSCF raises a clear error when ACTIVE array requests orbitals
    in an irrep with no available molecular orbitals.

    This is a regression test for issue #3096. Before the fix, this input would
    cause a cryptic "DSYEV diagonalizer failed" error deep in the DETCI code.
    After the fix, it raises an InputException with a clear message about which
    irrep has insufficient orbitals.
    """
    psi4.geometry("""
    0 1
    H
    H 1 0.74
    symmetry d2h
    """)

    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "rhf",
        "docc": [1, 0, 0, 0, 0, 0, 0, 0],
        # Invalid: irrep 5 (B2u in D2h) has 0 orbitals for H2/cc-pVDZ,
        # but we request 1 active orbital there
        "active": [1, 0, 0, 0, 0, 1, 0, 0],
    })

    with pytest.raises(RuntimeError) as e:
        psi4.energy("casscf")

    # Verify the error message mentions the problematic irrep and that it
    # exceeds available orbitals (the specific irrep number may vary based on
    # how D2h irreps are ordered, but the message pattern should match)
    assert "exceeds available orbitals" in str(e.value)
    assert "ACTIVE" in str(e.value)


def test_casscf_valid_active_succeeds():
    """Test that a valid CASSCF calculation with proper ACTIVE array succeeds.

    This complements test_casscf_invalid_active_irrep by verifying that a
    correctly specified active space works as expected.
    """
    psi4.geometry("""
    0 1
    H
    H 1 0.74
    symmetry d2h
    """)

    psi4.set_options({
        "basis": "sto-3g",  # Smaller basis for faster test
        "reference": "rhf",
        "docc": [1, 0, 0, 0, 0, 0, 0, 0],
        # Valid: only request active orbitals in irreps that have orbitals
        # For H2 in D2h with STO-3G: Ag has 1 orbital, B1u has 1 orbital
        "active": [1, 0, 0, 0, 0, 0, 0, 1],
    })

    # This should complete without raising an exception
    energy = psi4.energy("casscf")

    # Verify we got a reasonable energy (H2 ground state should be around -1 Hartree)
    assert energy < 0.0
    assert energy > -2.0

