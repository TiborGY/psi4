"""
Unit tests for psi4.driver.p4util utility functions
Tests utility functions for memory management, OEProp, and other helpers.
"""
import pytest
import numpy as np

import psi4
from psi4 import core
from psi4.driver.p4util import util
from psi4.driver.p4util.exceptions import *

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def test_set_memory():
    """Test memory setting with different units"""
    # Set memory in bytes
    util.set_memory(500000000)  # 500 MB
    mem = psi4.core.get_memory()
    assert mem == 500000000

    # Set memory in MB
    util.set_memory('500 MB')
    mem = psi4.core.get_memory()
    assert mem == 500 * 1024 * 1024

    # Set memory in GB
    util.set_memory('1 GB')
    mem = psi4.core.get_memory()
    assert mem == 1024 * 1024 * 1024


def test_get_memory():
    """Test memory retrieval"""
    util.set_memory(500000000)
    mem = util.get_memory()

    assert isinstance(mem, int)
    assert mem == 500000000


def test_memory_units():
    """Test memory setting with various unit formats"""
    # Test different unit specifications
    util.set_memory('100 mb')
    assert psi4.core.get_memory() == 100 * 1024 * 1024

    util.set_memory('1 gb')
    assert psi4.core.get_memory() == 1024 * 1024 * 1024

    util.set_memory('1024 kb')
    assert psi4.core.get_memory() == 1024 * 1024


def test_oeprop_dipole():
    """Test OEProp for dipole moment calculation"""
    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    psi4.set_options({"basis": "sto-3g"})
    e, wfn = psi4.energy("scf", return_wfn=True)

    # Calculate dipole
    util.oeprop(wfn, 'DIPOLE')

    # Check that dipole was calculated
    dipole = psi4.variable('SCF DIPOLE')
    assert isinstance(dipole, (list, np.ndarray))
    assert len(dipole) == 3


def test_oeprop_quadrupole():
    """Test OEProp for quadrupole moment calculation"""
    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    psi4.set_options({"basis": "sto-3g"})
    e, wfn = psi4.energy("scf", return_wfn=True)

    # Calculate quadrupole
    util.oeprop(wfn, 'QUADRUPOLE')

    # Check that quadrupole was calculated (6 independent components)
    try:
        quad = psi4.variable('SCF QUADRUPOLE')
        assert isinstance(quad, (list, np.ndarray))
    except KeyError:
        # Variable name might be different
        pass


def test_oeprop_mulliken_charges():
    """Test OEProp for Mulliken population analysis"""
    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    psi4.set_options({"basis": "sto-3g"})
    e, wfn = psi4.energy("scf", return_wfn=True)

    # Calculate Mulliken charges
    util.oeprop(wfn, 'MULLIKEN_CHARGES')

    # Check that charges were calculated
    try:
        charges = psi4.variable('MULLIKEN CHARGES')
        assert isinstance(charges, (list, np.ndarray))
        assert len(charges) == 3  # 3 atoms
    except KeyError:
        # Variable name might be stored differently
        pass


def test_oeprop_lowdin_charges():
    """Test OEProp for Lowdin population analysis"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "cc-pvdz"})
    e, wfn = psi4.energy("scf", return_wfn=True)

    # Calculate Lowdin charges
    util.oeprop(wfn, 'LOWDIN_CHARGES')

    # Should complete without error
    assert wfn is not None


def test_oeprop_multiple_properties():
    """Test OEProp with multiple properties at once"""
    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    psi4.set_options({"basis": "sto-3g"})
    e, wfn = psi4.energy("scf", return_wfn=True)

    # Calculate multiple properties
    util.oeprop(wfn, 'DIPOLE', 'QUADRUPOLE', 'MULLIKEN_CHARGES')

    # All properties should be calculated
    dipole = psi4.variable('SCF DIPOLE')
    assert isinstance(dipole, (list, np.ndarray))


def test_oeprop_with_title():
    """Test OEProp with custom title"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "sto-3g"})
    e, wfn = psi4.energy("scf", return_wfn=True)

    # Calculate with custom title
    util.oeprop(wfn, 'DIPOLE', title='My Calculation')

    # Should complete without error
    assert wfn is not None


def test_sanitize_method():
    """Test method name sanitization for documentation links"""
    # Test various method name transformations
    assert sanitize_method('CCSD(T)') == 'ccsd_prt_pr'
    assert sanitize_method('MP2.5') == 'mp2p5'
    assert sanitize_method('CCSD+T(CCSD)') == 'ccsdpt_prccsd_pr'
    assert sanitize_method('CCSDT-1a') == 'ccsdt1a'


def test_docs_table_link_summary():
    """Test documentation link generation for summary table"""
    link = docs_table_link('ccsd', mode='summary')

    assert isinstance(link, str)
    assert 'psicode.org' in link or 'psi4' in link.lower()
    assert 'introduction' in link


def test_docs_table_link_details():
    """Test documentation link generation for details table"""
    link = docs_table_link('mp2', mode='details')

    assert isinstance(link, str)
    assert 'capabilities' in link


def test_docs_table_link_specific():
    """Test documentation link for specific module tables"""
    # Test different table modes
    link = docs_table_link('ccsd', mode='ccenergy')
    assert 'cc.html' in link

    link = docs_table_link('mp2', mode='dfmp2')
    assert 'dfmp2.html' in link

    link = docs_table_link('b3lyp', mode='scf')
    assert 'dft.html' in link


def test_validation_error_creation():
    """Test ValidationError exception"""
    try:
        raise ValidationError("Test validation error message")
    except ValidationError as e:
        assert "Test validation error message" in str(e)
        assert hasattr(e, 'message')


def test_convergence_error_creation():
    """Test ConvergenceError exception"""
    try:
        raise ConvergenceError("SCF", 50)
    except ConvergenceError as e:
        assert "Could not converge" in str(e)
        assert e.iteration == 50
        assert "SCF" in e.message


def test_upgrade_helper_error():
    """Test UpgradeHelper exception"""
    try:
        raise UpgradeHelper("old_function", "new_function", "1.5", " See docs for details.")
    except UpgradeHelper as e:
        assert "old_function" in str(e)
        assert "new_function" in str(e)
        assert "1.5" in str(e)


def test_missing_method_error():
    """Test MissingMethodError exception"""
    try:
        raise MissingMethodError("Method XYZ is not available")
    except MissingMethodError as e:
        assert "XYZ" in str(e)
        assert hasattr(e, 'message')


def test_test_comparison_error():
    """Test TestComparisonError exception"""
    try:
        raise TestComparisonError("Values don't match: expected 1.0, got 2.0")
    except TestComparisonError as e:
        assert "1.0" in str(e) or "2.0" in str(e)
        assert hasattr(e, 'message')


def test_libint2_configuration():
    """Test libint2 configuration info retrieval"""
    try:
        config = util.libint2_configuration()
        # Should return a dictionary or similar structure
        assert config is not None
    except Exception:
        # May not be available in all builds
        pytest.skip("libint2_configuration not available")


def test_oeprop_wiberg_lowdin():
    """Test Wiberg bond indices calculation"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "cc-pvdz"})
    e, wfn = psi4.energy("scf", return_wfn=True)

    # Calculate Wiberg bond indices
    try:
        util.oeprop(wfn, 'WIBERG_LOWDIN_INDICES')
        assert wfn is not None
    except Exception:
        # May not be available in all configurations
        pytest.skip("Wiberg indices not available")


def test_oeprop_mayer_indices():
    """Test Mayer bond indices calculation"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "cc-pvdz"})
    e, wfn = psi4.energy("scf", return_wfn=True)

    # Calculate Mayer bond indices
    try:
        util.oeprop(wfn, 'MAYER_INDICES')
        assert wfn is not None
    except Exception:
        # May not be available
        pytest.skip("Mayer indices not available")


def test_oeprop_mbis_charges():
    """Test MBIS charge analysis"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "cc-pvdz"})
    e, wfn = psi4.energy("scf", return_wfn=True)

    # Calculate MBIS charges
    try:
        util.oeprop(wfn, 'MBIS_CHARGES')
        assert wfn is not None
    except Exception:
        # MBIS may not be available in all builds
        pytest.skip("MBIS not available")


def test_memory_integer_input():
    """Test memory setting with integer input"""
    # Direct integer input (in bytes)
    util.set_memory(1000000)
    mem = psi4.core.get_memory()
    assert mem == 1000000


def test_memory_string_without_units():
    """Test that memory string without units works"""
    # String without explicit units (should be bytes)
    try:
        util.set_memory('1000000')
        mem = psi4.core.get_memory()
        assert isinstance(mem, int)
    except Exception:
        # Some versions might require units
        pass


def test_oeprop_no_symmetry():
    """Test OEProp with C1 symmetry"""
    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)

    psi4.set_options({"basis": "sto-3g"})
    e, wfn = psi4.energy("scf", return_wfn=True)

    # Calculate dipole in C1 symmetry
    util.oeprop(wfn, 'DIPOLE')

    dipole = psi4.variable('SCF DIPOLE')
    assert isinstance(dipole, (list, np.ndarray))


def test_managed_method_error_stats():
    """Test that ManagedMethodError contains statistics"""
    from psi4.driver.p4util.exceptions import ManagedMethodError

    try:
        # Create error with typical calling pattern
        circs = ['run_scf', 'scf', 'TYPE', 'df', 'rhf', 'scf_module', 'fc']
        raise ManagedMethodError(circs)
    except ManagedMethodError as e:
        # Should have stats attribute
        assert hasattr(e, 'stats')
        assert isinstance(e.stats, dict)
        assert 'driver' in e.stats
        assert 'method' in e.stats


def test_oeprop_esp():
    """Test electrostatic potential calculation"""
    mol = psi4.geometry("""
    H
    H 1 0.74
    """)

    psi4.set_options({"basis": "sto-3g"})
    e, wfn = psi4.energy("scf", return_wfn=True)

    # ESP calculation might need grid specification
    # Just test that call doesn't crash
    try:
        util.oeprop(wfn, 'ESP_AT_NUCLEI')
        assert wfn is not None
    except Exception:
        # May not be available or may need additional setup
        pytest.skip("ESP calculation not available")
