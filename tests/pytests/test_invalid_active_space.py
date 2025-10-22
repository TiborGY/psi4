"""
Test that invalid ACTIVE space specifications are properly caught
and produce informative error messages instead of cryptic LAPACK failures.

This addresses the issue where mismatched ACTIVE orbital specifications
(e.g., requesting orbitals in irreps that don't exist or have insufficient orbitals)
would proceed to DETCI and fail with "DSYEV diagonalizer failed" errors.
"""

import pytest
import psi4


def test_invalid_active_exceeds_available_orbitals():
    """
    Test that ACTIVE specification exceeding available orbitals
    raises an informative InputException.

    This test uses H2 with cc-pVDZ and requests active orbitals in
    irreps that don't have enough (or any) orbitals, which was previously
    causing a DSYEV failure deep in the DETCI code.
    """
    psi4.core.clean()
    psi4.core.clean_options()

    mol = psi4.geometry("""
    H
    H 1 0.60
    """)

    psi4.set_options({
        'basis': 'cc-pVDZ',
        'docc': [1, 0, 0, 0, 0, 0, 0, 0],
        # This ACTIVE specification is invalid for H2/cc-pVDZ
        # as it requests orbitals in irreps that don't exist
        'active': [1, 0, 0, 0, 0, 1, 0, 0]
    })

    # Should raise InputException with informative error message
    with pytest.raises(psi4.InputException) as exc_info:
        psi4.energy('casscf')

    # Verify the error message is informative
    error_msg = str(exc_info.value)
    assert "irrep" in error_msg.lower(), "Error should mention irrep"
    assert ("orbital" in error_msg.lower() or
            "exceed" in error_msg.lower() or
            "insufficient" in error_msg.lower()), \
           "Error should mention orbitals/exceed/insufficient"


def test_invalid_active_insufficient_for_electrons():
    """
    Test that insufficient ACTIVE orbitals for occupied electrons
    raises an informative InputException.
    """
    psi4.core.clean()
    psi4.core.clean_options()

    mol = psi4.geometry("""
    O
    H 1 1.0
    H 1 1.0 2 104.5
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'restricted_docc': [1, 0, 0, 0],
        # ACTIVE space too small for the number of electrons
        'active': [1, 0, 0, 0]  # Only 1 orbital for 8 remaining electrons
    })

    # Should raise InputException with informative error message
    with pytest.raises(psi4.InputException) as exc_info:
        psi4.energy('casscf')

    # Verify the error message is informative
    error_msg = str(exc_info.value)
    assert ("insufficient" in error_msg.lower() or
            "electron" in error_msg.lower()), \
           "Error should mention insufficient orbitals for electrons"


def test_valid_active_space():
    """
    Test that a valid ACTIVE space specification works correctly.
    This is a control test to ensure our validation doesn't break valid inputs.
    """
    psi4.core.clean()
    psi4.core.clean_options()

    mol = psi4.geometry("""
    O
    H 1 1.0
    H 1 1.0 2 104.5
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'restricted_docc': [1, 0, 0, 0],
        'active': [3, 0, 1, 2],  # Valid active space
        'mcscf_maxiter': 2,  # Just a quick test, don't need convergence
        'mcscf_e_convergence': 1e-4,
        'mcscf_r_convergence': 1e-3
    })

    # Should not raise an exception
    try:
        energy = psi4.energy('casscf')
        # If it gets here, validation passed (even if CASSCF doesn't converge fully)
        assert True
    except psi4.InputException:
        # Should not get InputException for valid input
        pytest.fail("Valid ACTIVE space raised InputException")
