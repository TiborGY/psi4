"""
Test to verify that psi4.set_output_file() does not leak file handlers.

This addresses the issue where calling set_output_file() in a loop accumulates
logging FileHandlers, eventually hitting the OS file limit and causing a segfault.
"""
import pytest
import tempfile
import os
import logging

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def test_output_file_handler_cleanup():
    """Test that set_output_file() removes old logging FileHandlers."""

    # Get reference to psi4 logger
    logger = psi4.logger

    # Record initial handler counts
    initial_total = len(logger.handlers)
    initial_file_handlers = sum(1 for h in logger.handlers if isinstance(h, logging.FileHandler))

    with tempfile.TemporaryDirectory() as tmpdir:
        # Call set_output_file multiple times
        n_iterations = 20
        for i in range(n_iterations):
            outfile = os.path.join(tmpdir, f"test_{i}.dat")
            psi4.set_output_file(outfile)

            # Count handlers
            total_handlers = len(logger.handlers)
            file_handlers = sum(1 for h in logger.handlers if isinstance(h, logging.FileHandler))

            # Should have at most 1 more FileHandler than initially
            # (the handler for the current output file)
            assert file_handlers <= initial_file_handlers + 1, \
                f"Iteration {i}: Found {file_handlers} FileHandlers, expected at most {initial_file_handlers + 1}"

        # Final check - handler count should not have grown significantly
        final_total = len(logger.handlers)
        final_file_handlers = sum(1 for h in logger.handlers if isinstance(h, logging.FileHandler))

        # Should have at most 1 more FileHandler than initially
        assert final_file_handlers <= initial_file_handlers + 1, \
            f"Final: {final_file_handlers} FileHandlers, expected at most {initial_file_handlers + 1}"

        # Total handler count should not have grown by more than 1
        assert final_total <= initial_total + 1, \
            f"Handler count grew from {initial_total} to {final_total}, expected growth of at most 1"


def test_output_file_handler_cleanup_with_computation():
    """Test handler cleanup during actual computations in a loop."""

    # Get reference to psi4 logger
    logger = psi4.logger
    initial_file_handlers = sum(1 for h in logger.handlers if isinstance(h, logging.FileHandler))

    # Simple water molecule
    mol = psi4.geometry("""
        O
        H 1 0.96
        H 1 0.96 2 104.5
        symmetry c1
    """)

    psi4.set_options({
        'basis': 'sto-3g',
        'scf_type': 'pk',
        'e_convergence': 1e-8,
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        # Run multiple iterations with actual calculations
        n_iterations = 10
        for i in range(n_iterations):
            outfile = os.path.join(tmpdir, f"output_{i:03d}.dat")

            # Set output file and run calculation
            psi4.set_output_file(outfile)
            psi4.energy('scf')
            psi4.core.close_outfile()

            # Check that we haven't leaked handlers
            file_handlers = sum(1 for h in logger.handlers if isinstance(h, logging.FileHandler))
            assert file_handlers <= initial_file_handlers + 1, \
                f"Iteration {i}: Leaked handlers, found {file_handlers}"

        # Final check
        final_file_handlers = sum(1 for h in logger.handlers if isinstance(h, logging.FileHandler))
        assert final_file_handlers <= initial_file_handlers + 1, \
            f"After loop: {final_file_handlers} FileHandlers, expected at most {initial_file_handlers + 1}"


def test_output_file_switching():
    """Test switching between different output files."""

    logger = psi4.logger
    initial_file_handlers = sum(1 for h in logger.handlers if isinstance(h, logging.FileHandler))

    with tempfile.TemporaryDirectory() as tmpdir:
        file_a = os.path.join(tmpdir, "output_a.dat")
        file_b = os.path.join(tmpdir, "output_b.dat")

        # Switch back and forth between two files
        for i in range(10):
            if i % 2 == 0:
                psi4.set_output_file(file_a)
            else:
                psi4.set_output_file(file_b)

            file_handlers = sum(1 for h in logger.handlers if isinstance(h, logging.FileHandler))
            assert file_handlers <= initial_file_handlers + 1, \
                f"Iteration {i}: Found {file_handlers} FileHandlers after switching"

        # Cleanup
        psi4.core.close_outfile()


def test_output_file_append_mode():
    """Test that handler cleanup works in append mode."""

    logger = psi4.logger
    initial_file_handlers = sum(1 for h in logger.handlers if isinstance(h, logging.FileHandler))

    with tempfile.TemporaryDirectory() as tmpdir:
        outfile = os.path.join(tmpdir, "output.dat")

        # First write
        psi4.set_output_file(outfile, append=False)
        file_handlers = sum(1 for h in logger.handlers if isinstance(h, logging.FileHandler))
        assert file_handlers <= initial_file_handlers + 1

        # Multiple appends
        for i in range(5):
            psi4.set_output_file(outfile, append=True)
            file_handlers = sum(1 for h in logger.handlers if isinstance(h, logging.FileHandler))
            assert file_handlers <= initial_file_handlers + 1, \
                f"Append {i}: Found {file_handlers} FileHandlers"

        psi4.core.close_outfile()
