#!/usr/bin/env python
"""
Simple test to verify that logging handlers are properly cleaned up
when psi4.set_output_file() is called multiple times.
"""
import psi4
import tempfile
import os

def test_handler_cleanup():
    """Test that set_output_file() removes old handlers."""

    print("Testing logging handler cleanup...")

    # Get reference to psi4 logger
    logger = psi4.logger

    initial_handler_count = len(logger.handlers)
    print(f"Initial handler count: {initial_handler_count}")

    with tempfile.TemporaryDirectory() as tmpdir:

        # Call set_output_file multiple times
        for i in range(10):
            outfile = os.path.join(tmpdir, f"test_{i}.dat")
            psi4.set_output_file(outfile)

            # Check handler count
            handler_count = len(logger.handlers)
            file_handler_count = sum(1 for h in logger.handlers if isinstance(h, __import__('logging').FileHandler))

            print(f"Iteration {i+1:2d}: Total handlers = {handler_count}, FileHandlers = {file_handler_count}")

            # Should have at most 1 FileHandler at any time
            # (may have other handler types from initial setup)
            if file_handler_count > 1:
                print(f"  ERROR: Found {file_handler_count} FileHandlers, expected at most 1!")
                return False

    final_handler_count = len(logger.handlers)
    print(f"\nFinal handler count: {final_handler_count}")

    # The count shouldn't have grown by more than 1 (the current file handler)
    if final_handler_count > initial_handler_count + 1:
        print(f"FAIL: Handler count grew from {initial_handler_count} to {final_handler_count}")
        return False
    else:
        print(f"PASS: Handler count stable (grew by at most 1)")
        return True

if __name__ == "__main__":
    try:
        # Initialize psi4
        psi4.core.initialize()

        success = test_handler_cleanup()

        if success:
            print("\n✓ Test PASSED: Handlers properly cleaned up")
        else:
            print("\n✗ Test FAILED: Handler leak detected")
            exit(1)

    except Exception as e:
        print(f"\n✗ Test FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
    finally:
        psi4.core.finalize()
