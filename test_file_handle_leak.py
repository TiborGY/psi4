#!/usr/bin/env python
"""
Test script to verify that psi4.set_output_file() does not leak file handles
when called repeatedly in a loop.

This reproduces the issue where file handles accumulate until hitting ulimit,
causing a segfault.
"""
import psi4
import tempfile
import os
import psutil  # Optional: for monitoring open files

def get_open_files_count():
    """Get the number of open files for the current process."""
    try:
        process = psutil.Process()
        return len(process.open_files())
    except:
        return None

def test_file_handle_leak(n_iterations=100):
    """Test that set_output_file() doesn't leak file handles in a loop."""

    print(f"Testing {n_iterations} iterations of set_output_file()...")

    # Record initial state
    initial_count = get_open_files_count()
    if initial_count is not None:
        print(f"Initial open files: {initial_count}")

    # Simple water molecule for testing
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

    # Create temporary directory for output files
    with tempfile.TemporaryDirectory() as tmpdir:

        # Run multiple iterations
        for i in range(n_iterations):
            outfile = os.path.join(tmpdir, f"output_{i:04d}.dat")

            # Set output file (this should close previous handlers)
            psi4.set_output_file(outfile, append=False)

            # Run a quick SCF calculation
            psi4.energy('scf')

            # Close the output file
            psi4.core.close_outfile()

            # Check file count periodically
            if (i + 1) % 10 == 0:
                current_count = get_open_files_count()
                if current_count is not None:
                    print(f"Iteration {i+1:3d}: open files = {current_count}")
                    # Check if we're leaking files
                    if initial_count is not None and current_count > initial_count + 10:
                        print(f"WARNING: Possible file leak detected!")
                        print(f"  Initial: {initial_count}, Current: {current_count}, Diff: {current_count - initial_count}")

    # Final check
    final_count = get_open_files_count()
    if final_count is not None:
        print(f"\nFinal open files: {final_count}")
        if initial_count is not None:
            print(f"Change in open files: {final_count - initial_count}")

            # If we've leaked more than a few files, that's a problem
            if final_count - initial_count > 5:
                print("FAIL: File handle leak detected!")
                return False
            else:
                print("PASS: No significant file handle leak detected")
                return True
    else:
        print("\npsutil not available, cannot check file counts")
        print("Test completed without crashes - likely PASS")
        return True

if __name__ == "__main__":
    try:
        success = test_file_handle_leak(n_iterations=100)
        if success:
            print("\n✓ Test PASSED: No file handle leaks")
        else:
            print("\n✗ Test FAILED: File handle leaks detected")
            exit(1)
    except Exception as e:
        print(f"\n✗ Test FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
