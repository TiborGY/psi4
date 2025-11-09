#!/usr/bin/env python3
"""
Integration test that checks the actual source code for the bug.

This test FAILS if the buggy code is present and PASSES if fixed.

Usage:
    python3 test_source_code.py

Exit codes:
    0 - Test PASSED (bug is fixed, using &&)
    1 - Test FAILED (bug exists, using ||)
    2 - Could not determine (file not found or line changed)
"""

import sys
import os
import re

# Path to the buggy file (relative to this test directory)
SOURCE_FILE = "../../psi4/src/psi4/cc/ccenergy/t2.cc"
BUGGY_LINE_NUM = 48

def check_source_code():
    """
    Check if the bug exists in the source code.

    Returns:
        'BUGGY' if using || (bug exists)
        'FIXED' if using && (bug fixed)
        'UNKNOWN' if cannot determine
    """
    # Get the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    source_path = os.path.join(script_dir, SOURCE_FILE)

    if not os.path.exists(source_path):
        print(f"ERROR: Source file not found at {source_path}")
        return 'UNKNOWN'

    try:
        with open(source_path, 'r') as f:
            lines = f.readlines()

        if len(lines) < BUGGY_LINE_NUM:
            print(f"ERROR: File has fewer than {BUGGY_LINE_NUM} lines")
            return 'UNKNOWN'

        # Get the target line (0-indexed, so line 48 is index 47)
        target_line = lines[BUGGY_LINE_NUM - 1]

        print(f"Checking line {BUGGY_LINE_NUM}:")
        print(f"  {target_line.rstrip()}")
        print()

        # Check for buggy pattern: != "CC2" || != "EOM_CC2"
        # or fixed pattern: != "CC2" && != "EOM_CC2"

        # Look for the condition - it should have both "CC2" and "EOM_CC2"
        if '"CC2"' not in target_line or '"EOM_CC2"' not in target_line:
            print("WARNING: Line doesn't contain both 'CC2' and 'EOM_CC2'")
            print("The code may have been refactored.")
            return 'UNKNOWN'

        # Extract the condition part
        # Look for pattern: wfn != "CC2" [||/&&] wfn != "EOM_CC2"
        or_pattern = re.search(r'!= *"CC2" *\|\| *.* *= *"EOM_CC2"', target_line)
        and_pattern = re.search(r'!= *"CC2" *&& *.* *= *"EOM_CC2"', target_line)

        if or_pattern:
            print("❌ BUGGY pattern detected: Using || (OR) operator")
            print(f"   Pattern: {or_pattern.group()}")
            return 'BUGGY'
        elif and_pattern:
            print("✓ FIXED pattern detected: Using && (AND) operator")
            print(f"   Pattern: {and_pattern.group()}")
            return 'FIXED'
        else:
            print("WARNING: Could not detect || or && pattern")
            print("The code may have been refactored.")
            return 'UNKNOWN'

    except Exception as e:
        print(f"ERROR: {e}")
        return 'UNKNOWN'

def main():
    print("=" * 70)
    print("Integration Test: t2.cc:48 Boolean Logic Bug")
    print("=" * 70)
    print()
    print("This test checks the actual source code for the bug.")
    print()

    result = check_source_code()

    print()
    print("=" * 70)
    print("TEST RESULT:")
    print("=" * 70)

    if result == 'BUGGY':
        print()
        print("❌ TEST FAILED - Bug is present in source code")
        print()
        print("The source code uses || (OR) operator, which makes the condition")
        print("ALWAYS TRUE, even for CC2/EOM_CC2.")
        print()
        print("Fix required:")
        print(f"  File: {SOURCE_FILE}")
        print(f"  Line: {BUGGY_LINE_NUM}")
        print("  Change: || to &&")
        print()
        print("Before: if (params_.wfn != \"CC2\" || params_.wfn != \"EOM_CC2\")")
        print("After:  if (params_.wfn != \"CC2\" && params_.wfn != \"EOM_CC2\")")
        print()
        return 1

    elif result == 'FIXED':
        print()
        print("✓ TEST PASSED - Bug has been fixed!")
        print()
        print("The source code correctly uses && (AND) operator.")
        print("CC2 and EOM_CC2 will properly skip CCSD-specific terms.")
        print()
        return 0

    else:  # UNKNOWN
        print()
        print("⚠ TEST INCONCLUSIVE - Could not determine bug status")
        print()
        print("The source code may have been refactored or the file was not found.")
        print("Manual verification required.")
        print()
        return 2

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
