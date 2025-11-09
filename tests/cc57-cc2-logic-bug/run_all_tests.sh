#!/bin/bash
#
# Comprehensive test suite for the t2.cc:48 boolean logic bug
#
# This script runs multiple tests:
# 1. Unit test (C++ logic test) - always shows the bug
# 2. Source code test (Python) - FAILS if bug exists, PASSES if fixed
# 3. Demonstration script - shows the logic error
#
# Usage: ./run_all_tests.sh

set -e  # Exit on error

echo "======================================================================"
echo "Comprehensive Test Suite: t2.cc:48 Boolean Logic Bug"
echo "======================================================================"
echo ""

# Test 1: C++ Unit Test
echo "======================================================================"
echo "TEST 1: C++ Unit Test (Logic Demonstration)"
echo "======================================================================"
echo "This test compares buggy || logic vs fixed && logic"
echo ""

if [ ! -f test_t2_logic ]; then
    echo "Compiling test_t2_logic.cc..."
    g++ -std=c++11 -o test_t2_logic test_t2_logic.cc
    echo ""
fi

./test_t2_logic
echo ""

# Test 2: Source Code Check
echo "======================================================================"
echo "TEST 2: Source Code Integration Test"
echo "======================================================================"
echo "This test checks the actual source code and will:"
echo "  - FAIL (exit 1) if the bug exists in the source"
echo "  - PASS (exit 0) if the bug has been fixed"
echo ""

python3 test_source_code.py
TEST_RESULT=$?

echo ""
echo "======================================================================"

if [ $TEST_RESULT -eq 0 ]; then
    echo "✓ ALL TESTS PASSED - Bug has been fixed!"
    echo "======================================================================"
    exit 0
elif [ $TEST_RESULT -eq 1 ]; then
    echo "❌ TEST FAILED - Bug still exists in source code"
    echo "======================================================================"
    echo ""
    echo "To fix the bug:"
    echo "  1. Edit psi4/src/psi4/cc/ccenergy/t2.cc line 48"
    echo "  2. Change: if (params_.wfn != \"CC2\" || params_.wfn != \"EOM_CC2\")"
    echo "  3. To:     if (params_.wfn != \"CC2\" && params_.wfn != \"EOM_CC2\")"
    echo "  4. Re-run this test suite to verify the fix"
    echo ""
    exit 1
else
    echo "⚠ TEST INCONCLUSIVE - Could not determine bug status"
    echo "======================================================================"
    exit 2
fi
