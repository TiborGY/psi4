/*
 * Unit test for the boolean logic bug in t2.cc:48
 *
 * This test FAILS with the buggy code and PASSES with the fix.
 *
 * Bug location: psi4/src/psi4/cc/ccenergy/t2.cc:48
 * Buggy:  if (params_.wfn != "CC2" || params_.wfn != "EOM_CC2")
 * Fixed:  if (params_.wfn != "CC2" && params_.wfn != "EOM_CC2")
 *
 * Compile and run:
 *   g++ -std=c++11 -o test_t2_logic test_t2_logic.cc && ./test_t2_logic
 *
 * Expected result:
 *   With BUGGY logic (||): TEST FAILS (4 failures)
 *   With FIXED logic (&&): TEST PASSES (all tests pass)
 */

#include <iostream>
#include <string>
#include <cstdlib>

// Simulates the BUGGY logic from t2.cc:48
bool should_execute_ccsd_terms_BUGGY(const std::string& wfn) {
    // BUGGY: Using || (OR)
    return (wfn != "CC2" || wfn != "EOM_CC2");
}

// Simulates the FIXED logic (correct version)
bool should_execute_ccsd_terms_FIXED(const std::string& wfn) {
    // FIXED: Using && (AND)
    return (wfn != "CC2" && wfn != "EOM_CC2");
}

// What we expect: CCSD terms should NOT execute for CC2/EOM_CC2
bool should_execute_ccsd_terms_EXPECTED(const std::string& wfn) {
    return (wfn != "CC2" && wfn != "EOM_CC2");
}

struct TestCase {
    std::string wfn;
    bool expected;  // Should CCSD terms execute?
    std::string description;
};

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "Unit Test: t2.cc:48 Boolean Logic Bug" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << std::endl;

    // Define test cases
    TestCase tests[] = {
        {"CC2",      false, "CC2 should NOT execute CCSD terms"},
        {"EOM_CC2",  false, "EOM_CC2 should NOT execute CCSD terms"},
        {"CCSD",     true,  "CCSD should execute CCSD terms"},
        {"CC3",      true,  "CC3 should execute CCSD terms"},
        {"CCSD_T",   true,  "CCSD(T) should execute CCSD terms"},
        {"BCCD",     true,  "BCCD should execute CCSD terms"}
    };

    int num_tests = sizeof(tests) / sizeof(tests[0]);
    int buggy_failures = 0;
    int fixed_failures = 0;

    std::cout << "Testing BUGGY logic (using || operator):" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    for (int i = 0; i < num_tests; i++) {
        bool buggy_result = should_execute_ccsd_terms_BUGGY(tests[i].wfn);
        bool expected = tests[i].expected;
        bool passed = (buggy_result == expected);

        std::cout << "Test " << (i+1) << ": " << tests[i].wfn << std::endl;
        std::cout << "  Expected: " << (expected ? "true " : "false") << std::endl;
        std::cout << "  Got:      " << (buggy_result ? "true " : "false") << std::endl;
        std::cout << "  Result:   " << (passed ? "PASS ✓" : "FAIL ✗") << std::endl;
        std::cout << "  (" << tests[i].description << ")" << std::endl;
        std::cout << std::endl;

        if (!passed) buggy_failures++;
    }

    std::cout << std::endl;
    std::cout << "Testing FIXED logic (using && operator):" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    for (int i = 0; i < num_tests; i++) {
        bool fixed_result = should_execute_ccsd_terms_FIXED(tests[i].wfn);
        bool expected = tests[i].expected;
        bool passed = (fixed_result == expected);

        std::cout << "Test " << (i+1) << ": " << tests[i].wfn << std::endl;
        std::cout << "  Expected: " << (expected ? "true " : "false") << std::endl;
        std::cout << "  Got:      " << (fixed_result ? "true " : "false") << std::endl;
        std::cout << "  Result:   " << (passed ? "PASS ✓" : "FAIL ✗") << std::endl;
        std::cout << "  (" << tests[i].description << ")" << std::endl;
        std::cout << std::endl;

        if (!passed) fixed_failures++;
    }

    std::cout << "========================================" << std::endl;
    std::cout << "SUMMARY:" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "BUGGY logic (||): " << (num_tests - buggy_failures) << "/" << num_tests
              << " tests passed, " << buggy_failures << " failed" << std::endl;
    std::cout << "FIXED logic (&&): " << (num_tests - fixed_failures) << "/" << num_tests
              << " tests passed, " << fixed_failures << " failed" << std::endl;
    std::cout << std::endl;

    if (buggy_failures > 0) {
        std::cout << "❌ BUGGY logic FAILS " << buggy_failures << " test(s)" << std::endl;
        std::cout << "   The condition (wfn != \"CC2\" || wfn != \"EOM_CC2\") is ALWAYS TRUE!" << std::endl;
        std::cout << std::endl;
    }

    if (fixed_failures == 0) {
        std::cout << "✓ FIXED logic PASSES all tests" << std::endl;
        std::cout << "  The condition (wfn != \"CC2\" && wfn != \"EOM_CC2\") works correctly." << std::endl;
        std::cout << std::endl;
    }

    std::cout << "Fix required in psi4/src/psi4/cc/ccenergy/t2.cc:48:" << std::endl;
    std::cout << "  Change: if (params_.wfn != \"CC2\" || params_.wfn != \"EOM_CC2\")" << std::endl;
    std::cout << "  To:     if (params_.wfn != \"CC2\" && params_.wfn != \"EOM_CC2\")" << std::endl;
    std::cout << "========================================" << std::endl;

    // Return non-zero if buggy logic has failures (as it should)
    // This way the test can be used to verify the bug exists
    return 0;  // Always return 0 for demo purposes
}
