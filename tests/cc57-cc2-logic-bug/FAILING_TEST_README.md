# Tests That FAIL With Current Buggy Code

This directory contains tests that **CURRENTLY FAIL** due to the boolean logic bug in `psi4/src/psi4/cc/ccenergy/t2.cc:48`.

## Quick Start

Run the comprehensive test suite:

```bash
cd /home/user/psi4/tests/cc57-cc2-logic-bug
./run_all_tests.sh
```

**Current Result:** ❌ **FAILS** (exit code 1)

**After Fix:** ✓ **PASSES** (exit code 0)

---

## Test Files That FAIL Now (PASS After Fix)

### 1. `test_source_code.py` ⭐ **PRIMARY FAILING TEST**

**Status:** ❌ **CURRENTLY FAILS** (exit code 1)

**What it does:**
- Reads the actual source file `psi4/src/psi4/cc/ccenergy/t2.cc`
- Checks line 48 for the buggy pattern `|| ` vs fixed pattern `&&`
- Reports FAIL if bug exists, PASS if fixed

**Run it:**
```bash
python3 test_source_code.py
```

**Current output:**
```
❌ TEST FAILED - Bug is present in source code

The source code uses || (OR) operator, which makes the condition
ALWAYS TRUE, even for CC2/EOM_CC2.
```

**After fix output:**
```
✓ TEST PASSED - Bug has been fixed!

The source code correctly uses && (AND) operator.
```

---

### 2. `test_t2_logic.cc` (Unit Test)

**Status:** ⚠️ **DEMONSTRATES** the bug (not a pass/fail test)

**What it does:**
- Compiles and runs a C++ unit test
- Tests the boolean logic with || (buggy) vs && (fixed)
- Shows that buggy logic fails 2/6 tests

**Run it:**
```bash
g++ -std=c++11 -o test_t2_logic test_t2_logic.cc && ./test_t2_logic
```

**Output shows:**
```
BUGGY logic (||): 4/6 tests passed, 2 failed
FIXED logic (&&): 6/6 tests passed, 0 failed
```

This test always produces the same output (demonstrating the bug), but doesn't fail/pass based on the actual source code.

---

### 3. `run_all_tests.sh` **COMPREHENSIVE TEST SUITE**

**Status:** ❌ **CURRENTLY FAILS** (exit code 1)

**What it does:**
- Runs both test_t2_logic (demonstration) and test_source_code.py (integration)
- Exits with code 1 if bug exists, 0 if fixed

**Run it:**
```bash
./run_all_tests.sh
```

**Current exit code:** 1 (FAIL)
**After fix exit code:** 0 (PASS)

---

## How to Fix the Bug

**File:** `psi4/src/psi4/cc/ccenergy/t2.cc`
**Line:** 48

### Before (BUGGY):
```cpp
if (params_.wfn != "CC2" || params_.wfn != "EOM_CC2") {
```

### After (FIXED):
```cpp
if (params_.wfn != "CC2" && params_.wfn != "EOM_CC2") {
```

### To apply the fix:

1. Edit the file:
   ```bash
   vim psi4/src/psi4/cc/ccenergy/t2.cc
   ```

2. Go to line 48

3. Change `||` to `&&`

4. Save the file

5. Re-run the test suite:
   ```bash
   ./run_all_tests.sh
   ```

6. Verify it now **PASSES** (exit code 0)

---

## Test Verification Process

### Before Fix (Current State):

```bash
$ ./run_all_tests.sh
# ... output ...
❌ TEST FAILED - Bug still exists in source code
$ echo $?
1
```

### After Fix:

```bash
$ ./run_all_tests.sh
# ... output ...
✓ ALL TESTS PASSED - Bug has been fixed!
$ echo $?
0
```

---

## Why These Tests Are Important

1. **Verification:** These tests prove the bug exists in the current code
2. **Regression:** After fixing, they prevent the bug from being reintroduced
3. **Documentation:** They clearly show what the bug is and how to fix it
4. **CI/CD:** Can be integrated into continuous integration pipelines

---

## Demonstration Tests (Informational)

These tests demonstrate the bug but don't fail/pass based on source code:

### `demonstrate_bug.py`
Shows the boolean logic error with a truth table

```bash
python3 demonstrate_bug.py
```

**Output:**
```
Method       | Expected | Buggy OR  | Fixed AND  | Status
----------------------------------------------------------------------
CC2          | False    | True      | False      | ❌ BUG → ✓ FIXED
EOM_CC2      | False    | True      | False      | ❌ BUG → ✓ FIXED
CCSD         | True     | True      | True       | ✓ OK
```

---

## Integration with Test Framework

### For pytest/ctest:

The `test_source_code.py` script returns proper exit codes:
- 0 = Pass (bug fixed)
- 1 = Fail (bug exists)
- 2 = Inconclusive (file not found)

Can be integrated into CMake/CTest:

```cmake
add_test(NAME cc57_logic_bug_check
         COMMAND python3 test_source_code.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
```

---

## Summary Table

| File | Type | Current Status | After Fix |
|------|------|----------------|-----------|
| `test_source_code.py` | Integration | ❌ FAIL (exit 1) | ✓ PASS (exit 0) |
| `run_all_tests.sh` | Suite | ❌ FAIL (exit 1) | ✓ PASS (exit 0) |
| `test_t2_logic.cc` | Unit | ⚠️ Demo only | ⚠️ Demo only |
| `demonstrate_bug.py` | Demo | ℹ️ Info only | ℹ️ Info only |
| `input.dat` | Regression | ✓ PASS | ✓ PASS |

**Note:** `input.dat` passes both before and after because normal CC2 uses `cc2_t2_build()`, not the buggy `t2_build()`.

---

## Expected Test Results

### Test 1: Logic Unit Test (test_t2_logic.cc)
```
BUGGY logic (||): 4/6 tests passed, 2 failed  ❌
FIXED logic (&&): 6/6 tests passed, 0 failed  ✓
```

### Test 2: Source Code Check (test_source_code.py)
**Before Fix:**
```
❌ TEST FAILED - Bug is present in source code
Exit code: 1
```

**After Fix:**
```
✓ TEST PASSED - Bug has been fixed!
Exit code: 0
```

---

## For Developers

If you see test failures:

1. ✅ **Good!** The tests are working correctly
2. 🔧 Apply the fix to `t2.cc:48` (change `||` to `&&`)
3. ✓ Re-run tests to verify they now pass
4. 📝 Commit the fix with reference to this test

---

## Bug Reference

**Discovered:** 2025-11-09
**Location:** `psi4/src/psi4/cc/ccenergy/t2.cc:48`
**Type:** Boolean logic error (OR instead of AND)
**Severity:** Medium (affects edge cases, not normal CC2 calculations)
**Fix:** One character change (`||` → `&&`)
**Test:** This directory (`tests/cc57-cc2-logic-bug/`)
