# Quick Start Guide: Tests That FAIL Now

## ⭐ Run the Failing Test

```bash
cd /home/user/psi4/tests/cc57-cc2-logic-bug

# Option 1: Run comprehensive test suite (RECOMMENDED)
./run_all_tests.sh
# Current result: ❌ FAILS (exit code 1)
# After fix: ✓ PASSES (exit code 0)

# Option 2: Run just the source code check
python3 test_source_code.py
# Current result: ❌ FAILS (exit code 1)
# After fix: ✓ PASSES (exit code 0)

# Option 3: Run the C++ logic demonstration
g++ -std=c++11 -o test_t2_logic test_t2_logic.cc && ./test_t2_logic
# Shows: BUGGY logic fails 2/6 tests, FIXED logic passes 6/6
```

---

## What's in This Directory

### Tests That FAIL With Current Code ⭐

| File | What It Does | Current Status | After Fix |
|------|--------------|----------------|-----------|
| `test_source_code.py` | Checks actual source for bug | ❌ FAIL | ✓ PASS |
| `run_all_tests.sh` | Runs all tests | ❌ FAIL | ✓ PASS |
| `test_t2_logic.cc` | C++ unit test | Shows bug | Shows bug |

### Demonstrations (Informational)

| File | What It Does |
|------|--------------|
| `demonstrate_bug.py` | Truth table showing logic error |
| `input.dat` | Psi4 regression test (always passes) |

### Documentation

| File | What It Contains |
|------|------------------|
| `FAILING_TEST_README.md` | Detailed info about failing tests |
| `QUICK_START.md` | This file |
| `SUMMARY.md` | Comprehensive overview |
| `README.md` | Technical documentation |

---

## The Bug (One Line Summary)

**File:** `psi4/src/psi4/cc/ccenergy/t2.cc:48`

```cpp
// BUGGY (current):
if (params_.wfn != "CC2" || params_.wfn != "EOM_CC2")

// FIXED (correct):
if (params_.wfn != "CC2" && params_.wfn != "EOM_CC2")
```

**Change required:** `||` → `&&` (one character!)

---

## Quick Verification

### Before Fix:
```bash
$ python3 test_source_code.py
❌ TEST FAILED - Bug is present in source code
$ echo $?
1
```

### After Fix:
```bash
$ python3 test_source_code.py
✓ TEST PASSED - Bug has been fixed!
$ echo $?
0
```

---

## Why This Matters

The buggy OR logic is **always true**, even for CC2/EOM_CC2:
- `"CC2" != "EOM_CC2"` is true → condition is true ❌
- `"EOM_CC2" != "CC2"` is true → condition is true ❌

The correct AND logic properly excludes CC2/EOM_CC2:
- Both must be true for condition to be true ✓
- For CC2, first part is false → condition is false ✓

---

## Next Steps

1. **Verify the bug:** `./run_all_tests.sh` (should FAIL)
2. **Fix the bug:** Change `||` to `&&` in `t2.cc:48`
3. **Verify the fix:** `./run_all_tests.sh` (should PASS)
4. **Commit:** Include reference to this test case

See `FAILING_TEST_README.md` for complete documentation.
