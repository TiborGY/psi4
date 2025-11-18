# DIIS POC Testing Guide

This document describes how to test the libdiis-based DIIS proof of concept implementation for the ccenergy module.

## Overview

The POC implementation replaces ~250 lines of custom DIIS code in `diis_RHF.cc` with ~60 lines that leverage the centralized `libdiis/DIISManager` infrastructure. This testing phase validates that the new implementation:

1. **Functions correctly** - Produces correct energies
2. **Converges identically** - Same iteration count and convergence behavior
3. **Performs adequately** - No significant performance regression
4. **Integrates cleanly** - Works with existing DPD infrastructure

## Test Scripts

Four test scripts are provided:

### 1. `test_diis_poc.py` - Basic Functional Test
**Purpose**: Verify POC implementation works and produces correct energies

**Usage**:
```bash
python3 test_diis_poc.py [--verbose]
```

**What it does**:
- Runs RHF-CCSD/6-31G** calculation on H2O
- Checks for POC-specific output messages
- Compares final energy against reference
- Reports PASS/FAIL

**Success criteria**:
- Calculation completes without errors
- Energy matches reference within 1e-9 Hartree
- DIIS messages appear in output

---

### 2. `compare_diis_implementations.py` - Implementation Comparison
**Purpose**: Compare original vs POC implementation side-by-side

**Usage**:
```bash
python3 compare_diis_implementations.py [--molecule MOLECULE] [--basis BASIS]
```

**What it does**:
- Framework for comparing two builds (original vs POC)
- Tracks energy accuracy, iteration count, convergence rate
- Currently a template - requires two separate builds

**Note**: Full comparison requires:
1. Build Psi4 WITHOUT `-DUSE_LIBDIIS_POC` (original)
2. Build Psi4 WITH `-DUSE_LIBDIIS_POC` (POC)
3. Run same test with both builds
4. Compare outputs

---

### 3. `analyze_diis_convergence.py` - Convergence Analysis
**Purpose**: Parse and analyze CCSD convergence behavior

**Usage**:
```bash
python3 analyze_diis_convergence.py output.dat [--plot]
```

**What it does**:
- Extracts iteration-by-iteration data from output
- Prints convergence table
- Optionally plots convergence graph (requires matplotlib)

**Example output**:
```
CCSD CONVERGENCE ANALYSIS
================================================================================
 Iter             Energy              ΔE                    DIIS Info
--------------------------------------------------------------------------------
    1  -76.2289123456  0.00000000e+00     stored vector (subspace = 1)
    2  -76.2301234567 -1.21111111e-03     extrapolated with 2 vectors
    3  -76.2311567890 -1.03333333e-03     extrapolated with 3 vectors
   ...
```

---

### 4. `run_diis_poc_tests.sh` - Master Test Runner
**Purpose**: Orchestrate complete testing workflow

**Usage**:
```bash
./run_diis_poc_tests.sh [--quick|--full]
```

**Modes**:
- **Default**: Runs recommended tests (cc1, cc2, cc4a, cc13, cc54)
- **--quick**: Runs minimal test (cc1 only)
- **--full**: Runs all CC tests

**Workflow**:
1. Check build configuration
2. Run basic functional test
3. Run regression tests
4. Analyze results

---

## Testing Workflow

### Step 1: Build Psi4 with POC Flag

```bash
# From Psi4 root directory
mkdir objdir_poc
cd objdir_poc

# Configure with POC flag enabled
cmake .. \
    -DCMAKE_CXX_FLAGS="-DUSE_LIBDIIS_POC" \
    -DCMAKE_BUILD_TYPE=Release

# Build
make -j$(nproc)

# Install (optional)
make install
```

### Step 2: Run Basic Tests

```bash
# Quick validation
./run_diis_poc_tests.sh --quick
```

Expected output:
```
========================================================================
                    DIIS POC VALIDATION TEST SUITE
========================================================================

[1/4] Checking build configuration...
✓ POC code is present in source
✓ Build includes USE_LIBDIIS_POC flag

[2/4] Running basic functional test...
✓ Basic functional test PASSED

[3/4] Running regression tests...
✓ cc1 test PASSED

[4/4] Analyzing results...
Convergence analysis complete
```

### Step 3: Run Full Regression Suite

```bash
# Standard tests
./run_diis_poc_tests.sh

# Or run comprehensive suite
./run_diis_poc_tests.sh --full
```

### Step 4: Compare Implementations (Optional)

To compare original vs POC side-by-side:

```bash
# Build 1: Original implementation
mkdir objdir_original
cd objdir_original
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
cd ..

# Run test and save output
objdir_original/stage/bin/psi4 tests/cc1/input.dat -o cc1_original.out

# Build 2: POC implementation
mkdir objdir_poc
cd objdir_poc
cmake .. -DCMAKE_CXX_FLAGS="-DUSE_LIBDIIS_POC" -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
cd ..

# Run test and save output
objdir_poc/stage/bin/psi4 tests/cc1/input.dat -o cc1_poc.out

# Compare
diff cc1_original.out cc1_poc.out
python3 analyze_diis_convergence.py cc1_original.out
python3 analyze_diis_convergence.py cc1_poc.out
```

### Step 5: Analyze Results

Check for:

1. **Energy accuracy**: Final energies should match within 1e-9 Hartree
2. **Iteration count**: Should be identical
3. **Convergence pattern**: Energy changes per iteration should match
4. **DIIS messages**: POC should show "DIIS: extrapolated with N vectors"

---

## Success Criteria

The POC is considered successful if:

### ✓ Functional Requirements
- [ ] All regression tests pass (cc1, cc2, cc4a, cc13, cc54)
- [ ] Energies match original within 1e-9 Hartree
- [ ] Iteration counts are identical
- [ ] Convergence behavior is identical

### ✓ Technical Requirements
- [ ] DIIS extrapolation occurs (subspace_size >= 2)
- [ ] DPD buffers are handled correctly
- [ ] No memory leaks or errors
- [ ] POC-specific messages appear in output

### ✓ Performance Requirements
- [ ] Runtime within 5% of original implementation
- [ ] Memory usage comparable

---

## Interpreting Results

### Output Messages to Look For

**POC Initialization** (in output):
```
POC: Using libdiis for DIIS extrapolation
```

**DIIS Extrapolation** (in output):
```
DIIS: extrapolated with N vectors
```

**Normal Termination**:
```
*** Psi4 exiting successfully.
```

### Common Issues

**Issue**: No POC messages in output
- **Cause**: Build doesn't include `-DUSE_LIBDIIS_POC` flag
- **Fix**: Rebuild with flag enabled

**Issue**: Energy differs from reference
- **Cause**: DIIS implementation bug
- **Action**: Check error vector computation (R1, R2)

**Issue**: Crashes or segfaults
- **Cause**: DPD buffer handling issue
- **Action**: Check file2/buf4 initialization and cleanup

**Issue**: Iteration count differs
- **Cause**: DIIS behavior difference
- **Action**: Check DIIS parameters (maxvec, mindiis, removal policy)

---

## Recommended Test Suite

**Minimal** (2-3 minutes):
- cc1 (RHF-CCSD H2O)

**Standard** (10-15 minutes):
- cc1 (RHF-CCSD H2O)
- cc2 (RHF-CCSD H2O+)
- cc4a (RHF-CCSD CN radical)
- cc13 (RHF-CCSD+ CN)
- cc54 (RHF-CCSD(T) H2O)

**Comprehensive** (30+ minutes):
- All tests matching `^cc` pattern
- Covers RHF, ROHF, UHF (POC only affects RHF)

---

## Next Steps After Testing

If all tests pass:

1. **Document Results**: Create `POC_RESULTS.md` with:
   - Test outputs
   - Energy comparisons
   - Performance metrics
   - Any observations

2. **Code Review**: Prepare POC for review:
   - Ensure code is clean and well-commented
   - Verify all preprocessor guards are correct
   - Check for any edge cases

3. **Decision Point**: Evaluate whether to:
   - Proceed to Phase 2 (ROHF/UHF implementation)
   - Request feedback from core developers
   - Make adjustments based on findings

4. **Long-term Plan**: If POC validates approach:
   - Extend to ROHF and UHF
   - Remove preprocessor guards
   - Delete old implementation
   - Repeat for cclambda and ccresponse modules

---

## Additional Resources

- **POC Implementation Plan**: `CCENERGY_DIIS_POC_PLAN.md`
- **CC DIIS Analysis**: `CC_DIIS_CONSOLIDATION_ANALYSIS.md`
- **libdiis Documentation**: `psi4/src/psi4/libdiis/diismanager.h`
- **Original Implementation**: `psi4/src/psi4/cc/ccenergy/diis_RHF.cc`
- **POC Implementation**: `psi4/src/psi4/cc/ccenergy/diis_RHF_libdiis.cc`

---

## Contact

For questions or issues with the POC:
- Review implementation in `diis_RHF_libdiis.cc`
- Check `CC_DIIS_CONSOLIDATION_ANALYSIS.md` for technical details
- Consult libdiis documentation and examples in occ/dfocc modules

---

**Last Updated**: 2025-11-18
**Author**: Claude (Anthropic AI)
**Status**: Testing Phase (Phase 3 of 5)
