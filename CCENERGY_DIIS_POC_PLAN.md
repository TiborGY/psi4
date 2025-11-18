# ccenergy RHF DIIS Proof of Concept - Implementation Plan

## Overview

This document provides a detailed, step-by-step plan for implementing a proof of concept migration of ccenergy's RHF DIIS implementation from custom code to libdiis/DIISManager.

**Goal:** Validate that libdiis can successfully replace custom DIIS in cc modules
**Target:** ccenergy RHF only (smallest, lowest-risk scope)
**Success Metric:** Identical energies and convergence with reduced code

---

## Prerequisites

### Required Knowledge
- [ ] Familiarity with Psi4's DPD (Distributed Packed Data) library
- [ ] Understanding of CCSD amplitude iterations (T1, T2)
- [ ] Basic knowledge of DIIS algorithm
- [ ] C++ and Psi4 build system experience

### Required Access
- [ ] Psi4 development environment set up
- [ ] Ability to compile Psi4 from source
- [ ] Access to regression test suite
- [ ] Benchmarking tools available

### Recommended Reading
- Current implementation: `psi4/src/psi4/cc/ccenergy/diis_RHF.cc`
- libdiis interface: `psi4/src/psi4/libdiis/diismanager.h`
- libdiis Python backend: `psi4/driver/procrouting/diis.py`
- CC DIIS analysis: `CC_DIIS_CONSOLIDATION_ANALYSIS.md`

---

## Phase 1: Analysis & Setup (2 hours)

### Step 1.1: Understand Current Implementation

**Task:** Deep dive into existing RHF DIIS code

```bash
# Read the current implementation
cat psi4/src/psi4/cc/ccenergy/diis_RHF.cc
```

**Document:**
- [ ] Vector length calculation method
- [ ] Error vector construction (T1_new - T1_old, T2_new - T2_old)
- [ ] Amplitude vector storage
- [ ] B matrix construction approach
- [ ] Conditioning method (diis_invert_B)
- [ ] Extrapolation procedure
- [ ] PSIO file usage (PSIF_CC_DIIS_ERR, PSIF_CC_DIIS_AMP)

**Key Code Sections to Understand:**

1. **Vector length calculation** (lines 72-80):
```cpp
global_dpd_->file2_init(&T1, PSIF_CC_TAMPS, 0, 0, 1, "tIA");
global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
auto vector_length = 0;
for (int h = 0; h < nirreps; h++) {
    vector_length += T1.params->rowtot[h] * T1.params->coltot[h];
    vector_length += T2.params->rowtot[h] * T2.params->coltot[h];
}
```

2. **Error vector construction** (lines 88-124):
```cpp
// Flattens T1_new - T1_old into 1D array
for (int h = 0; h < nirreps; h++)
    for (int row = 0; row < T1a.params->rowtot[h]; row++)
        for (int col = 0; col < T1a.params->coltot[h]; col++)
            error[0][word++] = T1a.matrix[h][row][col] - T1b.matrix[h][row][col];
```

3. **B matrix and solve** (lines 166-215):
```cpp
for (int p = 0; p < nvector; p++) {
    for (int q = 0; q < p; q++) {
        product = C_DDOT(vector_length, vector[1], 1, vector[0], 1);
        B[p][q] = B[q][p] = product;
    }
}
diis_invert_B(B, C, nvector + 1, 1.0E-12);
```

**Deliverable:** Document with key insights about current implementation

---

### Step 1.2: Set Up Development Branch

**Task:** Create isolated development environment

```bash
# Create feature branch for POC
cd /home/user/psi4
git checkout claude/consolidate-diis-implementations-01Uw9XohC6D9jFN2riVN56EZ
git checkout -b poc/ccenergy-rhf-diis-libdiis

# Verify clean state
git status
```

**Checklist:**
- [ ] Branch created
- [ ] Based on latest consolidation work
- [ ] Clean working directory

---

### Step 1.3: Identify Test Cases

**Task:** Find regression tests that exercise RHF CCSD

```bash
# Search for CCSD tests
find psi4/tests -name "*.in" -type f | xargs grep -l "ccsd"
find psi4/tests -name "*.in" -type f | xargs grep -l "CCSD"

# Look specifically for RHF CCSD
find psi4/tests -name "*.in" -type f | xargs grep -l "reference.*rhf" | xargs grep -l "ccsd"
```

**Required Test Coverage:**
- [ ] Basic RHF CCSD (small molecule, e.g., H2O)
- [ ] Multiple iterations (tests DIIS subspace)
- [ ] Convergence test
- [ ] Energy accuracy test
- [ ] Different basis sets (if available)

**Example Tests to Identify:**
- `cc1` - Basic CCSD
- `cc2` - CCSD with different settings
- `cc4a` - CCSD gradient (uses DIIS heavily)
- `cc13` - CCSD(T)
- `cc54` - DF-CCSD

**Deliverable:** List of 5-10 test cases for validation

---

## Phase 2: Implementation (3-4 hours)

### Step 2.1: Create New Implementation File

**Task:** Add new file alongside existing implementation

```bash
# Create new file
touch psi4/src/psi4/cc/ccenergy/diis_RHF_libdiis.cc
```

**File Header:**
```cpp
/*
 * @BEGIN LICENSE
 * [Standard Psi4 license header]
 * @END LICENSE
 */

/*! \file
    \ingroup CCENERGY
    \brief RHF DIIS implementation using libdiis/DIISManager (POC)

    This is a proof-of-concept implementation to validate that
    libdiis can successfully replace the custom DIIS implementation
    for RHF CCSD amplitudes. If successful, this approach will be
    extended to ROHF and UHF references.

    Author: [Your name]
    Date: 2025-11-18
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/libdiis/diismanager.h"  // NEW: libdiis interface
#include "MOInfo.h"
#include "Params.h"
#include "psi4/cc/ccwave.h"

namespace psi {
namespace ccenergy {

void CCEnergyWavefunction::diis_RHF_libdiis(int iter) {
    // Implementation goes here
}

}  // namespace ccenergy
}  // namespace psi
```

**Checklist:**
- [ ] File created
- [ ] License header added
- [ ] Includes added
- [ ] Namespace structure correct
- [ ] Function signature matches existing `diis_RHF()`

---

### Step 2.2: Initialize DIISManager

**Task:** Set up DIISManager as class member

**Modify `psi4/src/psi4/cc/ccwave.h`:**

```cpp
// Around line 200-300, in protected members section
class CCEnergyWavefunction : public Wavefunction {
    // ... existing members ...

protected:
    // ... existing members ...

#ifdef USE_LIBDIIS_POC
    /// DIISManager for amplitude extrapolation (POC)
    std::shared_ptr<DIISManager> ccsd_diis_manager_;
#endif

    // ... rest of class ...
};
```

**Add initialization in constructor:**

Locate the CCEnergyWavefunction constructor and add:

```cpp
CCEnergyWavefunction::CCEnergyWavefunction(/* args */) {
    // ... existing initialization ...

#ifdef USE_LIBDIIS_POC
    // Initialize DIIS manager
    // Using same settings as original: 8 vectors, largest error removal, on disk
    ccsd_diis_manager_ = std::make_shared<DIISManager>(
        8,                                          // max_vecs (same as nvector in original)
        "CCSD DIIS RHF",                           // label
        DIISManager::RemovalPolicy::LargestError,  // same as original
        DIISManager::StoragePolicy::OnDisk         // same as original
    );
    // Note: No need to call set_error_vector_size or set_vector_size
    // libdiis auto-detects from first add_entry() call
#endif
}
```

**Checklist:**
- [ ] Header file modified
- [ ] DIISManager member added (conditional)
- [ ] Initialization added to constructor
- [ ] Compile flag defined (`USE_LIBDIIS_POC`)

---

### Step 2.3: Implement Core DIIS Logic

**Task:** Write the new `diis_RHF_libdiis()` function

**Implementation in `diis_RHF_libdiis.cc`:**

```cpp
void CCEnergyWavefunction::diis_RHF_libdiis(int iter) {
    // Early exit if not enough iterations
    if (iter < 2) {
        return;  // Need at least 2 iterations for DIIS
    }

    auto nirreps = moinfo_.nirreps;
    dpdfile2 T1_new, T1_old, R1;
    dpdbuf4 T2_new, T2_old, R2;

    // ==========================================
    // Step 1: Compute error vectors
    // ==========================================

    // Initialize error vector for T1: R1 = T1_new - T1_old
    global_dpd_->file2_init(&R1, PSIF_CC_OEI, 0, 0, 1, "R1_IA");
    global_dpd_->file2_init(&T1_new, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_init(&T1_old, PSIF_CC_TAMPS, 0, 0, 1, "tIA");

    // Compute R1 = T1_new - T1_old using DPD operations
    global_dpd_->file2_copy(&T1_new, PSIF_CC_OEI, "R1_IA");
    global_dpd_->file2_close(&T1_new);
    global_dpd_->file2_init(&R1, PSIF_CC_OEI, 0, 0, 1, "R1_IA");
    global_dpd_->file2_axpy(&T1_old, &R1, -1.0, 0);  // R1 = R1 + (-1.0) * T1_old
    global_dpd_->file2_close(&T1_old);

    // Initialize error vector for T2: R2 = T2_new - T2_old
    global_dpd_->buf4_init(&R2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "R2_IjAb");
    global_dpd_->buf4_init(&T2_new, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_init(&T2_old, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");

    // Compute R2 = T2_new - T2_old
    global_dpd_->buf4_copy(&T2_new, PSIF_CC_TAMPS, "R2_IjAb");
    global_dpd_->buf4_close(&T2_new);
    global_dpd_->buf4_init(&R2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "R2_IjAb");
    global_dpd_->buf4_axpy(&T2_old, &R2, -1.0);
    global_dpd_->buf4_close(&T2_old);

    // ==========================================
    // Step 2: Reopen amplitudes for extrapolation
    // ==========================================

    global_dpd_->file2_init(&T1_new, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->buf4_init(&T2_new, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    // ==========================================
    // Step 3: Add to DIIS subspace
    // ==========================================

    // libdiis expects error vectors first, then amplitude vectors
    bool added = ccsd_diis_manager_->add_entry(&R1, &R2, &T1_new, &T2_new);

    if (!added) {
        outfile->Printf("  Warning: DIIS add_entry failed\n");
    }

    // ==========================================
    // Step 4: Extrapolate if enough vectors
    // ==========================================

    if (ccsd_diis_manager_->subspace_size() >= 2) {
        // Perform extrapolation - updates T1_new and T2_new in place
        ccsd_diis_manager_->extrapolate(&T1_new, &T2_new);

        outfile->Printf("  DIIS: extrapolated using %d vectors\n",
                       ccsd_diis_manager_->subspace_size());
    } else {
        outfile->Printf("  DIIS: stored vector (subspace size = %d)\n",
                       ccsd_diis_manager_->subspace_size());
    }

    // ==========================================
    // Step 5: Cleanup
    // ==========================================

    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->file2_close(&T1_new);
    global_dpd_->buf4_close(&T2_new);
}
```

**Key Implementation Notes:**

1. **DPD Operations Used:**
   - `file2_copy()` - Copy file2 to new location
   - `file2_axpy()` - Compute Y = Y + alpha*X for file2
   - `buf4_copy()` - Copy buf4 to new location
   - `buf4_axpy()` - Compute Y = Y + alpha*X for buf4

2. **libdiis Handles:**
   - Vector flattening (DPD → Matrix internally)
   - PSIO storage management
   - B matrix construction
   - Linear system solving with conditioning
   - Amplitude extrapolation

3. **What We Eliminated:**
   - Manual vector length calculation (~10 lines)
   - Manual flattening to 1D arrays (~50 lines)
   - Manual PSIO read/write (~30 lines)
   - Manual B matrix construction (~40 lines)
   - Manual extrapolation loop (~30 lines)
   - Total: ~160 lines → ~60 lines

**Checklist:**
- [ ] Error vector computation implemented
- [ ] DIIS manager add_entry call added
- [ ] Extrapolation logic added
- [ ] Cleanup code added
- [ ] Comments explain each section

---

### Step 2.4: Add Compile-Time Switching

**Task:** Allow easy A/B testing between implementations

**Modify `psi4/src/psi4/cc/ccenergy/diis.cc`:**

```cpp
void CCEnergyWavefunction::diis(int iter) {
#ifdef USE_LIBDIIS_POC
    // Use new libdiis implementation (POC)
    if (params_.ref == 0) {
        diis_RHF_libdiis(iter);
    } else if (params_.ref == 1) {
        diis_ROHF(iter);  // Original implementation
    } else if (params_.ref == 2) {
        diis_UHF(iter);   // Original implementation
    }
#else
    // Use original implementation
    if (params_.ref == 0) {
        diis_RHF(iter);
    } else if (params_.ref == 1) {
        diis_ROHF(iter);
    } else if (params_.ref == 2) {
        diis_UHF(iter);
    }
#endif
}
```

**Add function declaration in `psi4/src/psi4/cc/ccwave.h`:**

```cpp
class CCEnergyWavefunction : public Wavefunction {
    // ... existing declarations ...

    void diis_RHF(int iter);
#ifdef USE_LIBDIIS_POC
    void diis_RHF_libdiis(int iter);  // POC implementation
#endif
    void diis_ROHF(int iter);
    void diis_UHF(int iter);

    // ...
};
```

**Update `psi4/src/psi4/cc/ccenergy/CMakeLists.txt`:**

```cmake
list(APPEND sources
  # ... existing sources ...
  diis.cc
  diis_RHF.cc
  diis_ROHF.cc
  diis_UHF.cc
  # POC: Add conditionally
  $<$<BOOL:${USE_LIBDIIS_POC}>:diis_RHF_libdiis.cc>
  # ... rest of sources ...
)
```

**Checklist:**
- [ ] Dispatcher modified for conditional compilation
- [ ] Function declaration added
- [ ] CMakeLists.txt updated
- [ ] Can toggle with compile flag

---

### Step 2.5: Build Configuration

**Task:** Set up build system for POC

**Create build script `build_poc.sh`:**

```bash
#!/bin/bash
# Build script for ccenergy DIIS POC

set -e

BUILD_DIR="build_diis_poc"
INSTALL_PREFIX="/tmp/psi4_diis_poc"

# Clean previous build
rm -rf ${BUILD_DIR}
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# Configure with POC flag enabled
cmake .. \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
    -DCMAKE_CXX_FLAGS="-DUSE_LIBDIIS_POC" \
    -DENABLE_OPENMP=ON \
    -DENABLE_XHOST=OFF \
    -DMAX_AM_ERI=6

# Build ccenergy module only (faster iteration)
make ccenergy -j4

echo "POC build complete!"
echo "To test: cd ${BUILD_DIR} && ctest -R cc"
```

**Make executable:**
```bash
chmod +x build_poc.sh
```

**Build both versions for comparison:**

```bash
# Build POC version
./build_poc.sh

# Build original version (for comparison)
mkdir -p build_original
cd build_original
cmake .. -DCMAKE_BUILD_TYPE=Debug -DENABLE_OPENMP=ON
make ccenergy -j4
cd ..
```

**Checklist:**
- [ ] Build script created
- [ ] POC version builds successfully
- [ ] Original version builds successfully
- [ ] Both available for testing

---

## Phase 3: Testing & Validation (2-3 hours)

### Step 3.1: Basic Functionality Test

**Task:** Verify POC doesn't crash

**Create minimal test file `test_diis_poc.py`:**

```python
#!/usr/bin/env python3
"""
Minimal test for ccenergy DIIS POC
"""

import psi4

# Simple RHF CCSD calculation
psi4.set_memory('1 GB')

mol = psi4.geometry("""
O
H 1 0.96
H 1 0.96 2 104.5
symmetry c1
""")

psi4.set_options({
    'basis': 'cc-pvdz',
    'scf_type': 'pk',
    'mp2_type': 'conv',
    'freeze_core': 'false',
    'e_convergence': 1e-8,
    'd_convergence': 1e-8,
    'r_convergence': 1e-7,
    'diis_max_vecs': 8,
})

# Run CCSD
e_ccsd = psi4.energy('ccsd')

print(f"\nCCSD Energy: {e_ccsd:.10f}")
print("If this printed, POC didn't crash!")
```

**Run test:**
```bash
cd build_diis_poc
export PYTHONPATH=$(pwd)/stage/lib:$PYTHONPATH
python3 ../test_diis_poc.py
```

**Success Criteria:**
- [ ] Calculation completes without errors
- [ ] DIIS messages printed
- [ ] Energy value obtained

---

### Step 3.2: Numerical Validation

**Task:** Compare energies between implementations

**Create comparison script `compare_implementations.py`:**

```python
#!/usr/bin/env python3
"""
Compare ccenergy DIIS implementations
"""

import psi4
import numpy as np

def run_ccsd(label, build_dir):
    """Run CCSD and return energy"""

    psi4.core.clean()
    psi4.set_memory('1 GB')

    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)

    psi4.set_options({
        'basis': 'cc-pvdz',
        'scf_type': 'pk',
        'freeze_core': 'false',
        'e_convergence': 1e-10,
        'd_convergence': 1e-10,
        'r_convergence': 1e-8,
    })

    e = psi4.energy('ccsd')
    print(f"{label}: {e:.12f}")
    return e

# Run both
print("Running Original Implementation...")
e_original = run_ccsd("Original", "build_original")

print("\nRunning POC Implementation...")
e_poc = run_ccsd("POC", "build_diis_poc")

# Compare
diff = abs(e_original - e_poc)
print(f"\nDifference: {diff:.2e}")

if diff < 1e-9:
    print("✅ PASS: Energies match to 1e-9 Eh")
elif diff < 1e-6:
    print("⚠️  WARNING: Small difference detected")
else:
    print("❌ FAIL: Significant difference!")

print(f"Relative error: {diff/abs(e_original):.2e}")
```

**Run comparison:**
```bash
python3 compare_implementations.py
```

**Success Criteria:**
- [ ] Energy difference < 1e-9 Hartree
- [ ] Relative error < 1e-11
- [ ] Both calculations converge in similar iterations

---

### Step 3.3: Regression Test Suite

**Task:** Run full ccenergy test suite

**Run tests with POC:**
```bash
cd build_diis_poc
ctest -R cc1 -V        # Basic CCSD
ctest -R cc2 -V        # CCSD variations
ctest -R cc4a -V       # CCSD gradient
ctest -R cc13 -V       # CCSD(T)
ctest -R cc54 -V       # DF-CCSD
```

**Run tests with original:**
```bash
cd ../build_original
ctest -R cc1 -V
ctest -R cc2 -V
ctest -R cc4a -V
ctest -R cc13 -V
ctest -R cc54 -V
```

**Create test comparison script `compare_tests.sh`:**

```bash
#!/bin/bash
# Compare test results between implementations

echo "Running tests with ORIGINAL implementation..."
cd build_original
ctest -R "cc(1|2|4a|13|54)" --output-on-failure > ../original_tests.log 2>&1
cd ..

echo "Running tests with POC implementation..."
cd build_diis_poc
ctest -R "cc(1|2|4a|13|54)" --output-on-failure > ../poc_tests.log 2>&1
cd ..

echo ""
echo "Comparing results..."
diff -u original_tests.log poc_tests.log || echo "Differences found - review manually"

# Extract energies from both
grep -i "ccsd.*energy" original_tests.log > original_energies.txt
grep -i "ccsd.*energy" poc_tests.log > poc_energies.txt

echo ""
echo "Energy comparison:"
paste original_energies.txt poc_energies.txt
```

**Success Criteria:**
- [ ] All tests pass with POC
- [ ] Same tests pass/fail in both
- [ ] Energies match to specified precision
- [ ] No new failures introduced

---

### Step 3.4: Convergence Behavior Analysis

**Task:** Verify DIIS is working correctly

**Create convergence analysis script `analyze_convergence.py`:**

```python
#!/usr/bin/env python3
"""
Analyze CCSD convergence behavior
"""

import psi4
import re

def run_and_analyze(label, psi4_path):
    """Run CCSD and extract iteration data"""

    psi4.core.clean()
    psi4.set_memory('1 GB')

    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    psi4.set_options({
        'basis': 'cc-pvdz',
        'freeze_core': 'false',
        'e_convergence': 1e-10,
        'd_convergence': 1e-10,
        'r_convergence': 1e-8,
        'cc_maxiter': 50,
    })

    # Capture output
    psi4.core.set_output_file(f'output_{label}.dat', False)

    e = psi4.energy('ccsd')

    # Parse output for iteration data
    with open(f'output_{label}.dat', 'r') as f:
        output = f.read()

    # Extract iteration data
    iterations = []
    for line in output.split('\n'):
        # Look for iteration lines
        match = re.search(r'\s+(\d+)\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)', line)
        if match:
            iter_num = int(match.group(1))
            e_corr = float(match.group(2))
            de = float(match.group(3))
            t1_rms = float(match.group(4))
            t2_rms = float(match.group(5))
            iterations.append((iter_num, e_corr, de, t1_rms, t2_rms))

    return e, iterations

# Compare convergence
print("Analyzing convergence patterns...\n")

e_orig, iter_orig = run_and_analyze("original", "build_original/stage")
e_poc, iter_poc = run_and_analyze("poc", "build_diis_poc/stage")

print(f"Original converged in {len(iter_orig)} iterations")
print(f"POC converged in {len(iter_poc)} iterations")
print(f"Energy difference: {abs(e_orig - e_poc):.2e}\n")

# Compare iteration-by-iteration
print("Iteration-by-iteration comparison:")
print(f"{'Iter':>4} {'Original DE':>15} {'POC DE':>15} {'Difference':>15}")
print("-" * 52)

for i in range(min(len(iter_orig), len(iter_poc))):
    orig_de = iter_orig[i][2]
    poc_de = iter_poc[i][2]
    diff = abs(orig_de - poc_de)
    print(f"{i+1:4d} {orig_de:15.8e} {poc_de:15.8e} {diff:15.8e}")
```

**Success Criteria:**
- [ ] Same number of iterations to converge
- [ ] Similar energy change per iteration
- [ ] DIIS kicks in at same iteration
- [ ] Final energies match

---

### Step 3.5: Performance Benchmarking

**Task:** Ensure no performance regression

**Create benchmark script `benchmark_diis.py`:**

```python
#!/usr/bin/env python3
"""
Benchmark DIIS performance
"""

import psi4
import time

def benchmark(label, num_runs=3):
    """Run CCSD multiple times and average"""

    times = []

    for run in range(num_runs):
        psi4.core.clean()
        psi4.set_memory('2 GB')

        # Larger molecule for meaningful timing
        mol = psi4.geometry("""
        O
        H 1 0.96
        H 1 0.96 2 104.5
        O 1 2.5 2 120.0 3 90.0
        H 4 0.96 1 120.0 2 0.0
        H 4 0.96 5 104.5 1 120.0
        """)

        psi4.set_options({
            'basis': 'cc-pvdz',
            'freeze_core': 'true',
            'e_convergence': 1e-8,
        })

        start = time.time()
        e = psi4.energy('ccsd')
        elapsed = time.time() - start

        times.append(elapsed)
        print(f"{label} run {run+1}: {elapsed:.2f}s")

    avg_time = sum(times) / len(times)
    std_time = (sum((t - avg_time)**2 for t in times) / len(times)) ** 0.5

    return avg_time, std_time

print("Benchmarking Original Implementation...")
orig_avg, orig_std = benchmark("Original")

print("\nBenchmarking POC Implementation...")
poc_avg, poc_std = benchmark("POC")

print("\n" + "="*50)
print(f"Original: {orig_avg:.2f} ± {orig_std:.2f}s")
print(f"POC:      {poc_avg:.2f} ± {poc_std:.2f}s")

overhead = (poc_avg - orig_avg) / orig_avg * 100
print(f"Overhead: {overhead:+.1f}%")

if abs(overhead) < 5.0:
    print("✅ PASS: Performance within 5%")
elif overhead < 10.0:
    print("⚠️  WARNING: Small performance regression")
else:
    print("❌ FAIL: Significant performance change")
```

**Success Criteria:**
- [ ] POC within 5% of original performance
- [ ] No major memory increase
- [ ] DIIS overhead acceptable

---

## Phase 4: Documentation & Review (1 hour)

### Step 4.1: Document Results

**Task:** Create comprehensive results report

**Create `POC_RESULTS.md`:**

```markdown
# ccenergy RHF DIIS POC Results

## Test Date
[Date]

## Implementation Summary
- Original lines: 258
- POC lines: ~60
- Reduction: 76.7%

## Numerical Validation

### Energy Comparison
| Test | Original (Eh) | POC (Eh) | Difference |
|------|---------------|----------|------------|
| H2O/cc-pVDZ | [value] | [value] | [diff] |
| ... | | | |

### Convergence Comparison
| Test | Original Iters | POC Iters | Match? |
|------|---------------|-----------|--------|
| H2O/cc-pVDZ | [n] | [n] | ✅/❌ |

## Performance Results

| Test | Original (s) | POC (s) | Overhead |
|------|-------------|---------|----------|
| H2O | [time] | [time] | [%] |

## Regression Tests

| Test | Original | POC | Status |
|------|----------|-----|--------|
| cc1 | PASS | PASS | ✅ |
| cc2 | PASS | PASS | ✅ |
| cc4a | PASS | PASS | ✅ |

## Conclusions

[Summary of findings]

## Recommendations

[Next steps]
```

**Fill in with actual results**

**Checklist:**
- [ ] All test results documented
- [ ] Energy comparisons included
- [ ] Performance data included
- [ ] Conclusions written

---

### Step 4.2: Code Review Prep

**Task:** Prepare code for review

**Review checklist:**
- [ ] Code follows Psi4 style guide
- [ ] Comments explain non-obvious parts
- [ ] No debug print statements left
- [ ] Error handling appropriate
- [ ] Memory management correct (DPD close calls)

**Create diff for review:**
```bash
git diff claude/consolidate-diis-implementations-01Uw9XohC6D9jFN2riVN56EZ > poc_changes.diff
```

**Key files to review:**
- `psi4/src/psi4/cc/ccenergy/diis_RHF_libdiis.cc` (new)
- `psi4/src/psi4/cc/ccenergy/diis.cc` (modified)
- `psi4/src/psi4/cc/ccwave.h` (modified)
- `psi4/src/psi4/cc/ccenergy/CMakeLists.txt` (modified)

---

### Step 4.3: Create Summary Presentation

**Task:** Prepare decision briefing

**Slides to create (`POC_PRESENTATION.md`):**

1. **Objective**
   - Validate libdiis for cc modules
   - Reduce code duplication
   - Maintain numerical accuracy

2. **Approach**
   - Implement RHF only (smallest scope)
   - Parallel implementation (safe rollback)
   - Comprehensive testing

3. **Technical Details**
   - Before/after code comparison
   - How libdiis handles DPD buffers
   - Key implementation decisions

4. **Results**
   - Numerical validation (energy plots)
   - Convergence comparison (iteration plots)
   - Performance benchmarks (timing charts)
   - Test suite results

5. **Code Reduction**
   - 258 lines → 60 lines
   - Eliminated manual buffer management
   - Simplified extrapolation logic

6. **Recommendations**
   - ✅ Proceed to Phase 2 (ROHF/UHF)
   - ⚠️ Address any issues found
   - ❌ Abandon if fundamental problems

---

## Phase 5: Decision & Next Steps (30 min)

### Step 5.1: Go/No-Go Decision

**Decision Criteria:**

**GO if:**
- ✅ All energies match to 1e-9 Eh
- ✅ All regression tests pass
- ✅ Performance within 5% of original
- ✅ Convergence behavior identical
- ✅ Code reduction achieved (~75%)

**NO-GO if:**
- ❌ Any numerical discrepancies > 1e-6
- ❌ Regression test failures
- ❌ Performance regression > 10%
- ❌ Convergence issues
- ❌ Implementation problems

**REVISE if:**
- ⚠️ Minor issues identified but fixable
- ⚠️ Performance needs optimization
- ⚠️ Edge cases need handling

---

### Step 5.2: If GO - Plan Phase 2

**Next Steps:**
1. Implement `diis_ROHF_libdiis()` (same pattern)
2. Implement `diis_UHF_libdiis()` (same pattern)
3. Run expanded test suite
4. Performance optimization if needed
5. Prepare for integration

**Timeline:** 1-2 weeks

---

### Step 5.3: If NO-GO - Document Blockers

**Document:**
- What specifically failed
- Why libdiis isn't suitable
- Alternative approaches
- Lessons learned

**Possible blockers:**
- DPD conversion overhead too high
- Numerical stability issues
- Missing features in libdiis
- Integration complexity

---

## Timeline Summary

| Phase | Duration | Tasks |
|-------|----------|-------|
| 1. Analysis & Setup | 2 hours | Understand code, create branch, identify tests |
| 2. Implementation | 3-4 hours | Write new code, build system, integration |
| 3. Testing | 2-3 hours | Validation, benchmarking, regression tests |
| 4. Documentation | 1 hour | Results report, code review prep |
| 5. Decision | 30 min | Go/no-go evaluation |
| **TOTAL** | **8-10 hours** | **Complete POC** |

---

## Risk Mitigation

### Technical Risks

| Risk | Mitigation |
|------|------------|
| DPD buffer conversion issues | Test with small molecules first |
| Numerical precision loss | Bit-for-bit comparison |
| Performance regression | Profile and optimize |
| Edge case failures | Comprehensive test suite |

### Implementation Risks

| Risk | Mitigation |
|------|------------|
| Breaking existing code | Parallel implementation with flag |
| Build system issues | Separate build directory |
| Test failures | Easy rollback via git |
| Integration problems | POC in isolated branch |

---

## Success Metrics

### Must Have (Hard Requirements)
- ✅ Energy accuracy < 1e-9 Eh
- ✅ All regression tests pass
- ✅ Convergence identical
- ✅ Code compiles cleanly
- ✅ No memory leaks

### Should Have (Soft Requirements)
- 📊 Performance within 5% of original
- 📊 Code reduction > 70%
- 📊 Clear documentation
- 📊 Easy to extend to ROHF/UHF

### Nice to Have (Bonus)
- 🎯 Performance improvement
- 🎯 Better error messages
- 🎯 Simplified debugging

---

## Deliverables Checklist

At end of POC, you should have:

- [ ] `diis_RHF_libdiis.cc` - New implementation
- [ ] Modified `diis.cc` - Dispatcher with flag
- [ ] Modified `ccwave.h` - New function declaration
- [ ] Modified `CMakeLists.txt` - Build integration
- [ ] `POC_RESULTS.md` - Results report
- [ ] `test_diis_poc.py` - Basic test
- [ ] `compare_implementations.py` - Comparison script
- [ ] `analyze_convergence.py` - Convergence analysis
- [ ] `benchmark_diis.py` - Performance benchmark
- [ ] Energy comparison data
- [ ] Test suite results
- [ ] Performance benchmarks
- [ ] Decision recommendation

---

## References

### Documentation
- CC DIIS Analysis: `CC_DIIS_CONSOLIDATION_ANALYSIS.md`
- libdiis interface: `psi4/src/psi4/libdiis/diismanager.h`
- DPD documentation: Psi4 wiki

### Code References
- Original implementation: `psi4/src/psi4/cc/ccenergy/diis_RHF.cc`
- libdiis Python backend: `psi4/driver/procrouting/diis.py`
- occ module example: `psi4/src/psi4/occ/occ_iterations.cc:501`

### Tests
- Basic CCSD: `psi4/tests/cc1`
- CCSD variations: `psi4/tests/cc2`
- CCSD gradient: `psi4/tests/cc4a`

---

## Notes

### Common Pitfalls
- Forgetting to close DPD files → memory leak
- Wrong PSIO file numbers → corrupted data
- Incorrect irrep handling → wrong energies
- Missing initialization → crash

### Debugging Tips
- Add `outfile->Printf()` for diagnostics
- Use `global_dpd_->file2_print()` to inspect arrays
- Check DIIS subspace size at each iteration
- Compare B matrices between implementations

### Performance Optimization Opportunities
- Reduce DPD conversions if needed
- Optimize error vector computation
- Consider InCore storage for small systems

---

*POC Plan prepared: 2025-11-18*
*Estimated total time: 8-10 hours*
*Expected outcome: Validation of libdiis for cc modules*
