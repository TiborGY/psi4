# fnocc DIIS Consolidation Analysis

**Date**: 2025-11-18
**Investigator**: Claude (Anthropic AI)
**Context**: Following successful ccenergy DIIS migration to libdiis

---

## Executive Summary

The fnocc module contains a custom DIIS implementation (~248 lines) that shares conceptual similarity with the ccenergy DIIS code we just consolidated. However, unlike ccenergy which uses DPD buffers, fnocc uses raw `double*` arrays.

**Key Question**: Can fnocc DIIS be migrated to libdiis like ccenergy?

**Short Answer**: **Partially** - Migration is technically feasible but may introduce performance overhead and requires wrapping raw arrays in Matrix objects.

**Recommendation**: **Lower Priority** - The benefit-to-effort ratio is less favorable than ccenergy migration due to data structure differences.

---

## Current Implementation

### File Structure
- **Location**: `psi4/src/psi4/fnocc/diis.cc`
- **Size**: 248 lines
- **Class**: `CoupledCluster` (defined in `fnocc/ccsd.h`)

### DIIS Functions (4 total)

1. **`DIIS(double* c, long int nvec, long int n, int replace_diis_iter)`** - 127 lines
   - Builds B matrix from stored error vector dot products
   - Solves linear system using LAPACK `DGESV`
   - Returns DIIS coefficients in array `c`

2. **`DIISOldVector(long int iter, int diis_iter, int replace_diis_iter)`** - 37 lines
   - Stores current amplitudes (T1 and T2) to PSIO
   - Manages circular buffer when subspace is full

3. **`DIISErrorVector(int diis_iter, int replace_diis_iter, int iter)`** - 35 lines
   - Computes and stores error vector
   - Returns error norm for convergence checking
   - Error vectors stored separately from amplitude vectors

4. **`DIISNewAmplitudes(int diis_iter, int& replace_diis_iter)`** - 47 lines
   - Reads stored amplitude vectors from PSIO
   - Computes linear combination using DIIS coefficients
   - Updates T1 and T2 with extrapolated values

### Data Structures

**Amplitudes (Raw Arrays)**:
```cpp
double* tb;  // T2 amplitudes: size = o*o*v*v
double* t1;  // T1 amplitudes: size = o*v
```

**Workspace Buffers**:
```cpp
double* tempt;   // Temporary storage for T amplitudes
double* tempv;   // Temporary storage for error vectors
```

**DIIS Coefficients**:
```cpp
double* diisvec;  // DIIS linear combination coefficients
long int maxdiis; // Maximum DIIS subspace size
```

### Storage Mechanism

**Uses Manual PSIO Management**:
- Error vectors: `PSIF_DCC_EVEC` file
  - Entries: `"evector1"`, `"evector2"`, ..., `"evectorN"`
  - Also stores `"error matrix"` (B matrix components)
- Amplitude vectors: `PSIF_DCC_OVEC` file
  - Entries: `"oldvector1"`, `"oldvector2"`, ..., `"oldvectorN"`

### Key Characteristics

1. **Flat Array Storage**: All data in contiguous `double*` arrays
2. **No Symmetry Blocking**: Unlike DPD, no irrep structure
3. **Manual Memory Management**: Explicit `malloc`/`free` calls
4. **BLAS Operations**: Uses `C_DDOT`, `C_DAXPY`, `C_DNRM2`, `C_DCOPY`
5. **LAPACK Solving**: Direct `DGESV` call for linear system

---

## Comparison with ccenergy

| Aspect | ccenergy (Before) | ccenergy (After libdiis) | fnocc (Current) |
|--------|-------------------|-------------------------|-----------------|
| **Data Structure** | DPD buffers (dpdbuf4, dpdfile2) | DPD buffers | Raw double* arrays |
| **Symmetry** | Yes (irrep blocking) | Yes (via DPD) | No (flat arrays) |
| **Storage** | PSIO + manual management | Automatic (libdiis) | PSIO + manual management |
| **B Matrix** | Manual construction | Automatic (libdiis) | Manual construction |
| **Linear Solve** | Manual DGESV | Automatic (libdiis) | Manual DGESV |
| **Code Lines** | ~258 (RHF) | ~60 (RHF) | ~248 (total) |
| **libdiis Support** | Native (DPD) | Yes | Requires wrapping |

---

## libdiis Compatibility Analysis

### What libdiis Supports

From `libdiis/diismanager.h` and `driver/procrouting/diis.py`:

✅ **Native Support**:
- `core.Matrix` - 2D matrices
- `core.Vector` - 1D vectors
- `core.dpdbuf4` - 4-index DPD buffers (converted to Matrix)
- `core.dpdfile2` - 2-index DPD files (converted to Matrix)
- `ambit.BlockedTensor` - Ambit tensors (if available)

❌ **NOT Natively Supported**:
- Raw `double*` arrays
- Custom C-style memory layouts

### Migration Approaches

#### Option 1: Wrap Arrays in Matrix Objects (Feasible but Overhead)

**Concept**: Create Matrix wrappers around existing double* arrays

**Implementation**:
```cpp
// Wrap T1 amplitudes
auto T1_matrix = std::make_shared<Matrix>("T1", o, v);
T1_matrix->set_pointer(t1);  // If such method exists

// Wrap T2 amplitudes (reshaped as 2D)
auto T2_matrix = std::make_shared<Matrix>("T2", o*o, v*v);
T2_matrix->set_pointer(tb);

// Add to DIIS
ccsd_diis_manager_->add_entry(error_vectors..., T1_matrix, T2_matrix);
```

**Issues**:
- Matrix class may not have `set_pointer()` method
- May require copying data: `double* → Matrix` (overhead)
- Need to copy back after extrapolation: `Matrix → double*`

#### Option 2: Copy Arrays to/from Matrix (Simple but Costly)

**Concept**: Copy data between double* and Matrix objects

**Implementation**:
```cpp
// Before DIIS
auto T1_matrix = std::make_shared<Matrix>("T1", o, v);
for (int i = 0; i < o; ++i)
    for (int a = 0; a < v; ++a)
        T1_matrix->set(i, a, t1[i*v + a]);

// After DIIS extrapolation
for (int i = 0; i < o; ++i)
    for (int a = 0; a < v; ++a)
        t1[i*v + a] = T1_matrix->get(i, a);
```

**Issues**:
- **Performance overhead**: O(n) copying for each DIIS call
- For large systems: T2 has o²v² elements (can be millions)
- Copy overhead may negate DIIS acceleration benefits

#### Option 3: Extend libdiis to Support Raw Arrays (Complex)

**Concept**: Add native support for double* to libdiis Python backend

**Implementation**:
- Modify `diis.py` to accept numpy arrays via pybind11
- Create C++ → Python bridge for double* arrays
- Handle memory management carefully

**Issues**:
- Requires modifying libdiis core
- Adds complexity to libdiis for a single module
- Maintenance burden

---

## Duplication Assessment

### Code Duplication: **MODERATE**

**Duplicated Concepts** (shared with ccenergy and libdiis):
- B matrix construction from error vector dot products
- DIIS linear system solving (Ax = b)
- Subspace management (circular buffer when full)
- Error vector storage and retrieval

**Unique to fnocc**:
- Raw array handling (vs DPD buffers)
- Flat memory layout (no symmetry blocking)
- Specific PSIO file naming scheme

### Estimated Code Reduction

**Best Case** (with efficient wrapping):
- Current: ~248 lines
- Potential: ~80-100 lines (60% reduction)
- **Savings**: ~150 lines

**Realistic** (with copying overhead):
- Current: ~248 lines
- Potential: ~100-120 lines (50% reduction)
- **Savings**: ~120-130 lines

Compare to ccenergy:
- RHF: 258 → 60 lines (83% reduction, ~200 lines saved)
- ROHF: 357 → 90 lines (75% reduction, ~267 lines saved)
- UHF: 365 → 110 lines (70% reduction, ~255 lines saved)
- **Total ccenergy**: ~720 lines saved

---

## Feasibility Analysis

### Technical Feasibility: **MEDIUM**

✅ **Pros**:
- DIIS algorithm is identical conceptually
- libdiis already handles similar operations
- Could provide code simplification and maintainability

❌ **Cons**:
- Requires data structure conversion (double* → Matrix)
- Performance overhead from copying
- fnocc uses flat arrays for performance (CC calculations are expensive)
- Matrix class may not support zero-copy wrapping of external pointers

### Performance Concerns: **MODERATE to HIGH**

**Critical Path**:
- CCSD iterations are already expensive (O(N⁶) scaling)
- DIIS called every iteration
- Copying large T2 arrays (o²v²) could be noticeable

**Example System**:
- Medium molecule: o=20, v=100
- T2 size: 20² × 100² = 40,000,000 elements
- Copy time (estimate): ~100ms per iteration
- DIIS calls per convergence: ~15-20 iterations
- Total overhead: ~2 seconds per CCSD calculation

**Impact**:
- For small molecules: Negligible
- For large molecules: Could add 5-10% to runtime

### Maintenance Benefit: **LOW to MEDIUM**

**Benefits**:
- Eliminates ~150 lines of custom DIIS code
- Centralizes DIIS in libdiis (single point of maintenance)
- Consistent DIIS behavior across modules

**Drawbacks**:
- fnocc is a specialized module (frozen natural orbitals CC)
- Less actively developed than ccenergy
- Code is already working and stable
- Migration effort may not be justified

---

## Recommendation

### Priority: **LOW**

**Rationale**:
1. **Modest Code Reduction**: ~150 lines vs ~720 for ccenergy
2. **Performance Risk**: Copying overhead for large systems
3. **Data Structure Mismatch**: Raw arrays vs Matrix objects
4. **Module Maturity**: fnocc is stable and less actively developed
5. **Effort vs Benefit**: Migration effort similar to ccenergy but lower payoff

### Suggested Approach

**Option A: Leave As-Is (Recommended)**
- fnocc DIIS is working correctly
- No known bugs or maintenance issues
- Performance is already optimized for flat arrays
- ~248 lines is manageable for a specialized module

**Option B: Document Only**
- Add comment referencing libdiis as conceptual template
- Note that migration was considered but deferred due to data structure differences
- Keeps door open for future migration if libdiis adds zero-copy array support

**Option C: Defer to Future**
- Wait for potential libdiis enhancements (e.g., native array support)
- Revisit if fnocc undergoes major refactoring
- Consider if performance profiling shows DIIS is a bottleneck

---

## Alternative: Focus on Higher-Value Targets

### cclambda and ccresponse (Higher Priority)

Based on the original CC_DIIS_CONSOLIDATION_ANALYSIS.md:

**cclambda DIIS**:
- File: `cc/cclambda/diis.cc` (likely ~400 lines)
- Uses: **DPD buffers** (like ccenergy)
- Migration: **High feasibility** (same as ccenergy)
- Reduction: **~300 lines** (75% estimate)

**ccresponse DIIS**:
- File: `cc/ccresponse/diis.cc` (likely ~300 lines)
- Uses: **DPD buffers** (like ccenergy)
- Migration: **High feasibility** (same as ccenergy)
- Reduction: **~220 lines** (73% estimate)

**Combined Potential**: ~520 lines saved with proven approach

vs fnocc: ~150 lines saved with uncertain approach

---

## Conclusion

The fnocc module's DIIS implementation, while duplicative in concept, is **NOT a priority candidate** for libdiis migration due to:

1. **Data structure mismatch** (raw arrays vs Matrix/DPD objects)
2. **Performance concerns** (copying overhead for large systems)
3. **Lower code reduction payoff** (~150 vs ~720 lines for ccenergy)
4. **Higher-value alternatives** (cclambda, ccresponse with DPD buffers)

**Recommendation**: **Defer fnocc migration**. Focus instead on:
1. ✅ **ccenergy** - COMPLETE (720 lines saved, validated)
2. 🎯 **cclambda** - Next target (~300 lines potential)
3. 🎯 **ccresponse** - Follow-up (~220 lines potential)
4. ⏸️ **fnocc** - Defer until higher priorities complete

This strategy maximizes code reduction with minimal risk and proven methodology.

---

## Technical Notes

### If Migration Were Pursued

**Estimated Effort**: 2-3 days
- 1 day: Implement Matrix wrapping/copying infrastructure
- 0.5 day: Modify DIIS calls
- 0.5 day: Testing and validation
- 0.5-1 day: Performance benchmarking

**Risk Level**: MEDIUM-HIGH
- Primary risk: Performance degradation for large systems
- Secondary risk: Increased code complexity from array↔Matrix conversion

**Testing Requirements**:
- Small molecule tests (o<10, v<50): Verify correctness
- Medium molecule tests (o~20, v~100): Check performance
- Large molecule tests (o>30, v>200): Benchmark overhead

---

**Analysis Complete**: 2025-11-18
**Next Action**: Document findings and focus on cclambda/ccresponse instead
