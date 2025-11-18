# halftrans Migration Assessment: ccenergy and cclambda

## Executive Summary

**Recommendation: DO NOT MIGRATE**

The `halftrans()` functions in ccenergy and cclambda serve a fundamentally different purpose than libtrans and **cannot be replaced** with libtrans functionality. These functions are essential for AO-basis algorithms and should be retained.

However, **consolidation between ccenergy and cclambda** is recommended to eliminate duplication.

---

## Analysis

### 1. What halftrans Does

**Purpose**: Transform only the last two indices of a 4-index tensor between MO and AO bases

**Signature** (ccenergy version):
```cpp
void halftrans(dpdbuf4 *Buf1,      // MO buffer
               int dpdnum1,         // DPD instance for Buf1
               dpdbuf4 *Buf2,       // SO/AO buffer
               int dpdnum2,         // DPD instance for Buf2
               double ***C1,        // Left transformation matrix
               double ***C2,        // Right transformation matrix
               int nirreps,
               int **mo_row,        // MO indexing array
               int **so_row,        // SO indexing array
               int *mospi_left,     // MO dims for left index
               int *mospi_right,    // MO dims for right index
               int *sospi,          // SO dims
               int type,            // 0=MO→SO, 1=SO→MO
               double alpha,        // Scaling factor for source
               double beta)         // Scaling factor for target
```

**Key functionality**:
- Transforms (ij|ab) ↔ (ij|pq) where indices 1,2 stay in MO, indices 3,4 transform
- Bidirectional: MO→AO **and** AO→MO
- Supports alpha/beta scaling for accumulation
- Works with DPD buffers (not IWL)

### 2. Usage Context

#### ccenergy/BT2_AO.cc (Lines 133, 202, etc.)
Purpose: AO-basis algorithm for computing (abcd) contribution to CCSD T2 amplitudes

Workflow:
```
1. Transform T2 amplitudes:     (ij|ab) → (ij|pq)  [MO→AO]
2. Contract with AO integrals:  (pq|rs) × (ij|pq) → (ij|rs)
3. Back-transform result:       (ij|rs) → (ij|cd)  [AO→MO]
```

**When used**: Optional algorithm activated by `AO_BASIS = DISK` or `DIRECT`
- Alternative to standard MO-basis algorithm
- Can be more efficient for large systems
- Avoids storing full (ab|cd) integral list in MO basis

#### cclambda/BL2_AO.cc (Lines 104, 171, etc.)
Purpose: AO-basis algorithm for Lambda equations

Workflow: Similar to BT2_AO
```
1. Transform Lambda:      (ij|ab) → (ij|pq)  [MO→AO]
2. Contract with AO ints: (pq|rs) × (ij|pq) → (ij|rs)
3. Back-transform:        (ij|rs) → (ij|cd)  [AO→MO]
```

### 3. Comparison with libtrans

#### libtrans Half-Transformations

**Purpose**: Break up full 4-index AO→MO transform for memory efficiency

**Process**:
```
Step 1: (pq|rs) → (ij|rs)  [First half:  indices 1,2]
Step 2: (ij|rs) → (ij|ab)  [Second half: indices 3,4]
```

**Direction**: AO → MO only (one direction)

**Methods**:
- `transform_tei_first_half(s1, s2)` - Transform indices 1,2
- `transform_tei_second_half(s1, s2, s3, s4)` - Transform indices 3,4
- `transform_tei(...)` with HalfTrans options

#### Key Differences

| Feature | halftrans (cc modules) | libtrans |
|---------|----------------------|----------|
| **Starting point** | MO integrals/amplitudes | AO integrals |
| **Ending point** | Mixed MO/AO format | Fully MO integrals |
| **Direction** | Bidirectional (MO↔AO) | Unidirectional (AO→MO) |
| **Which indices** | Last two (3,4) only | Both pairs (1,2 then 3,4) |
| **Purpose** | Create mixed MO/AO quantities | Memory-efficient full transform |
| **Use case** | AO-basis algorithms | Standard integral transforms |
| **Output format** | (ij\|pq) mixed MO/AO | (ij\|ab) fully MO |

### 4. Why libtrans Cannot Replace halftrans

1. **Different transformation goal**
   - halftrans: (MO|MO) → (MO|AO) for algorithm purposes
   - libtrans: (AO|AO) → (MO|MO) for integral generation

2. **Missing bidirectionality**
   - halftrans needs MO→AO **and** AO→MO (backtransform)
   - libtrans only does AO→MO (forward transform)

3. **No mixed MO/AO support**
   - libtrans is designed for fully-transformed integrals
   - Does not support (ij|pq) output format where ij≠MO, pq=AO

4. **Different input type**
   - halftrans operates on existing MO quantities (amplitudes, integrals)
   - libtrans operates on AO integrals from disk

5. **No public API for partial transformation**
   - libtrans `trans_one()` is protected, not part of public API
   - Would need significant refactoring to expose this capability

### 5. Code Duplication Analysis

#### Between ccenergy and cclambda

**ccenergy/halftrans.cc**: 142 lines
```cpp
void CCEnergyWavefunction::halftrans(
    dpdbuf4 *Buf1, int dpdnum1,
    dpdbuf4 *Buf2, int dpdnum2,
    double ***C1, double ***C2,  // ← Separate matrices
    int nirreps, int **mo_row, int **so_row,
    int *mospi_left, int *mospi_right,  // ← Different left/right dims
    int *sospi, int type,
    double alpha, double beta)
```

**cclambda/halftrans.cc**: 138 lines
```cpp
void halftrans(
    dpdbuf4 *Buf1, int dpdnum1,
    dpdbuf4 *Buf2, int dpdnum2,
    double ***C,  // ← Single matrix for both
    int nirreps, int **mo_row, int **so_row,
    int *mospi, int *sospi,  // ← Same dims for both
    int type,
    double alpha, double beta)
```

**Differences**:
- ccenergy version is more general: separate C1/C2 matrices, separate left/right dimensions
- cclambda version is specialized: single C matrix, same dimensions

**Core algorithm**: Identical BLAS operations (C_DGEMM calls)

#### Recommendation for Consolidation

✅ **Consolidate into a shared utility**

Location: `psi4/src/psi4/cc/ccwave.h` (CCWavefunction base class)

Rationale:
1. Both modules inherit from CCWavefunction
2. Algorithm is identical, only signature differs
3. ccenergy's more general signature can handle both use cases
4. Would eliminate ~138 lines of duplicated code

---

## Performance Assessment

### Current Implementation

**Characteristics**:
- Direct BLAS (DGEMM) calls for matrix multiplication
- Irrep-blocked for symmetry
- In-place operations on DPD buffers
- Memory efficient: transforms one irrep block at a time

**Performance**: Already optimized
- Uses Level-3 BLAS (DGEMM) - optimal for matrix operations
- Minimal temporary storage (one X matrix per irrep block)
- No unnecessary data movement

### If Migrated to libtrans (Hypothetical)

**Would require**:
1. Adding bidirectional transform capability to libtrans
2. Supporting mixed MO/AO output formats
3. Extending public API to expose partial transforms
4. Additional abstraction layers

**Performance impact**: Likely **negative**
- Additional function call overhead
- More abstraction → less control over memory layout
- Current implementation is already BLAS-optimized
- No clear performance benefit from migration

### Conclusion

Current implementation is already near-optimal. Migration would add complexity without performance benefit.

---

## Usage Frequency

### When is halftrans Used?

**ccenergy/BT2_AO**: Only when user sets `AO_BASIS = DISK` or `DIRECT`
- Not default behavior
- Optional algorithm choice
- Typically for large systems where MO integral storage is prohibitive

**cclambda/BL2_AO**: Similar - optional AO-basis algorithm for Lambda equations

**Frequency**: Relatively rare
- Default algorithms use MO-basis (BT2, not BT2_AO)
- AO-basis algorithms are specialized use cases
- Most users never trigger this code path

---

## Alternative: Consolidation Proposal

### Current Situation

```
ccenergy/halftrans.cc    (142 lines)
cclambda/halftrans.cc    (138 lines)
────────────────────────
Total: 280 lines of nearly identical code
```

### Proposed Solution

**Move to**: `psi4/src/psi4/cc/ccwave.h` + `ccwave.cc`

**Signature** (use ccenergy's more general version):
```cpp
class CCWavefunction {
protected:
    void halftrans(dpdbuf4 *Buf1, int dpdnum1,
                   dpdbuf4 *Buf2, int dpdnum2,
                   double ***C1, double ***C2,
                   int nirreps, int **mo_row, int **so_row,
                   int *mospi_left, int *mospi_right,
                   int *sospi, int type,
                   double alpha, double beta);
};
```

**Usage**:
```cpp
// ccenergy - unchanged
halftrans(&tau, 0, &tau1_AO, 1, C, C, nirreps, ...);

// cclambda - pass C twice
halftrans(&tau, 0, &tau1_AO, 1, C, C, nirreps, ...);
```

**Benefits**:
- Eliminates 138 lines of duplication
- Centralizes maintenance
- More general implementation (ccenergy version) available to both
- No performance impact
- Minimal code changes to call sites

---

## Detailed Recommendations

### 1. DO NOT Migrate to libtrans ❌

**Reasons**:
- Fundamentally different purpose
- Missing required functionality (bidirectional, mixed MO/AO)
- Would require extensive libtrans modifications
- No performance benefit
- Increased complexity for no gain
- These are specialized algorithms, not general transforms

### 2. DO Consolidate Between ccenergy and cclambda ✅

**Action Items**:

1. **Move ccenergy's halftrans to CCWavefunction base class**
   - Location: `psi4/src/psi4/cc/ccwave.h` (declaration)
   - Location: `psi4/src/psi4/cc/ccwave.cc` (implementation)
   - Make it protected member function

2. **Update ccenergy**
   - Remove `ccenergy/halftrans.cc` from sources
   - Update CMakeLists.txt
   - Change calls from `this->halftrans(...)` (should already work if base class member)

3. **Update cclambda**
   - Remove `cclambda/halftrans.cc` from sources
   - Update CMakeLists.txt
   - Change calls to pass C twice: `halftrans(..., C, C, ...)` instead of `halftrans(..., C, ...)`
   - Update signature to match base class

4. **Add documentation**
   - Document purpose in ccwave.h
   - Explain when this is used (AO-basis algorithms)
   - Note difference from libtrans half-transforms

### 3. Document the Distinction 📝

**Add to libtrans/README.txt**:

```
================================================================================
WHEN NOT TO USE LIBTRANS (CONTINUED)
================================================================================

Partial Back-Transformations (MO → AO):
----------------------------------------
Some specialized algorithms (e.g., AO-basis CCSD(T), AO-basis Lambda) need to:
1. Start with MO quantities (T2 amplitudes, Lambda)
2. Partially transform to mixed MO/AO format: (ij|ab) → (ij|pq)
3. Contract with AO integrals
4. Back-transform to MO: (ij|pq) → (ij|cd)

This workflow is fundamentally different from libtrans, which transforms
AO integrals → MO integrals (forward only).

See instead:
• psi4/src/psi4/cc/ccwave.h - CCWavefunction::halftrans()
  Transforms last two indices of dpdbuf4 between MO and AO bases

These are NOT redundant with libtrans - they serve different algorithmic needs.
```

### 4. Add Comment to halftrans Implementation 📝

```cpp
/**
 * halftrans(): Transform the last two indices of a dpdbuf4 between MO and SO bases
 *
 * This function is used in AO-basis algorithms for CCSD(T) and Lambda equations.
 * It enables workflows like:
 *   (ij|ab) MO → (ij|pq) AO → contract with (pq|rs) → (ij|rs) AO → (ij|cd) MO
 *
 * This is DIFFERENT from libtrans half-transformations:
 * - libtrans: (AO|AO) → (MO|AO) → (MO|MO) [full transform in two steps]
 * - halftrans: (MO|MO) ↔ (MO|AO) [partial transform for algorithm]
 *
 * Do NOT attempt to replace this with libtrans. The use cases are distinct.
 *
 * @param Buf1 MO-basis dpdbuf4 (already initialized)
 * @param dpdnum1 DPD instance number for Buf1
 * @param Buf2 SO-basis dpdbuf4 (already initialized)
 * @param dpdnum2 DPD instance number for Buf2
 * @param C1 Left transformation matrix (SO x MO, symmetry blocked)
 * @param C2 Right transformation matrix (SO x MO, symmetry blocked)
 * @param nirreps Number of irreps
 * @param mo_row MO index offsets for each (h, Gc) pair
 * @param so_row SO index offsets for each (h, Gc) pair
 * @param mospi_left Number of MOs per irrep for left index
 * @param mospi_right Number of MOs per irrep for right index
 * @param sospi Number of SOs per irrep
 * @param type 0 = MO → SO; 1 = SO → MO
 * @param alpha Scaling factor for source buffer
 * @param beta Scaling factor for target buffer (allows accumulation)
 */
```

---

## Testing Plan

If consolidation is implemented:

1. **Regression tests**
   - Run all CCSD(T) tests with `AO_BASIS = DISK`
   - Run all CC gradient tests (uses Lambda/BL2_AO)
   - Verify energies unchanged

2. **Performance tests**
   - Benchmark BT2_AO before/after consolidation
   - Ensure no performance regression
   - Should be identical (same BLAS operations)

3. **Build tests**
   - Verify ccenergy/halftrans.cc removed from build
   - Verify cclambda/halftrans.cc removed from build
   - Verify ccwave.cc added if new file created

---

## Implementation Effort

### DO NOT MIGRATE to libtrans
**Effort**: N/A (not recommended)

### DO CONSOLIDATE ccenergy/cclambda
**Effort**: Low (~2-4 hours)

**Steps**:
1. Copy ccenergy/halftrans.cc implementation to ccwave.cc (30 min)
2. Add declaration to ccwave.h (5 min)
3. Update ccenergy CMakeLists.txt (5 min)
4. Update cclambda call sites to pass C twice (30 min)
5. Update cclambda CMakeLists.txt (5 min)
6. Test builds (15 min)
7. Run regression tests (1-2 hours)
8. Add documentation/comments (30 min)

---

## Conclusion

The `halftrans()` functions in ccenergy and cclambda are **NOT redundant** with libtrans. They serve a distinct algorithmic purpose (partial back-transformation for AO-basis methods) that libtrans is not designed to handle.

**Recommended actions**:
1. ❌ **Do NOT** attempt migration to libtrans
2. ✅ **Do** consolidate ccenergy and cclambda versions into CCWavefunction base class
3. ✅ **Do** add documentation explaining the distinction from libtrans
4. ✅ **Do** add comments to prevent future confusion

**Outcome**:
- Eliminates duplication between ccenergy and cclambda (~138 lines)
- Preserves specialized functionality needed for AO-basis algorithms
- Clarifies relationship to libtrans
- Maintains optimal performance

---

## Appendix: Algorithm Comparison Table

| Aspect | libtrans | ccenergy/cclambda halftrans |
|--------|----------|----------------------------|
| **Input** | AO integrals from disk | MO amplitudes/integrals in memory |
| **Output** | Fully MO integrals | Mixed MO/AO quantities |
| **Transform** | (AO\|AO) → (MO\|MO) | (MO\|MO) ↔ (MO\|AO) |
| **Direction** | Forward only (AO→MO) | Bidirectional (MO↔AO) |
| **Indices** | All four (in two steps) | Last two only |
| **Purpose** | Generate MO integrals | Enable AO-basis algorithms |
| **Methods** | MP2, CCSD, CCSD(T) standard | CCSD(T) AO-basis optional |
| **Frequency** | Every correlation calculation | Rare (optional algorithm) |
| **Can replace?** | No - different use cases | No - different use cases |

---

**Document Prepared**: Analysis of halftrans migration feasibility
**Recommendation**: Consolidate between modules, do NOT migrate to libtrans
**Status**: Ready for implementation (consolidation only)
