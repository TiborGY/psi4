# Integral Transformation Redundancy Analysis

## Executive Summary

This document analyzes the redundancy in integral transformation code across the Psi4 codebase. Multiple modules implement similar functionality for AO→MO transformations, frozen core handling, and integral sorting. While `libtrans` provides a canonical framework, several modules reimplement portions of this functionality.

## Key Locations

### Primary Libraries

1. **psi4/src/psi4/libtrans/** - General AO→MO transformation library
   - Main class: `IntegralTransform` (integraltransform.h)
   - Comprehensive transformation framework with DPD and IWL support
   - Handles frozen core operator and energy calculation
   - Provides half-transformation infrastructure

2. **psi4/src/psi4/cctransort/** - CC-specific integral sorting
   - Uses `IntegralTransform` from libtrans
   - Adds specialized sorting patterns for CC methods
   - Sorts integrals into A, B, C, D, E, F buffers
   - Used by traditional CC methods via `RUN_CCTRANSORT` option

3. **psi4/src/psi4/psimrcc/** - MRCC transformation
   - Class: `CCTransform` (transform.h/transform.cc)
   - Independent implementation, does NOT use libtrans
   - Custom integral storage and transformation routines

### Method-Specific Implementations

4. **psi4/src/psi4/cc/ccenergy/** - CC energy module
   - Has custom `halftrans()` function (halftrans.cc)
   - Also uses `IntegralTransform` for density-fitted integrals (form_df_ints.cc)
   - Mixing of libtrans and custom implementations

5. **psi4/src/psi4/cc/cclambda/** - CC lambda module
   - Has own `halftrans()` implementation
   - Similar to ccenergy's version

6. **psi4/src/psi4/dct/** - DCT module
   - Has `half_transform()` function (half_transform.cc)
   - Based on ccenergy's halftrans but adapted for DCT
   - Does not use libtrans half-transformation routines

7. **psi4/src/psi4/occ/** - Orbital-optimized CC
   - Uses `IntegralTransform` from libtrans properly (trans_ints_rhf.cc, trans_ints_uhf.cc)
   - Good example of canonical usage

8. **psi4/src/psi4/dfocc/** - Density-fitted OCC
   - Does not use libtrans (density fitting uses different approach)
   - Not redundant - different methodology

## Detailed Redundancy Analysis

### 1. Frozen Core Operator Construction

**Canonical Implementation: libtrans/integraltransform_sort_so_tei.cc**

Location: Lines 191-214, 350-387

```cpp
// Constructs frozen core density matrix
for (int h = 0; h < nirreps_; ++h) {
    for (int p = soOffset; p < soOffset + sopi_[h]; ++p) {
        for (int q = 0; q <= p; ++q) {
            for (int i = 0; i < frzcpi_[h]; ++i)
                aFzcD[pq] += pCa[p][i] * pCa[q][i];
        }
    }
}

// Constructs frozen core operator: H_core + J - K from frozen orbitals
```

**Duplicate Implementations:**

- **None found** - This appears to be unique to libtrans
- Other modules rely on libtrans to provide the frozen core operator
- From README.txt: "some codes (mrcc, detci, and especially cc) request that tasks involving frozen core orbitals be either done by libtrans or converted into effective quantities free of frozen core orbitals"

**Status:** ✅ Canonical implementation used

### 2. Frozen Core Energy Calculation

**Canonical Implementation: libtrans/integraltransform_sort_so_tei.cc**

Location: Lines 354-375

```cpp
frozen_core_energy_ = 0.0;
for (int p = 0; p < nso_; p++) {
    for (int q = 0; q <= p; q++, pq++) {
        double prefact = p == q ? 1.0 : 2.0;
        frozen_core_energy_ += prefact * aFzcD[pq] * (aoH[pq] + aFzcOp[pq]);
    }
}
```

Accessor: `double get_frozen_core_energy() const` (integraltransform.h:188)

**Duplicate Implementations:**

- **psimrcc**: May have its own calculation (needs verification in sort.h/sort_*.cc)
- Most other modules use `ints->get_frozen_core_energy()` from libtrans

**Status:** ⚠️ Likely some duplication in psimrcc

### 3. Half-Transformations

**Canonical Implementation: libtrans**

Files:
- integraltransform_tei_1st_half.cc
- integraltransform_tei_2nd_half.cc
- integraltransform_tei.cc

Public interface:
```cpp
void transform_tei_first_half(const std::shared_ptr<MOSpace> s1,
                               const std::shared_ptr<MOSpace> s2);
void transform_tei_second_half(const std::shared_ptr<MOSpace> s1,
                                const std::shared_ptr<MOSpace> s2,
                                const std::shared_ptr<MOSpace> s3,
                                const std::shared_ptr<MOSpace> s4);
```

**Duplicate Implementations:**

#### ccenergy/halftrans.cc (Lines 40-60)
```cpp
/* halftrans(): Routine to transform the last two indices of a dpdbuf4
** between the MO and SO bases.
**
** dpdbuf4 *Buf1:    Pointer to the MO dpdbuf4 (already initialized)
** dpdbuf4 *Buf2:    Pointer to the SO dpdbuf4 (already initialized).
** ...
```

Function signature:
```cpp
void CCEnergyWavefunction::halftrans(dpdbuf4 *Buf1, int dpdnum1,
                                      dpdbuf4 *Buf2, int dpdnum2,
                                      double ***C1, double ***C2,
                                      int nirreps, int **mo_row, int **so_row,
                                      int *mospi_left, int *mospi_right,
                                      int *sospi, int type,
                                      double alpha, double beta)
```

#### cclambda/halftrans.cc
Similar implementation to ccenergy (separate copy)

#### dct/half_transform.cc (Lines 36-54)
```cpp
/* half_transform(): Routine to transform the last two indices of a dpdbuf4
 * between the MO and SO bases.
 * Based on code originally written by Daniel Crawford in ccenergy/halftrans.cc
```

Function signature:
```cpp
void DCTSolver::half_transform(dpdbuf4 *SO, dpdbuf4 *MO,
                                SharedMatrix &C1, SharedMatrix &C2,
                                int *mospi_left, int *mospi_right,
                                int **so_row, int **mo_row,
                                bool backwards,
                                double alpha, double beta)
```

**Analysis:**
- ccenergy, cclambda, and dct all have their own halftrans implementations
- These are specialized for transforming the last two indices of dpdbuf4 structures
- libtrans has more general transformation capabilities but may not expose this exact interface
- Code comment in dct acknowledges it's based on ccenergy's version

**Status:** ⚠️ **REDUNDANT** - Three separate implementations of similar functionality

### 4. Integral Sorting and Indexing Schemes

**libtrans approach:**
- Uses DPD (Distributed Packed Dense) format
- Supports IWL (Integrals With Labels) format
- Provides `DPD_ID()` methods for indexing
- Space-based transformations using MOSpace objects

**cctransort approach:**
- Builds on libtrans transformations
- Sorts integrals into specific buffers:
  - A: `<ij|kl>` (occ-occ | occ-occ)
  - B: `<ab|cd>` (vir-vir | vir-vir)
  - C: `<ia|jb>` (occ-vir | occ-vir)
  - D: `<ij|ab>` (occ-occ | vir-vir)
  - E: `<ai|jk>` (vir-occ | occ-occ)
  - F: `<ia|bc>` (occ-vir | vir-vir)

From sort_tei_rhf.cc:
```cpp
global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "kl", ...);
global_dpd_->buf4_sort(&K, PSIF_CC_AINTS, prqs, "ij", "kl", "A <ij|kl>");
```

**psimrcc approach:**
- Custom indexing using CCIndex class
- Uses `std::map<size_t, double> integral_map` for storage
- Block-based memory management
- Does NOT use libtrans infrastructure

From transform.h:
```cpp
CCIndex* oei_so_indexing;
CCIndex* tei_so_indexing;
CCIndex* tei_mo_indexing;
std::map<size_t, double> integral_map;
```

**Status:** 🔶 **MIXED** - cctransort complements libtrans; psimrcc is independent

## Usage Patterns

### Modules Using libtrans (Canonical)

1. **occ** - Exemplary usage
   ```cpp
   ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ,
                       IntegralTransform::HalfTrans::MakeAndKeep);
   ```

2. **ccenergy** (partial) - Uses for DF integrals
   - form_df_ints.cc includes libtrans headers

3. **detci** - References in comments suggest usage

### Modules NOT Using libtrans

1. **psimrcc** - Complete independent implementation
   - Rationale: May need different memory management or blocking
   - Custom CCTransform class

2. **dfocc** - Uses density-fitting approach
   - Rationale: Different integral technology

### Modules with Mixed Approach

1. **ccenergy/cclambda** - Uses libtrans AND custom halftrans
   - Likely for performance or historical reasons
   - Custom halftrans is more specialized

2. **dct** - Custom half_transform
   - Explicitly acknowledges derivation from ccenergy

## Recommendations

### 1. Document Canonical Usage Pattern ✅ HIGH PRIORITY

**Action:** Create comprehensive documentation for libtrans usage

**Location:** Expand `psi4/src/psi4/libtrans/README.txt`

**Include:**
- Standard transformation workflow examples
- When to use transform_tei vs transform_tei_first_half/second_half
- How to properly use frozen core functionality
- Space specification best practices
- Performance considerations (MakeAndKeep vs ReadAndNuke)

**Example to add:**
```cpp
// Standard transformation pattern for correlation methods
auto ints = std::make_shared<IntegralTransform>(
    wfn, spaces,
    IntegralTransform::TransformationType::Restricted,
    IntegralTransform::OutputType::DPDOnly,
    IntegralTransform::MOOrdering::QTOrder,
    IntegralTransform::FrozenOrbitals::OccAndVir
);

ints->update_orbitals();
ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir);
double fzc_energy = ints->get_frozen_core_energy();
```

### 2. Evaluate Custom Half-Transformation Necessity ⚠️ MEDIUM PRIORITY

**Analysis needed:**
- Are ccenergy/cclambda/dct halftrans functions performance-critical?
- Can they be replaced with libtrans equivalents?
- What specific features do they provide that libtrans doesn't?

**Suggested approach:**
1. Profile performance difference between custom halftrans and libtrans
2. If performance difference is negligible (<5%), deprecate custom versions
3. If custom versions are faster, document WHY and when to use them
4. Consider contributing optimizations back to libtrans

**Potential outcome:**
- Keep custom implementations if justified by performance
- Add documentation explaining the trade-offs
- Mark as "performance-optimized alternative to libtrans"

### 3. psimrcc Independence - Document Rationale 📝 LOW PRIORITY

**Current status:** psimrcc has complete independent implementation

**Action:** Add documentation to `psi4/src/psi4/psimrcc/transform.h`

```cpp
/**
 * CCTransform - PSIMRCC-specific integral transformation
 *
 * This class provides integral transformations for PSIMRCC methods.
 * It does NOT use libtrans because:
 * 1. Block-based memory management for large active spaces
 * 2. Custom integral storage optimized for MRCC algorithms
 * 3. Direct integration with PSIMRCC's CCIndex system
 *
 * For standard single-reference methods, prefer libtrans/IntegralTransform.
 */
```

### 4. Create Migration Guide for New Methods 📋 MEDIUM PRIORITY

**Location:** `psi4/docs/development/integral_transformations.md`

**Content:**
```markdown
# Integral Transformation Best Practices

## When to use libtrans/IntegralTransform

Use the standard IntegralTransform class if:
- Implementing single-reference correlation methods
- Need standard frozen core handling
- Want DPD or IWL formatted output
- Transform between well-defined orbital spaces

Example methods: OCC, CCSD, CCSD(T), MP2

## When to write custom transformations

Only implement custom transformations if:
- Using density-fitting (see dfocc)
- Multi-reference with custom blocking (see psimrcc)
- Performance profiling shows >10% improvement possible
- Require fundamentally different integral access patterns

## Hybrid Approaches

Some modules use libtrans + custom code:
- ccenergy: Uses libtrans for main transforms, custom halftrans for specialized operations
- Justify mixed approaches in code comments
```

### 5. Consolidate or Document Frozen Core Energy ⚠️ MEDIUM PRIORITY

**Action:**
1. Verify if psimrcc calculates frozen core energy separately
2. If yes, document why (e.g., "MRCC uses different frozen core definition")
3. If no clear reason, migrate to using libtrans value

**Add to psimrcc documentation:**
```cpp
// Frozen core energy handling:
// Option A: Use libtrans value if compatible
double fzc_e = ints->get_frozen_core_energy();

// Option B: If MRCC-specific calculation needed, document:
// "MRCC frozen core energy differs from libtrans because..."
```

### 6. Performance Benchmarking Study 🔬 LOW PRIORITY

**Objective:** Quantify performance differences

**Tasks:**
1. Benchmark libtrans vs custom halftrans on representative systems
2. Profile memory usage patterns
3. Identify specific bottlenecks
4. Document findings

**Outcome:**
- Data-driven decisions on what custom code to keep
- Potential optimizations to contribute back to libtrans
- Clear performance documentation

## Summary Table

| Module | Uses libtrans | Custom Implementation | Status | Priority |
|--------|--------------|----------------------|--------|----------|
| libtrans | ✅ (Core) | N/A | Canonical | - |
| cctransort | ✅ (Wraps) | Specialized sorting | Complementary | Low |
| psimrcc | ❌ | Complete custom | Independent | Document |
| ccenergy | 🔶 Partial | halftrans() | Mixed | Evaluate |
| cclambda | 🔶 Partial | halftrans() | Redundant | Evaluate |
| dct | ❌ | half_transform() | Derivative | Evaluate |
| occ | ✅ Full | None | Best practice | None |
| dfocc | ❌ | DF-specific | Different tech | None |

## Conclusion

The integral transformation infrastructure in Psi4 shows moderate redundancy:

**✅ Well-designed aspects:**
- libtrans provides comprehensive canonical framework
- Frozen core operator/energy centralized in libtrans
- OCC module demonstrates clean usage pattern
- cctransort properly complements rather than duplicates

**⚠️ Areas for improvement:**
- halftrans/half_transform duplicated across 3 modules (ccenergy, cclambda, dct)
- Lack of clear documentation on when to use libtrans vs custom
- No performance justification documented for custom implementations

**🔶 Acceptable independence:**
- psimrcc has valid reasons for custom implementation (MRCC-specific needs)
- dfocc uses different technology (density fitting)

**Priority actions:**
1. Document canonical libtrans usage patterns
2. Evaluate halftrans necessity with performance data
3. Add rationale documentation to custom implementations
4. Create development guide for new methods

The redundancy is manageable and some custom implementations may be justified, but better documentation is essential to prevent future duplication and guide developers toward appropriate choices.
