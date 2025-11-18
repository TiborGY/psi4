# Critical Code Duplications - Quick Reference

This document lists the most critical code duplications requiring immediate attention.

## 1. BLAS/LAPACK Headers (IDENTICAL - 4 files)

**Action Required:** Consolidate to single header

```
psi4/src/psi4/fnocc/blas.h (lines 29-165)
psi4/src/psi4/libqt/blas_intfc23_mangle.h (lines 29-164)
psi4/src/psi4/psimrcc/algebra_interface_mangle.h (lines 29-117)
psi4/src/psi4/mcscf/algebra_interface_mangle.h (lines 29-117)
```

**Effort:** Low | **Risk:** Low | **Impact:** High

---

## 2. MOInfo Structures (IDENTICAL - 7 files)

**Action Required:** Use existing `libmoinfo/moinfo.h`

```
psi4/src/psi4/cc/ccdensity/MOInfo.h
psi4/src/psi4/cc/ccenergy/MOInfo.h
psi4/src/psi4/cc/cceom/MOInfo.h
psi4/src/psi4/cc/cclambda/MOInfo.h
psi4/src/psi4/cc/cchbar/MOInfo.h
psi4/src/psi4/cc/ccresponse/MOInfo.h
psi4/src/psi4/cc/cctriples/MOInfo.h

✓ ALREADY EXISTS: psi4/src/psi4/libmoinfo/moinfo.h (UNUSED!)
```

**Effort:** Medium | **Risk:** Low | **Impact:** High

---

## 3. DIIS Implementations (8+ files)

**Action Required:** Mandate use of existing `libdiis/DIISManager`

```
✓ GENERAL IMPLEMENTATION EXISTS: psi4/src/psi4/libdiis/diismanager.{cc,h}

DUPLICATES TO REMOVE:
psi4/src/psi4/cc/ccenergy/diis_RHF.cc (9,860 bytes)
psi4/src/psi4/cc/ccenergy/diis_ROHF.cc (14,863 bytes)
psi4/src/psi4/cc/ccenergy/diis_UHF.cc (15,241 bytes)
psi4/src/psi4/cc/cclambda/diis.cc (38,787 bytes)
psi4/src/psi4/cc/ccresponse/diis.cc (11,291 bytes)
psi4/src/psi4/dfocc/diis.cc
psi4/src/psi4/fnocc/diis.cc
psi4/src/psi4/psimrcc/blas_diis.cc
psi4/src/psi4/mcscf/scf_diis.cc
```

**Total Duplication:** ~1,500-2,000 lines
**Effort:** Medium | **Risk:** Medium | **Impact:** High

---

## 4. CPHF/CPKS Solvers (4 complete implementations)

**Action Required:** Consolidate to single solver with specialized interfaces

```
psi4/src/psi4/libfock/apps.h:141 - class RCPHF
psi4/src/psi4/libfock/hamiltonian.h:142 - class CPHFRHamiltonian
psi4/src/psi4/fisapt/fisapt.h:213-267 - class CPHF_FISAPT
psi4/src/psi4/libsapt_solver/usapt0.h:338 - class CPKS_USAPT0
```

**Effort:** High | **Risk:** Medium | **Impact:** High

---

## 5. MP2 Implementations (9 files)

**Action Required:** Create MP2Base class with specialized implementations

```
psi4/src/psi4/cc/ccenergy/mp2_energy.cc
psi4/src/psi4/dfocc/cc_energy.cc (lines 69, 174)
psi4/src/psi4/occ/cc_energy.cc (line 77)
psi4/src/psi4/dfmp2/mp2.cc (line 187)
psi4/src/psi4/fnocc/mp2.cc
psi4/src/psi4/f12/mp2.cc
psi4/src/psi4/dlpno/mp2.cc
psi4/src/psi4/dct/dct_mp2_RHF.cc
psi4/src/psi4/dct/dct_mp2_UHF.cc
```

**Effort:** High | **Risk:** Medium | **Impact:** High

---

## 6. Triples (T) Implementations (9+ files)

**Action Required:** Unify around single (T) kernel with abstracted integral backend

```
psi4/src/psi4/cc/cctriples/triples.cc
psi4/src/psi4/fnocc/triples.cc (lines 387-390)
psi4/src/psi4/fnocc/lowmemory_triples.cc (lines 532-535)
psi4/src/psi4/dfocc/ccsd_triples.cc
psi4/src/psi4/dfocc/uccsd_triples_hm.cc (lines 271, 491)
psi4/src/psi4/dfocc/uccsd_triples_grad_hm.cc
psi4/src/psi4/dfocc/uccsdl_triples_hm.cc
psi4/src/psi4/dct/dct_triples.cc (lines 225-228)
psi4/src/psi4/psimrcc/mrcc_pert_triples.cc
psi4/src/psi4/psimrcc/mrccsd_t_compute*.cc

Plus spin-case variants: T3_AAA.cc, T3_AAB.cc, T3_BBB.cc, etc.
```

**Estimated Duplication:** ~3,000-5,000 lines
**Effort:** Very High | **Risk:** Medium-High | **Impact:** Very High

---

## 7. Tensor/Array Classes (3 nearly identical implementations)

**Action Required:** Create unified template-based tensor library

```
psi4/src/psi4/dfocc/tensors.cc (Tensor1d, Tensor2d, Tensor3d)
psi4/src/psi4/dfocc/arrays.cc (Array1d, Array2d, Array3d)
psi4/src/psi4/occ/arrays.cc (Array1d, Array2d, Array3d)
```

**What's duplicated:**
- Memory management (memalloc, release, init, zero)
- Vector operations (set, add, subtract, get, dot, rms)
- Matrix operations (gemv, to_lower_triangle, vector_dot)
- Printing functions

**Module sizes:**
- DFOCC: 123 files, 108,000 lines total
- OCC: 54 files, similar structure
- Estimated tensor duplication: 5,000-10,000 lines

**Effort:** Very High | **Risk:** Medium | **Impact:** Very High

---

## 8. get_moinfo() Functions (9 identical implementations)

**Action Required:** Create shared get_moinfo() in libccbase/

```
psi4/src/psi4/cc/ccenergy/get_moinfo.cc (lines 62-80)
psi4/src/psi4/cc/cclambda/get_moinfo.cc (lines 59-80)
psi4/src/psi4/cc/ccdensity/get_moinfo.cc
psi4/src/psi4/cc/cceom/get_moinfo.cc
psi4/src/psi4/cc/cchbar/get_moinfo.cc
psi4/src/psi4/cc/ccresponse/get_moinfo.cc
psi4/src/psi4/cc/cctriples/get_moinfo.cc
psi4/src/psi4/dfocc/get_moinfo.cc
psi4/src/psi4/occ/get_moinfo.cc
```

**Effort:** Medium | **Risk:** Low | **Impact:** Medium

---

## 9. Integral Sorting (3 major implementations)

**Action Required:** Create unified IntegralSorter class

```
psi4/src/psi4/dct/dct_integrals_RHF.cc (lines 46-243)
psi4/src/psi4/dct/dct_integrals_UHF.cc (lines 71-460)
psi4/src/psi4/fnocc/sortintegrals.cc (2,198 lines!)
psi4/src/psi4/psimrcc/sort_out_of_core.cc
psi4/src/psi4/psimrcc/sort_mrpt2.cc
```

**Duplicated sorting for:** OOOO, OOVV, OVOV, OVVV, VVVV, OOOV integrals
**Estimated Duplication:** ~2,000-3,000 lines
**Effort:** High | **Risk:** Medium | **Impact:** High

---

## 10. Davidson Solvers (3 implementations)

**Action Required:** Create unified DavidsonSolver class in libqt

```
psi4/src/psi4/detci/diag_h.cc (lines 143-278) - "Davidson/Liu"
psi4/src/psi4/detci/mitrush_iter.cc - "Mitrushenkov's Olsen-modified Davidson"
psi4/src/psi4/cc/cceom/diag.cc (line 34) - "Davidson-Liu for EOM"
```

**Effort:** High | **Risk:** Medium | **Impact:** Medium

---

## 11. MO Integral Transformation (4 implementations)

**Action Required:** Standardize on libtrans::IntegralTransform

```
✓ MAIN IMPLEMENTATION: psi4/src/psi4/libtrans/integraltransform.h

DUPLICATES:
psi4/src/psi4/dct/half_transform.cc
psi4/src/psi4/psimrcc/transform.h
psi4/src/psi4/fnocc/df_t1_transformation.cc

Referenced in 82 files total
```

**Effort:** High | **Risk:** Medium | **Impact:** High

---

## 12. Orbital Orthogonalization (6 implementations)

**Action Required:** Consolidate to libmints/orthog methods

```
✓ MAIN IMPLEMENTATION: psi4/src/psi4/libmints/orthog.{cc,h}

DUPLICATES:
psi4/src/psi4/libqt/schmidt.cc (lines 46-73)
psi4/src/psi4/mcscf/scf_S_inverse_sqrt.cc (lines 35-62)
psi4/src/psi4/libmints/matrix.cc:1584 - schmidt_orthog_columns()
psi4/src/psi4/libmints/matrix.cc:1951 - canonical_orthogonalization()
psi4/src/psi4/detci/civect.cc (lines 1857-1977) - Custom Gram-Schmidt
psi4/src/psi4/libscf_solver/rohf.cc:420 - prepare_canonical_orthogonalization()
```

**Effort:** Medium | **Risk:** Low | **Impact:** Medium

---

## 13. DCT Module RHF/UHF Split (22 paired files)

**Action Required:** Template or use policy-based design for spin cases

```
PATTERN: Nearly every DCT function has _RHF.cc and _UHF.cc versions:

psi4/src/psi4/dct/dct_compute_{RHF,UHF}.cc
psi4/src/psi4/dct/dct_scf_{RHF,UHF}.cc
psi4/src/psi4/dct/dct_oo_{RHF,UHF}.cc
psi4/src/psi4/dct/dct_integrals_{RHF,UHF}.cc
psi4/src/psi4/dct/dct_intermediates_{RHF,UHF}.cc
psi4/src/psi4/dct/dct_energy_{RHF,UHF}.cc
psi4/src/psi4/dct/dct_gradient_{RHF,UHF}.cc
psi4/src/psi4/dct/dct_lambda_{RHF,UHF}.cc
psi4/src/psi4/dct/dct_tau_{RHF,UHF}.cc
psi4/src/psi4/dct/dct_density_{RHF,UHF}.cc
psi4/src/psi4/dct/dct_mp2_{RHF,UHF}.cc
```

**Total:** 22 paired files (11 pairs)
**Effort:** Very High | **Risk:** High | **Impact:** High

---

## 14. Index/Utility Macros (Duplicate defines.h)

**Action Required:** Consolidate macro definitions

```
psi4/src/psi4/occ/defines.h
psi4/src/psi4/dfocc/defines.h

Both define:
- index2(i,j), index4(i,j,k,l)
- MIN0(a,b), MAX0(a,b)
- DIIS_MIN_DET 1.0E-16
- DIVERGE 1.0E+3

dfocc also has: index3, idx2, idx3, idx4, idx_asym
```

**Effort:** Low | **Risk:** Low | **Impact:** Low

---

## 15. Matrix Operations (Scattered duplications)

**Action Required:** Consolidate to common utilities

### to_lower_triangle():
```
psi4/src/psi4/libmints/matrix.cc:758
psi4/src/psi4/dfocc/tensors.cc:2441
psi4/src/psi4/dfocc/arrays.cc:757
```

### vector_dot():
```
psi4/src/psi4/libmints/matrix.cc (lines 1727, 1745, 2814)
psi4/src/psi4/libmints/vector.cc (~line 350)
psi4/src/psi4/dfocc/tensors.cc (lines 1847, 1856)
psi4/src/psi4/dfocc/arrays.cc
psi4/src/psi4/occ/arrays.cc
```

### rms():
```
psi4/src/psi4/dfocc/tensors.cc (lines 194, 3478)
psi4/src/psi4/occ/arrays.cc:160
psi4/src/psi4/dfocc/arrays.cc
psi4/src/psi4/libmints/matrix.h:740
psi4/src/psi4/libmints/vector.h:355
```

**Effort:** Medium | **Risk:** Low | **Impact:** Medium

---

## Quick Priority Summary

### Immediate (Low-Hanging Fruit):
1. **BLAS headers** - 4 files → 1
2. **MOInfo structures** - Use existing libmoinfo
3. **Index macros** - Consolidate defines.h
4. **String utilities** - Use libpsi4util consistently

### Short-term (High Impact, Manageable Risk):
5. **DIIS** - Mandate DIISManager usage
6. **get_moinfo()** - Create shared implementation
7. **Orthogonalization** - Use libmints/orthog
8. **Matrix operations** - Consolidate utilities

### Medium-term (Significant Refactoring):
9. **CPHF/CPKS** - Unify solvers
10. **MP2** - Create base class
11. **Integral sorting** - IntegralSorter class
12. **Davidson** - DavidsonSolver class
13. **MO transformation** - Standardize on libtrans

### Long-term (Major Projects):
14. **Triples** - Unify (T) implementations
15. **Tensor library** - Unify DFOCC/OCC
16. **DCT spin cases** - Template refactoring

---

## Estimated Impact

| Priority | Files Affected | Lines Duplicated | Effort | Risk |
|----------|---------------|------------------|--------|------|
| Immediate | 15-20 | 500-1,000 | Low | Low |
| Short-term | 30-50 | 3,000-5,000 | Medium | Low-Medium |
| Medium-term | 100-150 | 8,000-15,000 | High | Medium |
| Long-term | 200+ | 15,000-30,000 | Very High | Medium-High |

**Total Addressable Duplication:** 20,000-40,000 lines
**Potential Code Reduction:** 15,000-30,000 lines (50-75%)
