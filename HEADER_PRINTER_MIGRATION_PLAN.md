# HeaderPrinter Migration Plan

## Executive Summary

This document outlines the migration strategy for converting the remaining 40+ `print_header()` implementations across the Psi4 codebase to use the new `HeaderPrinter` utility class.

**Status**: 3 of 43 implementations migrated (7%)
**Estimated Total Effort**: 8-12 hours
**Recommended Timeline**: 3-4 phases over 2-3 weeks

## Migration Status

### ✅ Completed (3/43)
- `LinK::print_header()` - libfock/LinK.cc:86
- `MemDFJK::print_header()` - libfock/MemDFJK.cc:113
- `RDFMP2::print_header()` - dfmp2/mp2.cc:796

### 🔄 Remaining (40/43)

## Phase 1: Simple Arrow-Style Headers (High Priority)
**Estimated Effort**: 2-3 hours | **Difficulty**: Easy | **Files**: 15

These implementations follow the simple `==> Title <==` pattern with basic parameters. Low risk, high impact for code reduction.

### libfock - Simple JK Algorithms
| File | Line | Class | Parameters | Notes |
|------|------|-------|------------|-------|
| DirectJK.cc | 109 | DirectJK | 1 param (cutoff) | Similar to LinK |
| DiskJK.cc | 58 | DiskJK | 1 param (cutoff) | Similar to LinK |
| PKJK.cc | 72 | PKJK | 1 param (cutoff) | Similar to LinK |
| DirectDFJ.cc | 79 | DirectDFJ | 1 param (cutoff) | Similar to LinK |
| CDJK.cc | 145 | CDJK | 3 params | Slightly more complex |

### libfock - Specialized K Builders
| File | Line | Class | Parameters | Notes |
|------|------|-------|------------|-------|
| COSK.cc | 314 | COSK | 5 params | Medium complexity |
| snLinK.cc | 420 | snLinK | 2 params | Conditional compilation |

### libfock - Potentials and Solvers
| File | Line | Class | Parameters | Notes |
|------|------|-------|------------|-------|
| v.cc | 728 | VBase | 4 params | Base class |
| v.cc | 1118 | SAP | 2 params | Extends VBase |
| v.cc | 1260 | RV | 0 params | Delegates to VBase |
| v.cc | 2594 | UV | 0 params | Delegates to VBase |
| hamiltonian.cc | 90 | CPHFRHamiltonian | 1 param | Very simple (2 lines) |
| solver.cc | 110 | CGRSolver | 3 params | Simple solver |
| apps.cc | 144 | RCPHF | Multiple | Needs investigation |

### lib3index - Density Fitting Infrastructure
| File | Line | Class | Parameters | Notes |
|------|------|-------|------------|-------|
| dfhelper.cc | 278 | DFHelper | ~10 params | Configuration heavy |
| dftensor.cc | 123 | DFTensor | ~8 params | Configuration heavy |

**Success Criteria**:
- All simple arrow-style implementations migrated
- Code reduction of 50-70% in these files
- No changes to output formatting
- All implementations use `print_if(print_)` pattern

---

## Phase 2: Box-Style Headers with Authors (Medium Priority)
**Estimated Effort**: 3-4 hours | **Difficulty**: Medium | **Files**: 10

These implementations use the box format with author credits, similar to RDFMP2 (already migrated).

### dfmp2 - Remaining MP2 Variants
| File | Line | Class | Similarity to RDFMP2 | Notes |
|------|------|-------|----------------------|-------|
| mp2.cc | 2609 | UDFMP2 | Very high (95%) | Unrestricted variant |
| mp2.cc | 3235 | RODFMP2 | Very high (95%) | Restricted-open variant |
| corr_grad.cc | 108 | DFCorrGrad | High (80%) | Gradient implementation |

### libsapt_solver - SAPT Methods
| File | Line | Class | Pattern | Notes |
|------|------|-------|---------|-------|
| sapt0.cc | 171 | SAPT0 | BOX + authors | Foundation SAPT |
| usapt0.cc | 180 | USAPT0 | BOX + authors | Unrestricted SAPT0 |
| sapt2.cc | 182 | SAPT2 | BOX + authors | Second-order SAPT |
| sapt2p.cc | 143 | SAPT2p | BOX + authors | SAPT2+ variant |
| sapt2p3.cc | 166 | SAPT2p3 | BOX + authors | SAPT2+(3) variant |

### Other Box-Style Headers
| File | Line | Class | Pattern | Notes |
|------|------|-------|---------|-------|
| f12/mp2.cc | 92 | MP2F12 | BOX + authors | F12 correction |
| dlpno/mp2.cc | 304 | DLPNOMP2 | BOX + authors | Local correlation |

**Migration Strategy**:
1. Start with UDFMP2 and RODFMP2 (nearly identical to completed RDFMP2)
2. Extract common author lists into constants for SAPT family
3. Use `HeaderPrinter::BannerStyle::BOX` consistently
4. Maintain exact width (57 for MP2, varies for others)

**Success Criteria**:
- All box-style headers properly centered
- Author credits preserved exactly
- Subtitle formatting maintained
- Thread information included where applicable

---

## Phase 3: Complex Headers with Custom Formatting (Lower Priority)
**Estimated Effort**: 2-3 hours | **Difficulty**: Medium-Hard | **Files**: 8

These implementations have special formatting requirements or complex conditional logic.

### High-Complexity Implementations
| File | Line | Class | Complexity | Reason |
|------|------|-------|------------|--------|
| libscf_solver/hf.cc | 490 | HF | Very High | 47 lines, geometry, symmetry, algorithm info |
| libfock/DiskDFJK.cc | 270 | DiskDFJK | High | Complex conditional parameters |
| libfock/CompositeJK.cc | 201 | CompositeJK | High | Recursive printing of sub-algorithms |

### scfgrad - Gradient Implementations
| File | Line | Class | Complexity | Reason |
|------|------|-------|------------|--------|
| scfgrad/jk_grad.cc | 156 | DFJKGrad | Medium | Gradient-specific parameters |
| scfgrad/jk_grad.cc | 2213 | DirectJKGrad | Medium | Similar to DFJKGrad |

### Specialized Implementations
| File | Line | Class | Complexity | Reason |
|------|------|-------|------------|--------|
| fisapt/fisapt.cc | 96 | FISAPT | Medium | Fragment-specific info |
| fisapt/local2.cc | 97 | IBOLocalizer2 | Low | Simple localizer |
| libmints/thc_eri.cc | 59 | LS_THC_Computer | Medium | THC-specific parameters |

**Migration Strategy**:
1. **HF::print_header()**: Most complex - consider partial migration
   - Main banner → HeaderPrinter
   - Keep custom geometry/symmetry printing as-is
   - Migrate algorithm parameters section
2. **CompositeJK**: Use `print_to()` for sub-algorithm recursion
3. **DiskDFJK**: Handle conditional logic with intermediate variables
4. **Gradient implementations**: May need custom formatting for matrix data

**Success Criteria**:
- Complex conditional logic preserved
- Custom sections (geometry, symmetry) maintained
- No regression in output detail
- Code still readable despite complexity

---

## Phase 4: Utility and Visualization Classes (Lowest Priority)
**Estimated Effort**: 1-2 hours | **Difficulty**: Easy-Medium | **Files**: 7

These are peripheral utilities that could benefit from standardization but are lower priority.

### Orbital Localization
| File | Line | Class | Type | Notes |
|------|------|-------|------|-------|
| libmints/local.cc | 142 | BoysLocalizer | Simple | 3 params |
| libmints/local.cc | 341 | PMLocalizer | Simple | 3 params |

### Property Calculations
| File | Line | Class | Type | Notes |
|------|------|-------|------|-------|
| libmints/oeprop.cc | 783 | OEProp | Medium | Multiple properties |
| libcubeprop/cubeprop.cc | 99 | CubeProperties | Simple | Grid parameters |
| libcubeprop/csg.cc | 202 | CubicScalarGrid | Simple | Grid parameters |

**Migration Strategy**:
- Low priority but easy wins
- Can be done opportunistically during other work
- Good for new contributors as starter tasks

**Success Criteria**:
- Consistent formatting with core algorithms
- Parameter names standardized

---

## Implementation Guidelines

### Before Starting Each Migration

1. **Read the current implementation** completely
2. **Identify the banner style** (ARROW, BOX, NONE, or CUSTOM)
3. **List all parameters** and their types
4. **Note conditional logic** (if statements, compile-time flags)
5. **Check for custom sections** that should remain untouched

### Migration Checklist

- [ ] Add `#include "psi4/libpsi4util/header_printer.h"` to file
- [ ] Determine appropriate banner style
- [ ] Identify all parameters and their printf format strings
- [ ] Handle conditional parameters (if-wrapped)
- [ ] Preserve `print_if(print_)` or similar guards
- [ ] Test output matches original exactly
- [ ] Verify code reduction achieved
- [ ] Update local documentation if needed

### Common Patterns

#### Pattern 1: Simple Replacement
```cpp
// Before (7 lines)
if (print_) {
    outfile->Printf("\n");
    outfile->Printf("  ==> Algorithm <==\n\n");
    outfile->Printf("    Parameter:%11.0E\n", value);
}

// After (3 lines)
HeaderPrinter header("Algorithm");
header.add_parameter("Parameter", value)
      .print_if(print_);
```

#### Pattern 2: Conditional Parameters
```cpp
// Before
if (print_) {
    outfile->Printf("  ==> Algorithm <==\n\n");
    outfile->Printf("    Param1:%11d\n", p1);
    if (do_feature) outfile->Printf("    Param2:%11.3E\n", p2);
    outfile->Printf("    Param3:%11d\n", p3);
}

// After
if (print_) {
    HeaderPrinter header("Algorithm");
    header.add_parameter("Param1", p1);
    if (do_feature) {
        header.add_parameter("Param2", p2, "%11.3E");
    }
    header.add_parameter("Param3", p3)
          .print();
}
```

#### Pattern 3: Delegation
```cpp
// Before
void RV::print_header() const { VBase::print_header(); }
void UV::print_header() const { VBase::print_header(); }

// After - No change needed, these already delegate properly
```

#### Pattern 4: Box with Authors
```cpp
// Before (15+ lines)
outfile->Printf("\t -------------------\n");
outfile->Printf("\t        Title\n");
outfile->Printf("\t     Subtitle\n");
outfile->Printf("\t\n");
outfile->Printf("\t    Author Names\n");
outfile->Printf("\t -------------------\n");

// After (5-7 lines)
HeaderPrinter header("Title", HeaderPrinter::BannerStyle::BOX, 57);
header.subtitle("Subtitle")
      .add_authors({"Author Names"})
      .print();
```

---

## Risk Assessment

### Low Risk (Phases 1, 2, 4)
- **Impact**: Minimal
- **Testing**: Output comparison only
- **Rollback**: Easy (git revert)
- **Examples**: Simple JK algorithms, localizers

### Medium Risk (Phase 3 - most)
- **Impact**: Moderate
- **Testing**: Full regression tests recommended
- **Rollback**: Moderate effort
- **Examples**: DiskDFJK, gradient implementations

### High Risk (Phase 3 - HF only)
- **Impact**: High (affects all SCF calculations)
- **Testing**: Extensive testing required
- **Rollback**: May require careful cherry-picking
- **Examples**: HF::print_header()
- **Mitigation**: Consider partial migration or skip

---

## Testing Strategy

### For Each Migration

1. **Build Test**: Ensure code compiles
   ```bash
   cmake --build . --target psi4 -- -j4
   ```

2. **Output Comparison**: Run test case, compare output
   ```bash
   # Save reference output before migration
   psi4 test_case.in -o ref.out

   # After migration
   psi4 test_case.in -o new.out

   # Compare (should be identical)
   diff ref.out new.out
   ```

3. **Regression Tests**: Run relevant test suite
   ```bash
   ctest -R <module_name> -V
   ```

### Acceptance Criteria

- ✅ Code compiles without warnings
- ✅ Output format matches original exactly
- ✅ Code reduction of at least 30% in print_header() method
- ✅ All regression tests pass
- ✅ No performance degradation

---

## Metrics and Success Tracking

### Quantitative Metrics
- **Code Lines Reduced**: Target 400+ lines across all migrations
- **Files Modified**: 40 files
- **Average Reduction per File**: 50-70% for simple cases
- **Duplication Eliminated**: 100% (consolidate to single utility)

### Qualitative Metrics
- Improved code readability
- Easier maintenance
- Consistent formatting across modules
- Better developer experience

---

## Recommended Migration Order

### Week 1: Quick Wins (Phase 1)
**Days 1-2**: libfock simple algorithms (DirectJK, DiskJK, PKJK, DirectDFJ, CDJK)
**Days 3-4**: lib3index (DFHelper, DFTensor)
**Day 5**: libfock utilities (VBase, SAP, solvers)

### Week 2: Box-Style Migration (Phase 2)
**Days 1-2**: dfmp2 remaining (UDFMP2, RODFMP2, DFCorrGrad)
**Days 3-5**: libsapt_solver family (SAPT0, USAPT0, SAPT2, SAPT2p, SAPT2p3)

### Week 3: Complex and Utilities (Phases 3-4)
**Days 1-3**: Complex implementations (DiskDFJK, CompositeJK, gradients)
**Days 4-5**: Utility classes (localizers, properties, FISAPT)
**Decision Point**: Evaluate HF::print_header() - migrate or document why not

---

## Alternative Approaches

### Option A: All at Once (Not Recommended)
- **Pros**: Single PR, consistent across codebase immediately
- **Cons**: High risk, difficult to review, hard to debug issues
- **Verdict**: ❌ Too risky

### Option B: Phased Approach (Recommended)
- **Pros**: Incremental validation, easier review, lower risk
- **Cons**: Takes longer, temporary inconsistency
- **Verdict**: ✅ Best balance

### Option C: Only Migrate Simple Cases
- **Pros**: Low risk, quick wins
- **Cons**: Leaves complex cases duplicated
- **Verdict**: ⚠️ Acceptable fallback if time-constrained

---

## Dependencies and Blockers

### None Currently Identified
- HeaderPrinter utility is complete and tested
- No API changes required
- No external dependencies
- Can be done incrementally

### Potential Future Enhancements
- Add table formatting support for complex data
- Integration with Psi4 logging levels
- Template specializations for common types
- Performance profiling (though overhead should be negligible)

---

## Resources

- **HeaderPrinter API**: `psi4/src/psi4/libpsi4util/header_printer.h`
- **Usage Guide**: `psi4/src/psi4/libpsi4util/HEADER_PRINTER_USAGE.md`
- **Example Migrations**: LinK.cc, MemDFJK.cc, mp2.cc (RDFMP2)
- **Original Analysis**: See exploration task output for complete inventory

---

## Next Actions

1. **Review this plan** with team/maintainers
2. **Prioritize phases** based on project needs
3. **Assign ownership** for each phase
4. **Set up tracking** (GitHub issues, project board)
5. **Begin Phase 1** with DirectJK (simplest case)

---

## Appendix: Complete Implementation Inventory

### By Module
- **libfock**: 18 implementations (3 done, 15 remaining)
- **dfmp2**: 4 implementations (1 done, 3 remaining)
- **lib3index**: 2 implementations (0 done, 2 remaining)
- **libscf_solver**: 1 implementation (0 done, 1 remaining)
- **libsapt_solver**: 5 implementations (0 done, 5 remaining)
- **scfgrad**: 2 implementations (0 done, 2 remaining)
- **Other**: 11 implementations (0 done, 11 remaining)

### By Complexity
- **Trivial (1-3 lines)**: 8 implementations
- **Simple (4-10 lines)**: 18 implementations
- **Medium (11-20 lines)**: 12 implementations
- **Complex (21-50 lines)**: 5 implementations

### By Banner Style
- **ARROW**: ~30 implementations
- **BOX**: ~10 implementations
- **NONE**: ~2 implementations
- **CUSTOM/MIXED**: ~1 implementation

---

**Document Version**: 1.0
**Last Updated**: 2025-11-18
**Status**: Ready for Implementation
