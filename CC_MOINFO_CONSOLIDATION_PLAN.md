# CC Module MOInfo Consolidation Plan

## Executive Summary

This document outlines a detailed plan to consolidate 7 duplicate MOInfo structure definitions across Psi4's coupled cluster (CC) modules into a single, shared implementation. This will eliminate approximately 700 lines of duplicate code and improve maintainability.

## Current State Analysis

### Existing MOInfo Locations

1. `psi4/src/psi4/cc/ccenergy/MOInfo.h` (102 lines)
2. `psi4/src/psi4/cc/cclambda/MOInfo.h` (94 lines)
3. `psi4/src/psi4/cc/cceom/MOInfo.h` (97 lines)
4. `psi4/src/psi4/cc/ccdensity/MOInfo.h` (118 lines)
5. `psi4/src/psi4/cc/ccresponse/MOInfo.h` (94 lines)
6. `psi4/src/psi4/cc/cchbar/MOInfo.h` (78 lines)
7. `psi4/src/psi4/cc/cctriples/MOInfo.h` (84 lines)

### Common Fields (Present in all or most modules)

**Basic dimensions:**
- `int nirreps` - number of irreducible representations
- `int nmo` - number of molecular orbitals
- `int nso` - number of symmetry orbitals
- `Dimension orbspi` - MOs per irrep
- `Dimension clsdpi` - closed-shell orbitals per irrep (excluding frozen)
- `Dimension openpi` - open-shell orbitals per irrep
- `Dimension uoccpi` - unoccupied orbitals per irrep (excluding frozen virtual)
- `Dimension frdocc` - frozen core orbitals per irrep
- `Dimension fruocc` - frozen virtual orbitals per irrep
- `std::vector<std::string> labels` - irrep labels

**Orbital counts per irrep:**
- `Dimension occpi` - occupied orbitals (including open)
- `Dimension aoccpi` - alpha occupied orbitals
- `Dimension boccpi` - beta occupied orbitals
- `Dimension virtpi` - virtual orbitals
- `Dimension avirtpi` - alpha virtual orbitals
- `Dimension bvirtpi` - beta virtual orbitals

**Symmetry arrays (raw pointers):**
- `int* occ_sym` - occupied orbital symmetries
- `int* aocc_sym` - alpha occupied orbital symmetries
- `int* bocc_sym` - beta occupied orbital symmetries
- `int* vir_sym` - virtual orbital symmetries
- `int* avir_sym` - alpha virtual orbital symmetries
- `int* bvir_sym` - beta virtual orbital symmetries

**Offset arrays (raw pointers):**
- `int* occ_off` - occupied orbital offsets within each irrep
- `int* aocc_off` - alpha occupied orbital offsets
- `int* bocc_off` - beta occupied orbital offsets
- `int* vir_off` - virtual orbital offsets
- `int* avir_off` - alpha virtual orbital offsets
- `int* bvir_off` - beta virtual orbital offsets

**Energy values:**
- `double enuc` - nuclear repulsion energy
- `double escf` - SCF energy (from wavefunction)
- `double eref` - reference energy (from file)
- `double ecc` - coupled cluster energy

### Module-Specific Fields

**ccenergy (most comprehensive):**
- `int nao` - number of atomic orbitals
- `Dimension sopi` - SOs per irrep
- `int* sosym` - SO symmetry (Pitzer ordering)
- `int nvirt` - total number of virtual orbitals
- `int iter` - current CCSD iteration
- `double conv` - current convergence level
- `double emp2` - MP2 energy
- `double emp2_ss` - same-spin MP2 correlation energy
- `double emp2_os` - opposite-spin MP2 correlation energy
- `double emp2_s` - singles MP2 correlation energy
- `double ecc_ss` - same-spin CC correlation energy
- `double ecc_os` - opposite-spin CC correlation energy
- `double ecc_s` - singles CC correlation energy
- `double t1diag` - T1 diagnostic
- `double d1diag` - D1 diagnostic (Janssen-Nielsen)
- `double new_d1diag` - modified D1 diagnostic (Lee)
- `double d2diag` - D2 diagnostic (Nielsen-Janssen)
- `double*** Cv, ***Cav, ***Cbv` - virtual orbital transformation matrices
- `double*** Co, ***Cao, ***Cbo` - occupied orbital transformation matrices

**cclambda:**
- `int nao`, `Dimension sopi`, `int* sosym`, `int nvirt` (same as ccenergy)
- `int iter` - current lambda iteration
- `int sym` - symmetry of converged CCSD state
- `double conv` - convergence level
- `double lcc` - lambda pseudoenergy
- `double*** C, ***Ca, ***Cb` - orbital transformation matrices (different naming)

**cceom:**
- `int iopen` - open/closed shell flag (0=closed, >0=open)
- `int nvirt`, `Dimension sopi`, `int* sosym` (same as above)
- `int iter`, `int sym`, `double conv` (same as cclambda)
- `double t1diag`, `double d1diag` (diagnostics)
- `char** irr_labs_lowercase` - lowercase irrep labels
- `double*** C, ***Ca, ***Cb` - transformation matrices

**ccdensity:**
- `int nactive` - number of active orbitals
- `int nfzc`, `int nfzv` - frozen core/virtual counts
- `int nclsd`, `int nopen`, `int nuocc` - total counts
- `int sym` - state symmetry
- `double et` - (T) energy from cctriples
- Uses `std::vector<int>` instead of raw pointers for symmetry/offset arrays
- Extensive reordering arrays:
  - `std::vector<int> cc_occ, cc_aocc, cc_bocc, cc_vir, cc_avir, cc_bvir`
  - `std::vector<int> qt_occ, qt_aocc, qt_bocc, qt_vir, qt_avir, qt_bvir`
  - `std::vector<int> pitzer2qt, qt2pitzer`
- Matrix objects for densities:
  - `Matrix opdm, opdm_a, opdm_b` - one-particle density matrices
  - `Matrix ltd_mat, ltd_a_mat, ltd_b_mat` - left transition densities
  - `Matrix rtd_mat, rtd_a_mat, rtd_b_mat` - right transition densities
- `SharedMatrix Ca` - SCF orbitals
- `std::vector<SharedMatrix> L, nabla, dip`

**ccresponse:**
- `int nao`, `int noei`, `int ntri`, `int noei_ao` - various triangle sizes
- `int nactive`, `int nfzc`, `int nvirt` - orbital counts
- `Dimension sopi`, `Dimension actpi`, `Dimension act_occpi`
- `int* occ_sym, *aocc_sym, *bocc_sym, *vir_sym, *avir_sym, *bvir_sym` (raw pointers)
- `int* occ_off, *aocc_off, *bocc_off, *vir_off, *avir_off, *bvir_off`
- `std::shared_ptr<Matrix> Ca` - SCF orbitals
- `int* mu_irreps` - dipole component irreps (x,y,z)
- `int* l_irreps` - angular momentum component irreps
- `int natom` - number of atoms
- `double* zvals` - atomic Z values
- `double*** C` - virtual orbital transformation matrix

**cchbar (minimal):**
- Only basic structure: dimensions, symmetries, offsets
- No energy values or diagnostics

**cctriples:**
- `int iter`, `double conv` - iteration info
- `double enuc`, `double escf`, `double eref`, `double ecc` - energies
- Basic orbital dimensions and symmetries

### Why Not Use Existing libmoinfo?

The existing `psi4/src/psi4/libmoinfo/moinfo.h`:
- Is class-based with complex object-oriented design
- Designed for multi-reference methods with SlaterDeterminant support
- Uses different orbital space concepts (focc, docc, actv, extr, fvir, all)
- Would require massive refactoring of all CC modules
- Would break existing multi-reference functionality

**Decision:** Create a new, CC-specific unified MOInfo structure.

---

## Proposed Solution: Unified CC MOInfo

### Architecture Design

Create a new shared CC MOInfo implementation in: `psi4/src/psi4/cc/ccmoinfo/`

**Design Principles:**
1. **RAII (Resource Acquisition Is Initialization)**: Use smart pointers and RAII to manage memory automatically
2. **Backward compatibility**: Provide same interface as current structs where possible
3. **Extensibility**: Support module-specific fields through optional/conditional members
4. **Type safety**: Replace raw pointers with smart pointers or std::vector where appropriate
5. **Clear ownership**: Use unique_ptr for owned data, shared_ptr only when needed

### Proposed Class Structure

```cpp
// psi4/src/psi4/cc/ccmoinfo/CCMOInfo.h

namespace psi {
namespace ccmoinfo {

class CCMOInfo {
public:
    // Constructor from wavefunction
    CCMOInfo(std::shared_ptr<Wavefunction> wfn, int reference_type);

    // Destructor (handles cleanup automatically)
    ~CCMOInfo() = default;

    // Initialization method (reads from PSIO)
    void initialize();

    // Common basic dimensions (always present)
    int nirreps;
    int nmo;
    int nso;
    int nao;  // optional, only set if needed

    Dimension sopi;
    Dimension orbspi;
    Dimension clsdpi;
    Dimension openpi;
    Dimension uoccpi;
    Dimension frdocc;
    Dimension fruocc;

    std::vector<std::string> labels;

    // Orbital counts per irrep (always present)
    Dimension occpi;
    Dimension aoccpi;
    Dimension boccpi;
    Dimension virtpi;
    Dimension avirtpi;
    Dimension bvirtpi;

    int nvirt;  // total virtual orbitals

    // Symmetry arrays (use vectors instead of raw pointers)
    std::vector<int> occ_sym;
    std::vector<int> aocc_sym;
    std::vector<int> bocc_sym;
    std::vector<int> vir_sym;
    std::vector<int> avir_sym;
    std::vector<int> bvir_sym;
    std::vector<int> sosym;

    // Offset arrays (use vectors instead of raw pointers)
    std::vector<int> occ_off;
    std::vector<int> aocc_off;
    std::vector<int> bocc_off;
    std::vector<int> vir_off;
    std::vector<int> avir_off;
    std::vector<int> bvir_off;

    // Common energies
    double enuc;
    double escf;
    double eref;
    double ecc;

    // Iteration tracking (optional)
    int iter = 0;
    double conv = 0.0;
    int sym = 0;  // symmetry of state

    // MP2 energies (optional, for ccenergy)
    double emp2 = 0.0;
    double emp2_ss = 0.0;
    double emp2_os = 0.0;
    double emp2_s = 0.0;

    // CC energy components (optional, for ccenergy)
    double ecc_ss = 0.0;
    double ecc_os = 0.0;
    double ecc_s = 0.0;

    // Diagnostics (optional, for ccenergy/cceom)
    double t1diag = 0.0;
    double d1diag = 0.0;
    double new_d1diag = 0.0;
    double d2diag = 0.0;

    // Lambda pseudoenergy (optional, for cclambda)
    double lcc = 0.0;

    // Triples energy (optional, for ccdensity)
    double et = 0.0;

    // Open shell flag (optional, for cceom)
    int iopen = 0;

    // Active orbital counts (optional, for ccdensity/ccresponse)
    int nactive = 0;
    int nfzc = 0;
    int nfzv = 0;
    int nclsd = 0;
    int nopen_tot = 0;
    int nuocc = 0;

    // Triangle sizes (optional, for ccresponse)
    int noei = 0;
    int ntri = 0;
    int noei_ao = 0;

    // Active dimensions (optional, for ccresponse)
    Dimension actpi;
    Dimension act_occpi;

    // Orbital transformation matrices (optional, use smart pointers)
    // For modules that need them
    std::unique_ptr<double***> Cv;   // virtual transformation
    std::unique_ptr<double***> Cav;  // alpha virtual
    std::unique_ptr<double***> Cbv;  // beta virtual
    std::unique_ptr<double***> Co;   // occupied transformation
    std::unique_ptr<double***> Cao;  // alpha occupied
    std::unique_ptr<double***> Cbo;  // beta occupied

    // Alternative naming for transformation matrices (cclambda/cceom/ccresponse)
    std::unique_ptr<double***> C;    // virtual transformation
    std::unique_ptr<double***> Ca;   // alpha virtual
    std::unique_ptr<double***> Cb;   // beta virtual

    // Reordering arrays (optional, for ccdensity)
    std::vector<int> cc_occ, cc_aocc, cc_bocc;
    std::vector<int> cc_vir, cc_avir, cc_bvir;
    std::vector<int> qt_occ, qt_aocc, qt_bocc;
    std::vector<int> qt_vir, qt_avir, qt_bvir;
    std::vector<int> pitzer2qt, qt2pitzer;

    // Density matrices (optional, for ccdensity)
    std::shared_ptr<Matrix> opdm, opdm_a, opdm_b;
    std::shared_ptr<Matrix> ltd_mat, ltd_a_mat, ltd_b_mat;
    std::shared_ptr<Matrix> rtd_mat, rtd_a_mat, rtd_b_mat;

    // SCF orbitals (optional)
    std::shared_ptr<Matrix> Ca_matrix;
    std::vector<std::shared_ptr<Matrix>> L, nabla, dip;

    // Property-related (optional, for ccresponse)
    std::vector<int> mu_irreps;  // dipole component irreps
    std::vector<int> l_irreps;   // angular momentum component irreps
    int natom = 0;
    std::vector<double> zvals;

    // Lowercase labels (optional, for cceom)
    std::vector<std::string> irr_labs_lowercase;

private:
    std::shared_ptr<Wavefunction> wfn_;
    int reference_type_;

    // Helper methods for initialization
    void read_common_data();
    void read_orbital_info();
    void allocate_transformation_matrices();
};

} // namespace ccmoinfo
} // namespace psi
```

### Backward Compatibility Wrapper

To minimize code changes in existing modules, we can provide compatibility accessors:

```cpp
// For modules that use raw pointers, provide getters
int* get_occ_sym() { return occ_sym.data(); }
int* get_vir_sym() { return vir_sym.data(); }
// etc.
```

### Memory Management Strategy

**Current approach (problematic):**
- Manual allocation with `init_int_array()`
- Manual deallocation (often missing or error-prone)
- Memory leaks if exceptions occur

**New approach:**
- Use `std::vector<int>` for all arrays
- Use `std::unique_ptr` or `std::shared_ptr` for complex structures
- Automatic cleanup when CCMOInfo goes out of scope
- Exception-safe

---

## Implementation Plan

### Phase 1: Preparation and Analysis

**Tasks:**
1. ✓ Document current state (all MOInfo structures examined)
2. Create detailed field inventory with usage patterns
3. Identify dependencies between modules
4. Design test strategy to verify correctness

**Deliverables:**
- This planning document
- Field inventory spreadsheet/table
- Test plan document

**Estimated Effort:** 1 day

---

### Phase 2: Create Unified CCMOInfo Infrastructure

**Tasks:**

1. **Create directory structure**
   ```bash
   mkdir -p psi4/src/psi4/cc/ccmoinfo
   ```

2. **Create CCMOInfo.h**
   - Define unified class structure
   - Include all common fields
   - Include optional module-specific fields
   - Add appropriate comments and documentation

3. **Create CCMOInfo.cc**
   - Implement constructor
   - Implement initialization from wavefunction
   - Implement PSIO reading logic
   - Consolidate common code from all 7 get_moinfo.cc files
   - Add proper error handling

4. **Create CMakeLists.txt for ccmoinfo**
   ```cmake
   set(ccmoinfo_sources
       CCMOInfo.cc
   )

   add_library(ccmoinfo OBJECT ${ccmoinfo_sources})
   target_include_directories(ccmoinfo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
   ```

5. **Update psi4/src/psi4/cc/CMakeLists.txt**
   - Add ccmoinfo subdirectory
   - Link ccmoinfo to cc target
   ```cmake
   add_subdirectory(ccmoinfo)
   target_sources(cc PRIVATE $<TARGET_OBJECTS:ccmoinfo>)
   ```

**Deliverables:**
- `psi4/src/psi4/cc/ccmoinfo/CCMOInfo.h`
- `psi4/src/psi4/cc/ccmoinfo/CCMOInfo.cc`
- `psi4/src/psi4/cc/ccmoinfo/CMakeLists.txt`
- Updated `psi4/src/psi4/cc/CMakeLists.txt`

**Estimated Effort:** 2-3 days

---

### Phase 3: Module Migration (Incremental)

Migrate modules one at a time, starting with simplest. This allows testing at each step.

**Migration Order (simplest to most complex):**

1. **cchbar** (simplest - minimal fields)
2. **cctriples** (simple - few special fields)
3. **cclambda** (moderate - has lcc and sym)
4. **cceom** (moderate - has diagnostics)
5. **ccresponse** (complex - many special fields)
6. **ccenergy** (most complex - all diagnostics and MP2)
7. **ccdensity** (special - uses vectors and matrices differently)

**For Each Module:**

**Step 3.1: Update includes**
```cpp
// Old:
#include "MOInfo.h"

// New:
#include "psi4/cc/ccmoinfo/CCMOInfo.h"
using psi::ccmoinfo::CCMOInfo;
```

**Step 3.2: Update type references**
```cpp
// Old:
struct MOInfo { ... };

// New:
// (use CCMOInfo from ccmoinfo)
```

**Step 3.3: Update initialization code**
- Modify module's main class to use CCMOInfo
- Update get_moinfo() or equivalent to use CCMOInfo::initialize()
- Migrate any module-specific initialization to CCMOInfo

**Step 3.4: Update field references**
- If using raw pointers, add `.data()` calls for vector access
- Update transformation matrix access if needed
- Search-and-replace field name changes

**Step 3.5: Update CMakeLists.txt**
```cmake
# Add dependency on ccmoinfo
target_include_directories(ccXXXXX PRIVATE ${CMAKE_SOURCE_DIR}/psi4/cc/ccmoinfo)
```

**Step 3.6: Remove old files**
- Delete old `MOInfo.h`
- Delete old `get_moinfo.cc` (if fully consolidated)
- OR keep get_moinfo.cc but simplify to just call CCMOInfo methods

**Step 3.7: Compile and test**
```bash
cd build
make ccXXXXX
ctest -R ccXXXXX -V
```

**Deliverables (per module):**
- Updated source files
- Updated CMakeLists.txt
- Passing tests
- Git commit for that module

**Estimated Effort:** 1-2 days per module = 7-14 days total

---

### Phase 4: Testing and Validation

**Test Strategy:**

1. **Unit Tests**
   - Create unit tests for CCMOInfo initialization
   - Test all field population
   - Test different reference types (RHF, ROHF, UHF)

2. **Integration Tests**
   - Run existing CC test suite
   - Verify numerical results unchanged
   - Check for memory leaks with valgrind

3. **Regression Tests**
   - Compare outputs before/after refactoring
   - Ensure bit-identical results for standard tests

**Test Commands:**
```bash
# Run all CC tests
ctest -R cc -V

# Run specific module tests
ctest -R ccenergy -V
ctest -R cclambda -V
# etc.

# Memory leak check (if available)
valgrind --leak-check=full ./ccenergy_test
```

**Deliverables:**
- Test results documentation
- Memory leak analysis
- Performance comparison (if significant changes)

**Estimated Effort:** 2-3 days

---

### Phase 5: Documentation and Cleanup

**Tasks:**

1. **Code Documentation**
   - Add Doxygen comments to CCMOInfo class
   - Document each field's purpose
   - Add usage examples

2. **Update User Documentation**
   - Update developer documentation
   - Add migration guide for developers

3. **Final Cleanup**
   - Remove all old MOInfo.h files
   - Remove old get_moinfo.cc files (if consolidated)
   - Clean up any dead code
   - Run code formatter (if project uses one)

4. **Create Summary Document**
   - Lines of code removed
   - Modules affected
   - Breaking changes (if any)
   - Migration notes

**Deliverables:**
- Documented CCMOInfo class
- Updated developer docs
- Migration guide
- Summary report

**Estimated Effort:** 1-2 days

---

### Phase 6: Code Review and Merge

**Tasks:**

1. **Create Pull Request**
   - Write comprehensive PR description
   - Link to this planning document
   - Highlight testing results

2. **Code Review**
   - Address reviewer feedback
   - Make requested changes
   - Re-test as needed

3. **Final Merge**
   - Merge to development branch
   - Monitor CI/CD pipeline
   - Address any post-merge issues

**Deliverables:**
- Pull request
- Code review responses
- Merged code

**Estimated Effort:** 1-3 days (depends on review feedback)

---

## Timeline Estimate

| Phase | Duration | Cumulative |
|-------|----------|------------|
| Phase 1: Preparation | 1 day | 1 day |
| Phase 2: Infrastructure | 2-3 days | 3-4 days |
| Phase 3: Migration (7 modules) | 7-14 days | 10-18 days |
| Phase 4: Testing | 2-3 days | 12-21 days |
| Phase 5: Documentation | 1-2 days | 13-23 days |
| Phase 6: Review & Merge | 1-3 days | 14-26 days |

**Total Estimated Effort: 14-26 days (3-5 weeks)**

Note: This assumes one person working full-time. Parallel development could reduce timeline.

---

## Risk Assessment

### High-Priority Risks

1. **Risk:** Breaking existing functionality during migration
   - **Mitigation:** Incremental migration with testing at each step
   - **Contingency:** Keep old code until all tests pass

2. **Risk:** Module-specific fields not properly identified
   - **Mitigation:** Thorough analysis phase; make fields optional
   - **Contingency:** Add fields as discovered during testing

3. **Risk:** Performance degradation from abstraction
   - **Mitigation:** Use inline accessors; benchmark critical paths
   - **Contingency:** Optimize or provide low-level accessors

### Medium-Priority Risks

4. **Risk:** Memory management issues with pointer conversions
   - **Mitigation:** Use smart pointers; thorough memory leak testing
   - **Contingency:** Provide raw pointer accessors if needed

5. **Risk:** Build system complications
   - **Mitigation:** Simple CMake structure; test on multiple platforms
   - **Contingency:** Keep build changes minimal

### Low-Priority Risks

6. **Risk:** Incomplete documentation
   - **Mitigation:** Document as you go; peer review
   - **Contingency:** Add documentation in follow-up PR

---

## Success Criteria

1. ✅ All 7 CC modules use unified CCMOInfo
2. ✅ All existing tests pass
3. ✅ No memory leaks introduced
4. ✅ At least 600 lines of duplicate code removed
5. ✅ Code properly documented
6. ✅ Build time not significantly increased
7. ✅ No performance regressions

---

## Alternative Approaches Considered

### Alternative 1: Minimal Refactoring
**Approach:** Keep separate MOInfo.h files but extract common initialization code

**Pros:**
- Less invasive
- Lower risk

**Cons:**
- Doesn't solve duplication problem
- Still have to maintain 7 separate structures
- Future changes require 7 updates

**Decision:** Rejected - doesn't achieve goal

### Alternative 2: Template-Based Approach
**Approach:** Use C++ templates to share common structure with module-specific extensions

**Pros:**
- Compile-time customization
- No runtime overhead

**Cons:**
- More complex to understand
- Harder to debug
- Overkill for this use case

**Decision:** Rejected - unnecessary complexity

### Alternative 3: Inheritance-Based Approach
**Approach:** Create base MOInfo class with derived classes for each module

**Pros:**
- Clear separation of common/specific fields
- Follows OOP principles

**Cons:**
- More complex class hierarchy
- Potential virtual function overhead
- Each module still has separate class

**Decision:** Rejected - adds complexity without major benefit

### Alternative 4: Extend libmoinfo
**Approach:** Add CC-specific functionality to existing libmoinfo

**Pros:**
- Uses existing infrastructure
- One MOInfo for entire codebase

**Cons:**
- libmoinfo designed for multi-reference methods
- Would require major refactoring of libmoinfo
- Risk breaking multi-reference codes
- Architectural mismatch (class-based vs struct-based)

**Decision:** Rejected - too risky and incompatible designs

---

## Appendix A: Field Inventory Table

| Field | ccenergy | cclambda | cceom | ccdensity | ccresponse | cchbar | cctriples | Type | Notes |
|-------|----------|----------|-------|-----------|------------|--------|-----------|------|-------|
| nirreps | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | int | Always present |
| nmo | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | int | Always present |
| nso | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | int | Common |
| nao | ✓ | ✓ | - | - | ✓ | - | - | int | ccenergy, cclambda, ccresponse |
| sopi | ✓ | ✓ | ✓ | - | ✓ | - | - | Dimension | Most modules |
| sosym | ✓ | ✓ | ✓ | - | - | - | - | int* | ccenergy, cclambda, cceom |
| orbspi | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | Dimension | Always present |
| clsdpi | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | Dimension | Always present |
| openpi | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | Dimension | Always present |
| uoccpi | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | Dimension | Always present |
| frdocc | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | Dimension | Always present |
| fruocc | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | Dimension | Always present |
| nvirt | ✓ | ✓ | ✓ | - | ✓ | - | - | int | Most modules |
| labels | ✓ | ✓ | - | ✓ | ✓ | ✓ | ✓ | vector<string> | Most (cceom uses irr_labs) |
| irr_labs | - | - | ✓ | - | - | - | - | vector<string> | cceom only |
| irr_labs_lowercase | - | - | ✓ | - | - | - | - | char** | cceom only |
| occpi | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | Dimension | Always present |
| aoccpi | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | Dimension | Always present |
| boccpi | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | Dimension | Always present |
| virtpi | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | Dimension | Always present |
| avirtpi | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | Dimension | Always present |
| bvirtpi | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | Dimension | Always present |
| occ_sym | ✓ | ✓ | ✓ | ✓(vec) | ✓ | ✓ | ✓ | int* | Always (ccdensity uses vector) |
| aocc_sym | ✓ | ✓ | ✓ | ✓(vec) | ✓ | ✓ | ✓ | int* | Always (ccdensity uses vector) |
| bocc_sym | ✓ | ✓ | ✓ | ✓(vec) | ✓ | ✓ | ✓ | int* | Always (ccdensity uses vector) |
| vir_sym | ✓ | ✓ | ✓ | ✓(vec) | ✓ | ✓ | ✓ | int* | Always (ccdensity uses vector) |
| avir_sym | ✓ | ✓ | ✓ | ✓(vec) | ✓ | ✓ | ✓ | int* | Always (ccdensity uses vector) |
| bvir_sym | ✓ | ✓ | ✓ | ✓(vec) | ✓ | ✓ | ✓ | int* | Always (ccdensity uses vector) |
| occ_off | ✓ | ✓ | ✓ | ✓(vec) | ✓ | ✓ | ✓ | int* | Always (ccdensity uses vector) |
| aocc_off | ✓ | ✓ | ✓ | ✓(vec) | ✓ | ✓ | ✓ | int* | Always (ccdensity uses vector) |
| bocc_off | ✓ | ✓ | ✓ | ✓(vec) | ✓ | ✓ | ✓ | int* | Always (ccdensity uses vector) |
| vir_off | ✓ | ✓ | ✓ | ✓(vec) | ✓ | ✓ | ✓ | int* | Always (ccdensity uses vector) |
| avir_off | ✓ | ✓ | ✓ | ✓(vec) | ✓ | ✓ | ✓ | int* | Always (ccdensity uses vector) |
| bvir_off | ✓ | ✓ | ✓ | ✓(vec) | ✓ | ✓ | ✓ | int* | Always (ccdensity uses vector) |
| iter | ✓ | ✓ | ✓ | - | - | - | ✓ | int | Iteration counter |
| conv | ✓ | ✓ | ✓ | - | - | - | ✓ | double | Convergence level |
| sym | - | ✓ | ✓ | ✓ | - | - | - | int | State symmetry |
| enuc | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | double | Nuclear repulsion |
| escf | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | double | SCF energy |
| eref | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | double | Reference energy |
| ecc | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | double | CC energy |
| emp2 | ✓ | - | - | - | - | - | - | double | ccenergy only |
| emp2_ss | ✓ | - | - | - | - | - | - | double | ccenergy only |
| emp2_os | ✓ | - | - | - | - | - | - | double | ccenergy only |
| emp2_s | ✓ | - | - | - | - | - | - | double | ccenergy only |
| ecc_ss | ✓ | - | - | - | - | - | - | double | ccenergy only |
| ecc_os | ✓ | - | - | - | - | - | - | double | ccenergy only |
| ecc_s | ✓ | - | - | - | - | - | - | double | ccenergy only |
| t1diag | ✓ | - | ✓ | - | - | - | - | double | Diagnostics |
| d1diag | ✓ | - | ✓ | - | - | - | - | double | Diagnostics |
| new_d1diag | ✓ | - | - | - | - | - | - | double | ccenergy only |
| d2diag | ✓ | - | - | - | - | - | - | double | ccenergy only |
| lcc | - | ✓ | - | - | - | - | - | double | Lambda pseudoenergy |
| et | - | - | - | ✓ | - | - | - | double | (T) energy |
| iopen | - | - | ✓ | - | - | - | - | int | cceom only |
| nactive | - | - | - | ✓ | ✓ | - | - | int | ccdensity, ccresponse |
| Cv, Cav, Cbv | ✓ | - | - | - | - | - | - | double*** | Virtual orb transform |
| Co, Cao, Cbo | ✓ | - | - | - | - | - | - | double*** | Occupied orb transform |
| C, Ca, Cb | - | ✓ | ✓ | - | ✓ | - | - | double*** | Alternative naming |

(Additional ccdensity/ccresponse specific fields omitted for brevity - see main document)

---

## Appendix B: Code Migration Checklist

Use this checklist for each module migration:

- [ ] Create feature branch
- [ ] Update includes to use CCMOInfo.h
- [ ] Update type declarations
- [ ] Modify initialization code
- [ ] Update all field references
- [ ] Add `.data()` for vector-to-pointer conversions if needed
- [ ] Update CMakeLists.txt dependencies
- [ ] Compile module successfully
- [ ] Run module unit tests
- [ ] Run module integration tests
- [ ] Check for memory leaks
- [ ] Verify numerical results unchanged
- [ ] Remove old MOInfo.h
- [ ] Clean up get_moinfo.cc
- [ ] Update module documentation
- [ ] Code review
- [ ] Commit changes
- [ ] Merge to development branch

---

## Appendix C: Build Commands Reference

```bash
# Create build directory (if not exists)
mkdir -p build && cd build

# Configure with CMake
cmake ..

# Build specific module
make ccenergy
make cclambda
make cceom
make ccdensity
make ccresponse
make cchbar
make cctriples

# Build all CC modules
make cc

# Run tests for specific module
ctest -R ccenergy -V
ctest -R cclambda -V
# etc.

# Run all CC tests
ctest -R cc -V

# Clean build
make clean

# Rebuild from scratch
rm -rf * && cmake .. && make -j4
```

---

## Questions for Stakeholders

Before proceeding with implementation, please address:

1. **Backward Compatibility:** Do we need to maintain API compatibility for external code that might use these MOInfo structures?

2. **Naming Convention:** Should we keep CCMOInfo or use a different name (e.g., UnifiedMOInfo, SharedMOInfo)?

3. **Testing Resources:** What level of testing is expected? Do we have access to reference data for validation?

4. **Performance Requirements:** Are there any performance-critical code paths that we need to benchmark?

5. **Release Timeline:** Is there a target release for this refactoring?

6. **Code Review Process:** Who should review the code at each phase?

---

**Document Version:** 1.0
**Author:** Claude (Anthropic AI Assistant)
**Date:** 2025-11-18
**Status:** DRAFT - Awaiting Approval
