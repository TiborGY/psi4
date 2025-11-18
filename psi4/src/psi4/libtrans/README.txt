================================================================================
                    LIBTRANS: Integral Transformation Library
                              Usage Guide and Best Practices
================================================================================

TABLE OF CONTENTS
1. Overview
2. When to Use libtrans
3. Quick Start Guide
4. IntegralTransform Class Reference
5. MOSpace: Defining Orbital Spaces
6. Common Usage Patterns
7. Performance Optimization
8. Advanced Features
9. Frozen Core Handling
10. Common Pitfalls and Solutions
11. When NOT to Use libtrans

================================================================================
1. OVERVIEW
================================================================================

The libtrans library provides a comprehensive framework for transforming
one- and two-electron integrals from the atomic orbital (AO) basis to the
molecular orbital (MO) basis. It is the CANONICAL implementation for integral
transformations in Psi4.

Key Features:
• Full AO → MO transformation for arbitrary orbital spaces
• Automated frozen core handling (operator construction and energy)
• Support for restricted, unrestricted, and semicanonical orbitals
• Half-transformation capabilities for memory efficiency
• Multiple output formats (DPD, IWL)
• Flexible orbital ordering (QT order, Pitzer order)

Main Classes:
• IntegralTransform - Primary transformation interface
• MOSpace - Defines orbital subspaces for transformations

Files:
• integraltransform.h/.cc - Main transformation class
• integraltransform_tei.cc - Two-electron integral transformations
• integraltransform_tei_1st_half.cc - First half of transformation
• integraltransform_tei_2nd_half.cc - Second half of transformation
• integraltransform_oei.cc - One-electron integral transformations
• integraltransform_sort_so_tei.cc - SO presorting and frozen core setup
• mospace.h/.cc - Orbital space definitions

================================================================================
2. WHEN TO USE LIBTRANS
================================================================================

USE libtrans if you are:
✓ Implementing single-reference correlation methods (MP2, CCSD, etc.)
✓ Need standard frozen core handling
✓ Transforming integrals between well-defined orbital spaces
✓ Want DPD (Distributed Packed Dense) formatted output
✓ Want IWL (Integrals With Labels) formatted output
✓ Need to handle restricted/unrestricted/semicanonical orbitals

DO NOT use libtrans if:
✗ Using density-fitting (see dfocc for DF-specific implementations)
✗ Implementing multi-reference methods with custom blocking (see psimrcc)
✗ Need fundamentally different integral access patterns
✗ Have verified >10% performance improvement with custom code

If you're unsure, START with libtrans. Only write custom code if profiling
shows a clear performance bottleneck that cannot be addressed within libtrans.

================================================================================
3. QUICK START GUIDE
================================================================================

Minimal Example (Restricted, Correlation Method):
-------------------------------------------------

#include "psi4/libtrans/integraltransform.h"
#include "psi4/libtrans/mospace.h"

// Define orbital spaces
std::vector<std::shared_ptr<MOSpace>> spaces;
spaces.push_back(MOSpace::occ);  // Occupied orbitals
spaces.push_back(MOSpace::vir);  // Virtual orbitals

// Create transformation object
auto ints = std::make_shared<IntegralTransform>(
    wfn,                                            // Wavefunction
    spaces,                                         // Orbital spaces
    IntegralTransform::TransformationType::Restricted,
    IntegralTransform::OutputType::DPDOnly,
    IntegralTransform::MOOrdering::QTOrder,
    IntegralTransform::FrozenOrbitals::OccOnly,
    false                                           // Don't initialize yet
);

// Configure options
ints->set_print(1);                     // Print level
ints->set_dpd_id(0);                    // DPD instance
ints->set_keep_dpd_so_ints(false);      // Don't keep SO ints after use
ints->initialize();                      // Now initialize

// Transform integrals
ints->update_orbitals();                // Sync orbitals from wavefunction
ints->transform_tei(MOSpace::occ, MOSpace::occ,
                    MOSpace::vir, MOSpace::vir);

// Get frozen core energy if needed
double E_fzc = ints->get_frozen_core_energy();

================================================================================
4. INTEGRALTRANSFORM CLASS REFERENCE
================================================================================

Constructor Parameters:
----------------------

IntegralTransform(
    std::shared_ptr<Wavefunction> wfn,
    SpaceVec spaces,
    TransformationType transformationType = Restricted,
    OutputType outputType = DPDOnly,
    MOOrdering moOrdering = QTOrder,
    FrozenOrbitals frozenOrbitals = OccAndVir,
    bool initialize = true
);

Parameters Explained:
--------------------

1. TransformationType:
   • Restricted    - Same alpha and beta orbitals (RHF reference)
   • Unrestricted  - Different alpha and beta orbitals (UHF reference)
   • SemiCanonical - Start from restricted, diagonalize occ/vir blocks
                     separately (ROHF → semicanonical)

2. OutputType:
   • DPDOnly    - Write integrals to DPD structure only (RECOMMENDED)
   • IWLOnly    - Write integrals to IWL structure only
   • IWLAndDPD  - Write to both formats (use sparingly, doubles storage)

3. MOOrdering:
   • QTOrder     - Order by class: fzc < occ < socc < vir < fzv
                   then by irrep within classes
   • PitzerOrder - Order by irrep, then by orbital energy within irreps

   Note: Only affects IWL output. DPD output is always space-based.

4. FrozenOrbitals:
   • None      - No orbitals excluded
   • OccOnly   - Exclude frozen occupied orbitals only
   • VirOnly   - Exclude frozen virtual orbitals only
   • OccAndVir - Exclude both frozen occ and frozen vir (STANDARD for CC)

   IMPORTANT: Orbitals are only frozen if the wavefunction detected them.
   Requesting frozen-core with no frozen orbitals → all-electron transform.

Key Methods:
-----------

void initialize()
  • Call after setting all options (if initialize=false in constructor)

void update_orbitals()
  • Sync orbitals from wavefunction to transformation object
  • ALWAYS call before each transformation if orbitals have changed

void presort_so_tei()
  • Presort SO basis integrals for efficient transformation
  • Constructs frozen core operator and computes frozen core energy
  • Called automatically by transform_tei() if needed

void transform_tei(s1, s2, s3, s4, HalfTrans=MakeAndNuke)
  • Transform (s1 s2|s3 s4) integrals
  • See section 6 for HalfTrans options

void transform_tei_first_half(s1, s2)
  • Transform first two indices only: (s1 s2|AO AO) → intermediates
  • Use when sharing first half across multiple final transforms

void transform_tei_second_half(s1, s2, s3, s4)
  • Complete transformation using existing first-half intermediates
  • Must match the s1, s2 from a previous first_half call

double get_frozen_core_energy()
  • Returns energy contribution from frozen core orbitals only
  • Computed during presort_so_tei()
  • Use to adjust total energies for correlation methods

Setters (call before initialize()):
----------------------------------

set_print(int level)               - Print verbosity (0-6)
set_dpd_id(int id)                 - DPD instance number (usually 0)
set_memory(size_t MB)              - Available memory in MB
set_keep_dpd_so_ints(bool)         - Keep SO integrals after transformation
set_keep_iwl_so_ints(bool)         - Keep IWL SO integrals
set_keep_ht_ints(bool)             - Keep half-transformed intermediates
set_tei_already_presorted(bool)    - Skip presorting if already done

================================================================================
5. MOSPACE: DEFINING ORBITAL SPACES
================================================================================

Predefined Spaces:
-----------------

MOSpace::occ  - Occupied orbitals (label 'O')
                Handles frozen orbitals per FrozenOrbitals setting
                Restricted: DOCC + SOCC
                Unrestricted: Occupied alpha or beta

MOSpace::vir  - Virtual orbitals (label 'V')
                Handles frozen orbitals per FrozenOrbitals setting
                Restricted: SOCC + Virtual
                Unrestricted: Virtual alpha or beta

MOSpace::fzc  - Frozen core (label 'o')
MOSpace::fzv  - Frozen virtual (label 'v')
MOSpace::all  - All active MOs (label 'A')
MOSpace::nil  - AO basis (label 'n')

Custom Spaces:
-------------

For specialized orbital subsets, create custom MOSpace objects:

MOSpace(char label,
        const std::vector<int> orbsPI)  // Orbitals per irrep

MOSpace(char label,
        const std::vector<int> aOrbs,    // Alpha orbital indices
        const std::vector<int> bOrbs,    // Beta orbital indices
        const std::vector<int> aIndex,   // Alpha reindexing
        const std::vector<int> bIndex)   // Beta reindexing

Labels:
  Use unique single characters. Reserved labels: o O v V A n d

Example - Active space for CASSCF:
  std::vector<int> active_per_irrep = {2, 0, 1, 1};
  auto active = std::make_shared<MOSpace>('a', active_per_irrep);

================================================================================
6. COMMON USAGE PATTERNS
================================================================================

Pattern 1: Simple Transform (All at Once)
------------------------------------------

// For methods needing only (ij|ab) type integrals
ints->update_orbitals();
ints->transform_tei(MOSpace::occ, MOSpace::occ,
                    MOSpace::vir, MOSpace::vir);

// Access via DPD:
dpdbuf4 K;
global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0,
                       ID("[O,O]"), ID("[V,V]"),
                       ID("[O>=O]+"), ID("[V>=V]+"), 0,
                       "MO Ints (OO|VV)");
// ... use integrals ...
global_dpd_->buf4_close(&K);


Pattern 2: Multiple Transforms with Half-Transform Reuse
--------------------------------------------------------

// First transform: (OO|OO) - MakeAndKeep the first half
ints->transform_tei(MOSpace::occ, MOSpace::occ,
                    MOSpace::occ, MOSpace::occ,
                    IntegralTransform::HalfTrans::MakeAndKeep);

// Reuse (OO| first half for (OO|OV)
ints->transform_tei(MOSpace::occ, MOSpace::occ,
                    MOSpace::occ, MOSpace::vir,
                    IntegralTransform::HalfTrans::ReadAndKeep);

// Reuse again for (OO|VV), then delete
ints->transform_tei(MOSpace::occ, MOSpace::occ,
                    MOSpace::vir, MOSpace::vir,
                    IntegralTransform::HalfTrans::ReadAndNuke);

// New first half for (OV| transforms
ints->transform_tei(MOSpace::occ, MOSpace::vir,
                    MOSpace::occ, MOSpace::vir,
                    IntegralTransform::HalfTrans::MakeAndKeep);

// Reuse for (OV|VV), then delete
ints->transform_tei(MOSpace::occ, MOSpace::vir,
                    MOSpace::vir, MOSpace::vir,
                    IntegralTransform::HalfTrans::ReadAndNuke);

HalfTrans Options:
  • MakeAndKeep - Compute first-half intermediates, keep for reuse
  • ReadAndKeep - Read existing first-half, keep for more reuse
  • MakeAndNuke - Compute first-half, delete after final transform
  • ReadAndNuke - Read first-half, delete after final transform (DEFAULT)


Pattern 3: Explicit First/Second Half Transforms
-------------------------------------------------

// When you need fine control or custom intermediate processing:

ints->transform_tei_first_half(MOSpace::occ, MOSpace::occ);
// ... potentially do something with intermediates ...
ints->transform_tei_second_half(MOSpace::occ, MOSpace::occ,
                                MOSpace::vir, MOSpace::vir);


Pattern 4: Unrestricted Transform
----------------------------------

std::vector<std::shared_ptr<MOSpace>> spaces;
spaces.push_back(MOSpace::occ);
spaces.push_back(MOSpace::vir);

auto ints = std::make_shared<IntegralTransform>(
    wfn, spaces,
    IntegralTransform::TransformationType::Unrestricted,
    IntegralTransform::OutputType::DPDOnly,
    IntegralTransform::MOOrdering::QTOrder,
    IntegralTransform::FrozenOrbitals::OccOnly,
    false
);

ints->initialize();
ints->update_orbitals();

// Transform all spin cases needed:
// Alpha-Alpha
ints->transform_tei(MOSpace::occ, MOSpace::occ,
                    MOSpace::occ, MOSpace::occ,
                    IntegralTransform::HalfTrans::MakeAndKeep);

// Alpha-Beta
ints->transform_tei(MOSpace::occ, MOSpace::occ,
                    MOSpace::occ, MOSpace::vir,
                    IntegralTransform::HalfTrans::ReadAndKeep);
// ... etc for other spin combinations


Pattern 5: Transformation for Orbital-Optimized Methods
-------------------------------------------------------

// For methods that update orbitals iteratively (e.g., OMP2, OCEPA)

auto ints = std::make_shared<IntegralTransform>(
    wfn, spaces,
    IntegralTransform::TransformationType::Restricted,
    IntegralTransform::OutputType::DPDOnly,
    IntegralTransform::MOOrdering::QTOrder,
    IntegralTransform::FrozenOrbitals::OccOnly,
    false
);

// Keep SO integrals to avoid re-reading from disk
ints->set_keep_dpd_so_ints(true);
ints->set_keep_iwl_so_ints(true);
ints->initialize();

// In each orbital optimization iteration:
for (int iter = 0; iter < max_iter; ++iter) {
    // ... optimize orbitals, update wavefunction ...

    ints->update_orbitals();  // Sync new orbitals
    ints->transform_tei(...); // Transform with new orbitals

    // ... compute energy and gradient ...
}


Pattern 6: Getting Frozen Core Information
------------------------------------------

// After any transformation that calls presort_so_tei():
double E_frozen_core = ints->get_frozen_core_energy();

// Frozen core operator is automatically written to PSIF_OEI
// and incorporated into correlation calculations that use libtrans.
// You typically don't need to manually load it.


Pattern 7: Custom DPD Buffer Names
-----------------------------------

// If you need specific buffer names (e.g., for multiple transformations)
ints->set_aa_int_name("My AA Ints");
ints->set_ab_int_name("My AB Ints");
ints->set_bb_int_name("My BB Ints");
ints->transform_tei(...);

// Access with custom name:
global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0,
                       ID("[O,O]"), ID("[V,V]"),
                       ID("[O>=O]+"), ID("[V>=V]+"), 0,
                       "My AA Ints");

================================================================================
7. PERFORMANCE OPTIMIZATION
================================================================================

Memory Management:
-----------------

1. Set appropriate memory limit:
   ints->set_memory(2000);  // 2000 MB

   Larger memory → fewer passes → faster transformation
   Default uses Psi4's global memory setting

2. Reuse half-transformed intermediates:
   Use ReadAndKeep when transforming multiple integrals with same first half
   Example: (OO|OO), (OO|OV), (OO|VV) share (OO| first half

3. SO integral caching:
   For iterative methods (orbital optimization):
   - set_keep_dpd_so_ints(true)  // Keep presorted SO integrals
   - set_keep_iwl_so_ints(true)  // Keep IWL SO integrals

   Trade-off: Uses more disk but avoids re-reading/presorting

4. Clean up when done:
   Use MakeAndNuke / ReadAndNuke for last use of a first-half
   Don't keep intermediates you won't reuse


Transform Ordering:
------------------

Order your transforms to maximize half-transform reuse:

GOOD:
  (OO|OO) MakeAndKeep
  (OO|OV) ReadAndKeep
  (OO|VV) ReadAndNuke  ← Last use of (OO|
  (OV|OV) MakeAndKeep
  (OV|VV) ReadAndNuke  ← Last use of (OV|
  (VV|VV) MakeAndNuke  ← Only use of (VV|

BAD:
  (OO|OO) MakeAndNuke  ← Wastes (OO| first half
  (OV|OV) MakeAndNuke  ← Wastes (OV| first half
  (OO|OV) MakeAndNuke  ← Recomputes (OO| first half!


Presort Optimization:
--------------------

If doing multiple transformations with same SO integrals:

ints->presort_so_tei();  // Do once
ints->set_tei_already_presorted(true);

// Now all transforms skip re-presorting
ints->transform_tei(...);
ints->transform_tei(...);


DPD vs IWL:
----------

Prefer DPD output (OutputType::DPDOnly):
  ✓ More efficient for most correlation methods
  ✓ Better memory access patterns
  ✓ Automatic symmetry handling
  ✓ Direct integration with libdpd

Use IWL only if:
  ✗ Legacy code requires it
  ✗ Specific algorithm demands sequential access

================================================================================
8. ADVANCED FEATURES
================================================================================

Backtransformation (Density Matrices):
--------------------------------------

For gradient calculations, you can backtransform density matrices:

ints->backtransform_density();
ints->backtransform_tpdm_restricted();   // or _unrestricted()


Custom PSIO Object:
------------------

To use a specific PSIO instance:

ints->set_psio(my_psio);  // Before initialize()


DPD Instance Management:
-----------------------

Multiple DPD instances can coexist:

auto ints1 = std::make_shared<IntegralTransform>(...);
ints1->set_dpd_id(0);

auto ints2 = std::make_shared<IntegralTransform>(...);
ints2->set_dpd_id(1);

// Switch between them:
dpd_set_default(ints1->get_dpd_id());
// ... work with ints1 buffers ...

dpd_set_default(ints2->get_dpd_id());
// ... work with ints2 buffers ...


File Number Configuration:
-------------------------

Default files are set automatically, but can be overridden:

ints->set_dpd_int_file(PSIF_LIBTRANS_DPD);  // DPD integrals file
ints->set_so_tei_file(PSIF_SO_TEI);         // SO TEI input file

================================================================================
9. FROZEN CORE HANDLING
================================================================================

The libtrans library provides CANONICAL frozen core handling for Psi4.
This is a key design feature: correlation modules can ignore frozen orbitals
entirely, with libtrans handling all bookkeeping.

What libtrans Does Automatically:
---------------------------------

1. Detects frozen orbitals from wavefunction object
   - Frozen occupied: frzcpi_
   - Frozen virtual: frzvpi_

2. Constructs frozen core operator
   Location: integraltransform_sort_so_tei.cc

   The frozen core operator is:
   F_core = H_core + J_core - K_core

   Where J_core and K_core are Coulomb and exchange contributions from
   frozen core electron density. This is written to PSIF_OEI as "MO-basis
   Alpha Fock matrix" and acts as an effective one-electron operator for
   correlation calculations.

3. Computes frozen core energy
   E_fzc = Σ_ij D_ij^core (H_ij + F_ij^core)

   This is the energy contribution from frozen core orbitals only:
   - Core kinetic energy
   - Core-nuclear attraction
   - Core-core repulsion

   Retrieved via: get_frozen_core_energy()

4. Excludes frozen orbitals from transformed integrals
   Only active orbitals appear in output integrals per FrozenOrbitals setting


Using Frozen Core Information:
------------------------------

// In your correlation method:
double E_fzc = ints->get_frozen_core_energy();
double E_correlation = /* your correlation energy */;
double E_total = E_ref + E_fzc + E_correlation;

// Frozen core operator is automatically used in:
// - One-electron integral transformations
// - Fock matrix construction
// You don't need to manually apply it


Frozen Orbital Options:
-----------------------

FrozenOrbitals::None
  - No orbitals excluded (all-electron calculation)
  - Use for benchmarking or when explicitly treating all electrons

FrozenOrbitals::OccOnly
  - Exclude frozen occupied only
  - Common for gradient calculations
  - Virtual space includes all virtuals (frozen + active)

FrozenOrbitals::VirOnly
  - Exclude frozen virtual only
  - Uncommon, for specialized calculations

FrozenOrbitals::OccAndVir (STANDARD)
  - Exclude both frozen occ and frozen vir
  - Standard for coupled cluster methods
  - Most efficient: smallest integral space


Checking Frozen Orbital Setup:
------------------------------

// After initialize(), check dimensions:
outfile->Printf("  Frozen core orbitals: %d\n", ints->nirrep());
// Use wavefunction accessors for details:
outfile->Printf("  Frozen core per irrep: %s\n",
                wfn->frzcpi().to_string().c_str());
outfile->Printf("  Frozen virt per irrep: %s\n",
                wfn->frzvpi().to_string().c_str());

================================================================================
10. COMMON PITFALLS AND SOLUTIONS
================================================================================

Pitfall 1: Forgetting update_orbitals()
---------------------------------------
SYMPTOM: Transformations use old/wrong orbitals

SOLUTION: Always call before transformation:
  ints->update_orbitals();
  ints->transform_tei(...);


Pitfall 2: Wrong Half-Transform Management
------------------------------------------
SYMPTOM: "Cannot read half-transformed integrals" error

SOLUTION: Match MakeAndKeep → ReadAndKeep → ... → ReadAndNuke
  First use: MakeAndKeep
  Middle uses: ReadAndKeep
  Last use: ReadAndNuke


Pitfall 3: Mixing DPD Instances
-------------------------------
SYMPTOM: Cannot find integrals, DPD errors

SOLUTION: Track DPD instance consistently:
  ints->set_dpd_id(0);
  ints->initialize();
  dpd_set_default(ints->get_dpd_id());  // Before DPD calls


Pitfall 4: Insufficient Memory
------------------------------
SYMPTOM: Very slow transformation with many passes

SOLUTION: Increase memory allocation:
  ints->set_memory(4000);  // MB


Pitfall 5: Wrong Orbital Ordering Assumption
--------------------------------------------
SYMPTOM: Integrals appear in unexpected order

SOLUTION: Remember:
  - QTOrder: By orbital class (occ, vir, etc.)
  - PitzerOrder: By irrep
  - Only affects IWL output
  - DPD output is space-based (use MOSpace to access)


Pitfall 6: Not Initializing Before Use
--------------------------------------
SYMPTOM: Segfault or "Object not initialized" error

SOLUTION: If initialize=false in constructor:
  auto ints = std::make_shared<IntegralTransform>(..., false);
  ints->set_print(1);      // Configure options
  ints->set_dpd_id(0);
  ints->initialize();      // NOW initialize


Pitfall 7: Transform Order Causes Recomputation
-----------------------------------------------
SYMPTOM: Unexpectedly slow, seeing repeated first-half transforms

SOLUTION: Group transforms by first-half:
  All (OO|XX) together
  All (OV|XX) together
  All (VV|XX) together


Pitfall 8: Expecting Frozen Core with No Frozen Orbitals
---------------------------------------------------------
SYMPTOM: get_frozen_core_energy() returns 0.0, but expected non-zero

SOLUTION: Check if wavefunction actually has frozen orbitals:
  if (wfn->frzcpi().sum() == 0) {
      // No frozen core orbitals defined
      // Add to input: set freeze_core true
  }


Pitfall 9: Using IWLAndDPD Unnecessarily
----------------------------------------
SYMPTOM: Running out of disk space

SOLUTION: Use DPDOnly unless you specifically need IWL:
  OutputType::DPDOnly  // Not IWLAndDPD


Pitfall 10: Not Cleaning Up Intermediates
-----------------------------------------
SYMPTOM: Disk fills up with intermediate files

SOLUTION: Use ReadAndNuke for last use of a first-half:
  ints->transform_tei(..., HalfTrans::ReadAndNuke);  // Last use

  Or set:
  ints->set_keep_ht_ints(false);  // Don't keep any intermediates

================================================================================
11. WHEN NOT TO USE LIBTRANS
================================================================================

libtrans is designed for standard single-reference correlation methods.
Consider alternatives if:

Density-Fitting Methods:
------------------------
Modern DF/CD methods use 3-index integrals (ia|Q) and specialized algorithms.

See instead:
• psi4/src/psi4/dfocc/ - DF orbital-optimized CC
• psi4/lib3index/ - 3-index integral handling

libtrans is designed for 4-index AO→MO transforms, not optimal for DF.


Multi-Reference Methods with Custom Blocking:
---------------------------------------------
MRCC methods need specialized integral storage and blocking strategies.

See instead:
• psi4/src/psi4/psimrcc/ - PSIMRCC transformation (CCTransform class)

PSIMRCC uses:
- Block-based memory management for large active spaces
- Custom integral maps: std::map<size_t, double>
- Direct integration with MRCC's CCIndex system


Performance-Critical Custom Needs:
----------------------------------
If profiling shows >10% time in transformations AND you have:
- Very specific access patterns
- Specialized symmetry exploitation
- Target-specific optimizations (GPU, etc.)

Then custom implementation may be justified. But:
1. Profile first with libtrans
2. Document why custom is necessary
3. Consider contributing optimizations back to libtrans


Half-Transformations for Specialized Algorithms:
------------------------------------------------
Some algorithms (e.g., AO-based RI-MP2, Laplace transform methods) need
MO-basis on one side, AO on the other.

libtrans provides first_half and second_half methods, but if you need
very specific intermediate formats, custom code may be clearer.

However, ccenergy/halftrans.cc, cclambda/halftrans.cc, and
dct/half_transform.cc show that similar needs might still benefit from
consolidation into libtrans. Evaluate carefully!

================================================================================

EXAMPLE: Complete Correlation Method Setup
===========================================

// Example: MP2 energy using libtrans

#include "psi4/libtrans/integraltransform.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/libdpd/dpd.h"

void compute_mp2_energy(SharedWavefunction ref_wfn) {
    // 1. Define orbital spaces
    std::vector<std::shared_ptr<MOSpace>> spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);

    // 2. Create transformation object
    auto ints = std::make_shared<IntegralTransform>(
        ref_wfn, spaces,
        IntegralTransform::TransformationType::Restricted,
        IntegralTransform::OutputType::DPDOnly,
        IntegralTransform::MOOrdering::QTOrder,
        IntegralTransform::FrozenOrbitals::OccAndVir,
        false  // Delay initialization
    );

    // 3. Configure options
    ints->set_print(1);
    ints->set_dpd_id(0);
    ints->set_keep_dpd_so_ints(false);
    ints->initialize();

    // 4. Transform integrals
    ints->update_orbitals();
    ints->transform_tei(MOSpace::occ, MOSpace::occ,
                        MOSpace::vir, MOSpace::vir);

    // 5. Get frozen core energy
    double E_fzc = ints->get_frozen_core_energy();

    // 6. Compute MP2 correlation energy
    dpd_set_default(ints->get_dpd_id());

    dpdbuf4 K, D;
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0,
                           ID("[O,O]"), ID("[V,V]"),
                           ID("[O>=O]+"), ID("[V>=V]+"), 0,
                           "MO Ints (OO|VV)");

    // ... compute MP2 energy using K integrals ...
    double E_mp2_corr = /* MP2 correlation computation */;

    global_dpd_->buf4_close(&K);

    // 7. Total energy
    double E_ref = ref_wfn->energy();
    double E_total = E_ref + E_fzc + E_mp2_corr;

    outfile->Printf("  Frozen core energy:   %20.14f\n", E_fzc);
    outfile->Printf("  MP2 correlation:      %20.14f\n", E_mp2_corr);
    outfile->Printf("  Total MP2 energy:     %20.14f\n", E_total);
}

================================================================================

For further information, consult:
• integraltransform.h - Full class interface
• occ/trans_ints_rhf.cc - Exemplary usage in OCC module
• occ/trans_ints_uhf.cc - Unrestricted example
• DPD documentation in libdpd/

Questions or improvements? Contact the Psi4 development team.

================================================================================
