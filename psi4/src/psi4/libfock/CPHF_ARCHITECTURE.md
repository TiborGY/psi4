# CPHF/CPKS Architecture Documentation

## Overview

This document describes the consolidated CPHF (Coupled-Perturbed Hartree-Fock) and CPKS (Coupled-Perturbed Kohn-Sham) architecture in Psi4. The consolidation reduces code duplication across multiple response theory implementations while maintaining backward compatibility.

## Motivation

Prior to consolidation, Psi4 contained 4 independent CPHF/CPKS solver implementations:

1. **RCPHF** (`libfock/apps.h:141`) - General-purpose restricted CPHF
2. **CPHFRHamiltonian** (`libfock/hamiltonian.h:142`) - Hamiltonian for RCPHF
3. **CPHF_FISAPT** (`fisapt/fisapt.h:213`) - F-SAPT two-monomer CPHF
4. **CPKS_USAPT0** (`libsapt_solver/usapt0.h:338`) - Unrestricted SAPT0 CPKS

These implementations shared 80-90% of their code, particularly:
- Preconditioned Conjugate Gradient (PCG) iterative solver
- Jacobi preconditioning with orbital energy differences
- JK-based Hamiltonian-vector products
- Convergence checking

The consolidation extracts shared logic into reusable Hamiltonian classes while preserving specialized behavior.

## Architecture Design

### Class Hierarchy

```
Hamiltonian (base)
  └── RHamiltonian (restricted base)
        ├── CPHFRHamiltonian (general CPHF, configurable JK scaling)
        ├── CPHFFISAPTHamiltonian (F-SAPT two-monomer)
        └── CPKSUSAPTHamiltonian (USAPT unrestricted)
```

### Separation of Concerns

The architecture separates two orthogonal concerns:

1. **What to compute** (Hamiltonian classes)
   - Defines the response operator: `H·x` product
   - Orbital spaces (occupied/virtual)
   - JK scaling factors
   - Diagonal approximation for preconditioning

2. **How to compute** (Solver classes)
   - PCG iteration algorithm
   - Convergence criteria
   - Memory management
   - Output formatting

### Key Interfaces

#### Base Hamiltonian Interface

```cpp
class RHamiltonian : public Hamiltonian {
public:
    // Returns diagonal for preconditioning/guess
    virtual std::shared_ptr<Vector> diagonal() = 0;

    // Forms Hamiltonian-vector product: b = H·x
    virtual void product(const std::vector<std::shared_ptr<Vector>>& x,
                        std::vector<std::shared_ptr<Vector>>& b) = 0;
};
```

#### CPHFRHamiltonian - General CPHF

The base class for restricted CPHF problems with configurable JK scaling.

**Product formula:**
```
b = j_scale*J - k_scale_1*K - k_scale_2*K^T + (eps_vir - eps_occ)*x
```

**Default scaling:** `j_scale=4.0, k_scale_1=1.0, k_scale_2=1.0` (RHF CPHF)

**Key features:**
- Symmetry support via irreps
- Pack/unpack between Matrix and Vector formats
- Configurable via `set_jk_scales()`

**Example:**
```cpp
auto H = std::make_shared<CPHFRHamiltonian>(
    jk, Caocc, Cavir, eps_aocc, eps_avir);

// Customize for different formalism (e.g., UHF-style scaling)
H->set_jk_scales(2.0, 1.0, 1.0);

auto solver = CGRSolver::build_solver(options, H);
solver->solve();
```

#### CPHFFISAPTHamiltonian - Two-Monomer CPHF

Specialized for F-SAPT calculations with two interacting monomers.

**Problem structure:**
- Two independent monomers A and B
- Each responds to perturbation from the other
- Map-based interface: keys "A" and "B"

**Product formula:** Same as CPHFRHamiltonian (4J - K - K^T + diagonal)

**Example:**
```cpp
auto H = std::make_shared<CPHFFISAPTHamiltonian>(
    jk,
    Cocc_A, Cvir_A, eps_occ_A, eps_vir_A,
    Cocc_B, Cvir_B, eps_occ_B, eps_vir_B);

// Prepare perturbations
std::map<std::string, SharedMatrix> b;
b["A"] = perturbation_on_A;  // [no_A x nv_A]
b["B"] = perturbation_on_B;  // [no_B x nv_B]

// Compute products
auto products = H->product_map(b);
// products["A"] contains H_A·x_A
// products["B"] contains H_B·x_B
```

#### CPKSUSAPTHamiltonian - Unrestricted Two-Monomer CPKS

Specialized for unrestricted SAPT0 with spin-resolved responses.

**Problem structure:**
- Two monomers (A and B)
- Each has alpha and beta spins
- Map-based interface: keys "Aa", "Ab", "Ba", "Bb"

**Product formula:**
```
b = 2*J_total - K - K^T + (eps_vir - eps_occ)*x
where J_total = J_alpha + J_beta
```

**Key difference:** `2J` scaling (vs `4J`) due to unrestricted formalism

**Example:**
```cpp
auto H = std::make_shared<CPKSUSAPTHamiltonian>(
    jk,
    Cocca_A, Cvira_A, eps_occa_A, eps_vira_A,
    Coccb_A, Cvirb_A, eps_occb_A, eps_virb_A,
    Cocca_B, Cvira_B, eps_occa_B, eps_vira_B,
    Coccb_B, Cvirb_B, eps_occb_B, eps_virb_B);

// Prepare perturbations for all spin combinations
std::map<std::string, SharedMatrix> b;
b["Aa"] = perturbation_A_alpha;  // [noa_A x nva_A]
b["Ab"] = perturbation_A_beta;   // [nob_A x nvb_A]
b["Ba"] = perturbation_B_alpha;  // [noa_B x nva_B]
b["Bb"] = perturbation_B_beta;   // [nob_B x nvb_B]

// Compute products
auto products = H->product_map(b);
```

## Implementation Details

### JK Scaling Rationale

Different electronic structure theories use different JK scaling factors:

| Theory | J Scale | K Scale | Formula |
|--------|---------|---------|---------|
| RHF CPHF | 4.0 | 1.0 | 4J - K - K^T |
| UHF CPKS | 2.0 | 1.0 | 2J - K - K^T |
| Custom | configurable | configurable | j*J - k₁*K - k₂*K^T |

The base `CPHFRHamiltonian` class supports all via `set_jk_scales()`.

### Preconditioning

All implementations use Jacobi preconditioning with orbital energy differences:

```cpp
z[i][a] = r[i][a] / (eps_vir[a] - eps_occ[i])
```

This is the diagonal approximation to the full CPHF Hessian and provides excellent convergence acceleration.

### Memory Efficiency

The map-based interface (`product_map()`) allows computing only required products:
- If only monomer A needs updating, pass only `b["A"]`
- JK object computes only necessary J/K matrices
- Reduces memory footprint and computation time

### Backward Compatibility

Original solver classes (`CPHF_FISAPT`, `CPKS_USAPT0`) maintain their public interfaces:
- `compute_cphf()` / `compute_cpks()` methods unchanged
- Input/output data structures unchanged
- Internal implementation delegates to new Hamiltonians
- Numerical results bit-identical (same algorithms)

## Creating Custom Hamiltonians

To create a new CPHF variant:

### 1. Derive from CPHFRHamiltonian

For standard CPHF-like problems:

```cpp
class MyCustomHamiltonian : public CPHFRHamiltonian {
public:
    MyCustomHamiltonian(std::shared_ptr<JK> jk,
                        SharedMatrix Caocc, SharedMatrix Cavir,
                        std::shared_ptr<Vector> eps_aocc,
                        std::shared_ptr<Vector> eps_avir)
        : CPHFRHamiltonian(jk, Caocc, Cavir, eps_aocc, eps_avir) {
        // Set custom scaling
        set_jk_scales(3.0, 0.5, 0.5);  // Example: 3J - 0.5K - 0.5K^T
    }

    void print_header() const override {
        outfile->Printf("  ==> My Custom CPHF <==\n\n");
    }
};
```

### 2. Derive from RHamiltonian

For completely custom products:

```cpp
class CompletelyCustomHamiltonian : public RHamiltonian {
protected:
    // Your custom data members
    SharedMatrix my_custom_matrix_;

public:
    CompletelyCustomHamiltonian(std::shared_ptr<JK> jk, /* ... */)
        : RHamiltonian(jk) { /* ... */ }

    std::shared_ptr<Vector> diagonal() override {
        // Return diagonal approximation for preconditioning
        // Typically: eps_vir - eps_occ
    }

    void product(const std::vector<std::shared_ptr<Vector>>& x,
                 std::vector<std::shared_ptr<Vector>>& b) override {
        // Implement custom H·x product
        // Can use jk_->compute() for Coulomb/exchange
    }

    void print_header() const override {
        outfile->Printf("  ==> Completely Custom Hamiltonian <==\n\n");
    }
};
```

### 3. Use with CGRSolver

All Hamiltonians work with the existing solver infrastructure:

```cpp
auto solver = CGRSolver::build_solver(options, my_hamiltonian);
solver->b().push_back(perturbation_vector);
solver->initialize();
solver->solve();
auto solution = solver->x();
```

## Design Decisions

### Why Map-Based Interface for FISAPT/USAPT?

The map interface provides:
1. **Flexibility:** Easy to handle variable number of response vectors
2. **Clarity:** Keys ("A", "B", "Aa", etc.) self-document the problem
3. **Efficiency:** Compute only requested products
4. **Extensibility:** Easy to add more monomers/spins in future

### Why Keep Original Solver Classes?

Maintaining `CPHF_FISAPT` and `CPKS_USAPT0`:
1. **Backward compatibility:** Existing code continues to work
2. **Minimal risk:** PCG loop unchanged reduces integration risk
3. **Gradual migration:** Can test Hamiltonians independently
4. **Clear interfaces:** Users familiar with existing APIs

Future work could fully migrate to CGRSolver, but current approach balances innovation with stability.

### Why Not Merge All into One Class?

Separation provides:
1. **Single Responsibility:** Each class handles one problem type
2. **Testability:** Unit test each Hamiltonian independently
3. **Maintainability:** Changes to FISAPT don't affect USAPT
4. **Performance:** No runtime overhead from unused features
5. **Clarity:** Clear which Hamiltonian for which problem

## Testing Strategy

### Unit Tests

Each Hamiltonian should have tests verifying:
1. **Diagonal correctness:** `diagonal()` returns `eps_vir - eps_occ`
2. **Product correctness:** `product()` matches analytical formula
3. **Scaling:** JK scaling factors applied correctly
4. **Symmetry:** Results respect point group symmetry (if applicable)

### Integration Tests

Solver tests should verify:
1. **Convergence:** Solutions converge to specified tolerance
2. **Correctness:** Compare against known reference results
3. **Performance:** Convergence rate comparable to original

### Regression Tests

Existing FISAPT and SAPT0 tests automatically verify:
1. **Numerical equivalence:** Energies match to 1e-10 Hartree
2. **Convergence behavior:** Iterations match ±2 of original
3. **No regressions:** All test cases pass

## Performance Considerations

### Computational Complexity

All implementations have identical computational complexity:
- **Dominant cost:** JK builds - `O(N³)` to `O(N²)` depending on algorithm
- **Preconditioner:** `O(N)` - negligible
- **Other operations:** `O(N²)` - matrix multiplications

### Memory Usage

Memory requirements:
- **JK object:** Dominates (integral storage)
- **Response vectors:** `O(n_occ × n_vir)` per perturbation
- **Temporary matrices:** `O(n_basis²)` for products
- **Overhead from abstraction:** Negligible (few pointers)

### Optimization Opportunities

The new architecture enables:
1. **Shared JK object:** Reuse across multiple CPHF calls
2. **Batched products:** Compute multiple perturbations together
3. **Lazy evaluation:** Skip inactive spin components
4. **Template specialization:** Compile-time optimization for specific cases

## Migration Guide

### For Developers Using RCPHF

No changes needed - existing code continues to work.

Optional: Customize JK scaling if needed:
```cpp
auto H = std::make_shared<CPHFRHamiltonian>(/*...*/);
H->set_jk_scales(2.0, 1.0, 1.0);  // UHF-like scaling
```

### For Developers Using CPHF_FISAPT

No changes needed - interface unchanged.

Internal: Now uses `CPHFFISAPTHamiltonian` for products.

### For Developers Using CPKS_USAPT0

No changes needed - interface unchanged.

Internal: Now uses `CPKSUSAPTHamiltonian` for products.

### For New CPHF Implementations

Instead of copying an existing solver:
1. Create Hamiltonian class (see "Creating Custom Hamiltonians")
2. Use existing `CGRSolver` for iteration
3. Test against analytical results

## Future Enhancements

Potential improvements:

### 1. Full CGRSolver Migration
Migrate FISAPT/USAPT to use `CGRSolver` directly instead of custom PCG loops.

**Benefits:**
- Further code reduction
- Access to advanced solver features (DIIS, line search, etc.)
- Unified convergence behavior

**Challenges:**
- More invasive changes
- Need to adapt map interface to vector interface
- Testing burden

### 2. GPU Acceleration
Current JK objects support GPU. Future work:
- Ensure Hamiltonians work with GPU-accelerated JK
- Batch multiple perturbations for better GPU utilization
- Consider GPU-native matrix operations

### 3. Parallel Solvers
Enable parallel CPHF for multiple perturbations:
- Solve independent perturbations concurrently
- Requires thread-safe Hamiltonian classes
- Significant speedup for polarizability tensors

### 4. Additional Hamiltonians
Extend framework to:
- Time-dependent response (TDHF/TDDFT)
- Orbital-optimized methods (OO-MP2, Brueckner)
- Extended systems (periodic CPHF)

## References

### Code Locations

- **Base classes:** `psi4/src/psi4/libfock/hamiltonian.{h,cc}`
- **FISAPT Hamiltonian:** `psi4/src/psi4/libfock/hamiltonian_fisapt.{h,cc}`
- **USAPT Hamiltonian:** `psi4/src/psi4/libfock/hamiltonian_usapt.{h,cc}`
- **Solver infrastructure:** `psi4/src/psi4/libfock/solver.{h,cc}`
- **RCPHF application:** `psi4/src/psi4/libfock/apps.{h,cc}`

### Related Documentation

- **Psi4 Developer Guide:** General coding standards
- **Response Theory:** Mathematical background for CPHF/CPKS
- **JK Object Documentation:** Using Coulomb and exchange builders

### Commit History

- **Phase 1 & 2:** Extended base class, created FISAPT Hamiltonian, refactored CPHF_FISAPT
- **Phase 3 & 4:** Created USAPT Hamiltonian, refactored CPKS_USAPT0

### Contact

For questions or contributions, see Psi4 development documentation or file an issue on the Psi4 GitHub repository.

---

*Last updated: 2025-01-18*
*Architecture version: 1.0 (initial consolidation)*
