# Orthogonalization Utility Library

## Overview

The orthogonalization utility library (`ortho_util.h/cc`) provides a unified interface for various orthogonalization algorithms used throughout Psi4. This library consolidates orthogonalization functionality previously scattered across multiple modules (CC, EOM, DETCI, DFOCC) into a single, reusable API.

## Scope and Distinction from Basis Set Orthogonalization

**Important:** This library provides **vector orthogonalization for iterative subspace methods**. It is distinct from basis set orthogonalization used in SCF calculations.

### OrthoUtil (This Library)
- **Purpose:** Orthogonalize vectors using Gram-Schmidt algorithms
- **Use Case:** Davidson/Olsen eigensolvers, CI subspace expansion, CC equation solvers
- **Frequency:** Called thousands of times per calculation (performance-critical)
- **Data Structure:** Raw `double**` arrays for minimal overhead
- **Algorithm:** Incremental vector orthogonalization against growing basis sets
- **Example:** Orthogonalize new trial vector against existing Krylov subspace

### BasisSetOrthogonalization (psi4/libmints/orthog.h)
- **Purpose:** Transform atomic orbital basis to orthogonal molecular orbital basis
- **Use Case:** SCF initialization, basis set setup
- **Frequency:** Called once per calculation
- **Data Structure:** `SharedMatrix` objects (high-level Psi4 abstractions)
- **Algorithm:** Eigendecomposition of overlap matrix S, compute X = S^{-1/2}
- **Example:** Convert non-orthogonal AO basis to orthogonal MO basis for SCF

### Why Are They Separate?

While both involve "orthogonalization", they solve fundamentally different mathematical problems:

1. **Different Operations:**
   - OrthoUtil: v' = v - Σᵢ⟨vᵢ|v⟩vᵢ (Gram-Schmidt projection)
   - BasisSetOrthogonalization: X = U s^{-1/2} U^T where S = U s U^T (overlap matrix transformation)

2. **Different Abstraction Levels:**
   - OrthoUtil: Low-level numerical library (libqt)
   - BasisSetOrthogonalization: High-level quantum chemistry library (libmints)

3. **Different Performance Profiles:**
   - OrthoUtil: Hot path in iterative loops (must be highly optimized)
   - BasisSetOrthogonalization: One-time initialization (optimization less critical)

**For SCF and basis set work, use BasisSetOrthogonalization. For iterative solvers and subspace methods, use OrthoUtil.**

## Features

- **Classical Gram-Schmidt Orthogonalization**: Standard Gram-Schmidt process for orthogonalizing vectors
- **Modified Gram-Schmidt Orthogonalization**: Numerically stable variant for nearly-dependent vectors
- **Vector Orthogonalization**: Orthogonalize a vector against a list of basis vectors
- **Schmidt-Add**: Orthogonalize and conditionally add vectors to a basis with threshold checking
- **Utility Functions**: Normalization, dot products, and orthonormality checking

## Usage

### Including the Library

```cpp
#include "psi4/libqt/qt.h"
// ortho_util.h is automatically included via qt.h
```

### Classical Gram-Schmidt

```cpp
// Orthogonalize rows of a matrix
double** matrix = block_matrix(n_rows, n_cols);
// ... fill matrix ...
OrthoUtil::classical_gram_schmidt(matrix, n_rows, n_cols);
```

### Modified Gram-Schmidt (More Numerically Stable)

```cpp
// Orthogonalize rows of a matrix with better numerical stability
double** matrix = block_matrix(n_rows, n_cols);
// ... fill matrix ...
OrthoUtil::modified_gram_schmidt(matrix, n_rows, n_cols);
```

### Orthogonalize a Vector Against a Basis

```cpp
double* new_vector = new double[dim];
double** basis = block_matrix(num_basis_vecs, dim);
// ... initialize vectors ...

// Orthogonalize new_vector against all basis vectors
OrthoUtil::orthogonalize_vector(new_vector, basis, num_basis_vecs, dim);

// Normalize the result
double norm = OrthoUtil::normalize_vector(new_vector, dim);
```

### Schmidt-Add with Threshold

```cpp
double** basis_list = block_matrix(max_basis, dim);
int num_basis = 0;
double* candidate = new double[dim];
// ... initialize candidate vector ...

// Attempt to add vector to basis (only if norm > threshold after orthogonalization)
double threshold = 1.0e-6;
bool added = OrthoUtil::schmidt_add(candidate, basis_list, &num_basis,
                                    max_basis, dim, threshold);

if (added) {
    // Vector was added to basis
    std::cout << "Vector added. New basis size: " << num_basis << std::endl;
}
```

### Column-Based Orthogonalization

```cpp
// Orthogonalize columns instead of rows
double** matrix = block_matrix(n_rows, n_cols);
// ... fill matrix ...
OrthoUtil::modified_gram_schmidt_columns(matrix, n_rows, n_cols);
```

### Check Orthonormality

```cpp
double** vectors = block_matrix(num_vecs, dim);
// ... fill with supposedly orthonormal vectors ...

double tolerance = 1.0e-10;
bool is_orthonormal = OrthoUtil::check_orthonormal(vectors, num_vecs, dim, tolerance);

if (!is_orthonormal) {
    std::cout << "Warning: Vectors are not orthonormal!" << std::endl;
}
```

## Migration Guide

### From CC/CCEOM Module

**Before:**
```cpp
// In cceom/schmidt_add.cc
for (i = 0; i < *numCs; i++) {
    // ... manual orthogonalization with DPD objects ...
    dotval = global_dpd_->file2_dot(RIA, &CME);
    global_dpd_->file2_axpy(&CME, RIA, -1.0 * dotval, 0);
}
norm = norm_C(RIA, Ria, RIJAB, Rijab, RIjAb);
```

**After (for regular vectors):**
```cpp
#include "psi4/libqt/qt.h"

// Convert to regular vectors if possible, then:
OrthoUtil::orthogonalize_vector(v, basis, num_basis, dim);
double norm = OrthoUtil::normalize_vector(v, dim);
```

### From DFOCC Module

**Before:**
```cpp
// In dfocc/arrays.cc - Array2d::mgs()
for (int k = 0; k < dim1_; k++) {
    rmgs1 = 0.0;
    for (int i = 0; i < dim1_; i++) {
        rmgs1 += A2d_[i][k] * A2d_[i][k];
    }
    rmgs1 = sqrt(rmgs1);
    // ... manual orthogonalization ...
}
```

**After:**
```cpp
#include "psi4/libqt/qt.h"

// Use the common utility
OrthoUtil::modified_gram_schmidt(A2d_, dim1_, dim2_);
```

### From DETCI Module

**Before:**
```cpp
// In detci/civect.cc - CIvect::schmidt_add()
for (cvect = 0; cvect < L; cvect++) {
    tval = 0.0;
    // ... manual dot product and orthogonalization ...
}
```

**After:**
```cpp
#include "psi4/libqt/qt.h"

// For regular vector arrays:
bool success = OrthoUtil::schmidt_add(v, basis_list, &num_basis,
                                      max_basis, dim, threshold);
```

## API Reference

### OrthoUtil Class

All methods are static and can be called without instantiating the class.

#### classical_gram_schmidt
```cpp
static void classical_gram_schmidt(double** A, int rows, int cols,
                                   std::string out_fname = "outfile");
```
Orthogonalizes matrix rows using Classical Gram-Schmidt. Wrapper around existing `schmidt()`.

#### modified_gram_schmidt
```cpp
static void modified_gram_schmidt(double** A, int rows, int cols);
```
Orthogonalizes matrix rows using Modified Gram-Schmidt. More numerically stable than classical GS.

#### orthogonalize_vector
```cpp
static void orthogonalize_vector(double* v, double** basis, int num_basis,
                                int dim, double* overlaps = nullptr);
```
Orthogonalizes vector `v` against all vectors in `basis`. Optionally stores overlap integrals.

#### normalize_vector
```cpp
static double normalize_vector(double* v, int dim);
```
Normalizes vector to unit length. Returns the original norm.

#### dot_product
```cpp
static double dot_product(const double* v1, const double* v2, int dim);
```
Computes dot product of two vectors using BLAS.

#### schmidt_add
```cpp
static bool schmidt_add(double* v, double** basis_list, int* num_basis,
                       int max_basis, int dim, double threshold);
```
Orthogonalizes `v` and adds to basis if norm > threshold. Returns true if added.

#### modified_gram_schmidt_columns
```cpp
static void modified_gram_schmidt_columns(double** A, int rows, int cols);
```
Orthogonalizes matrix columns instead of rows.

#### check_orthonormal
```cpp
static bool check_orthonormal(double** A, int num_vecs, int dim,
                             double tolerance = 1.0e-10);
```
Verifies orthonormality of vector set within tolerance.

## Algorithm Details

### Classical Gram-Schmidt
- For each vector v_i:
  1. Normalize v_i
  2. For each subsequent vector v_j (j > i):
     - Project v_j onto v_i
     - Subtract projection from v_j

### Modified Gram-Schmidt
- For each vector v_i:
  1. Normalize v_i
  2. For each subsequent vector v_j (j > i):
     - Project v_j onto v_i
     - Immediately subtract projection from v_j
- More stable because projections are computed on already-orthogonalized vectors

## Performance Considerations

- All functions use BLAS operations (C_DDOT, C_DAXPY, C_DNRM2) for performance
- Modified Gram-Schmidt has the same computational complexity as Classical GS (O(n²m) for n vectors of dimension m)
- For large-scale problems, consider block orthogonalization or QR decomposition

## Testing

Unit tests for the orthogonalization utilities can be found in the libqt test suite.

## References

1. Golub, G. H.; Van Loan, C. F. "Matrix Computations", 4th Edition, Johns Hopkins, 2013.
2. Sherrill, C. D. "Notes on Classical and Modified Gram-Schmidt", 1994.

## Contributing

When adding new orthogonalization methods, please:
1. Add appropriate documentation in the header file
2. Include usage examples in this README
3. Add unit tests
4. Update migration guides for affected modules
