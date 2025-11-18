/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*!
** \file ortho_util_example.cc
** \brief Example usage of the orthogonalization utility library
** \ingroup QT
**
** This file demonstrates various use cases for the OrthoUtil class.
** It is not compiled into the library but serves as documentation.
*/

#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/psi4-dec.h"
#include <iostream>

namespace psi {

// Example 1: Basic Modified Gram-Schmidt orthogonalization
void example_modified_gram_schmidt() {
    outfile->Printf("\n==> Example 1: Modified Gram-Schmidt <==\n\n");

    int n_vecs = 3;
    int dim = 3;

    // Create a matrix with non-orthogonal vectors
    double** A = block_matrix(n_vecs, dim);
    A[0][0] = 1.0; A[0][1] = 0.0; A[0][2] = 1.0;
    A[1][0] = 1.0; A[1][1] = -1.0; A[1][2] = 0.0;
    A[2][0] = 0.0; A[2][1] = 0.0; A[2][2] = 0.5;

    outfile->Printf("Original matrix:\n");
    print_mat(A, n_vecs, dim, "outfile");

    // Orthogonalize using Modified Gram-Schmidt
    OrthoUtil::modified_gram_schmidt(A, n_vecs, dim);

    outfile->Printf("\nOrthogonalized matrix (MGS):\n");
    print_mat(A, n_vecs, dim, "outfile");

    // Verify orthonormality
    bool is_orthonormal = OrthoUtil::check_orthonormal(A, n_vecs, dim);
    outfile->Printf("\nOrthonormality check: %s\n",
                   is_orthonormal ? "PASSED" : "FAILED");

    free_block(A);
}

// Example 2: Orthogonalizing a vector against a basis
void example_vector_orthogonalization() {
    outfile->Printf("\n==> Example 2: Vector Orthogonalization <==\n\n");

    int dim = 4;
    int n_basis = 2;

    // Create orthonormal basis vectors
    double** basis = block_matrix(n_basis, dim);
    basis[0][0] = 1.0; basis[0][1] = 0.0; basis[0][2] = 0.0; basis[0][3] = 0.0;
    basis[1][0] = 0.0; basis[1][1] = 1.0; basis[1][2] = 0.0; basis[1][3] = 0.0;

    // Create a new vector to orthogonalize
    double* new_vec = new double[dim];
    new_vec[0] = 1.0; new_vec[1] = 1.0; new_vec[2] = 1.0; new_vec[3] = 0.0;

    outfile->Printf("Original vector: [%.3f, %.3f, %.3f, %.3f]\n",
                   new_vec[0], new_vec[1], new_vec[2], new_vec[3]);

    // Orthogonalize against basis
    double* overlaps = new double[n_basis];
    OrthoUtil::orthogonalize_vector(new_vec, basis, n_basis, dim, overlaps);

    outfile->Printf("After orthogonalization: [%.3f, %.3f, %.3f, %.3f]\n",
                   new_vec[0], new_vec[1], new_vec[2], new_vec[3]);
    outfile->Printf("Overlaps with basis: [%.3f, %.3f]\n",
                   overlaps[0], overlaps[1]);

    // Normalize
    double norm = OrthoUtil::normalize_vector(new_vec, dim);
    outfile->Printf("Norm before normalization: %.6f\n", norm);
    outfile->Printf("Normalized vector: [%.3f, %.3f, %.3f, %.3f]\n",
                   new_vec[0], new_vec[1], new_vec[2], new_vec[3]);

    delete[] new_vec;
    delete[] overlaps;
    free_block(basis);
}

// Example 3: Schmidt-add with threshold
void example_schmidt_add() {
    outfile->Printf("\n==> Example 3: Schmidt-Add with Threshold <==\n\n");

    int dim = 3;
    int max_basis = 5;
    double threshold = 1.0e-6;

    double** basis_list = block_matrix(max_basis, dim);
    int num_basis = 0;

    // First vector
    double* v1 = new double[dim];
    v1[0] = 1.0; v1[1] = 0.0; v1[2] = 0.0;
    bool added1 = OrthoUtil::schmidt_add(v1, basis_list, &num_basis,
                                        max_basis, dim, threshold);
    outfile->Printf("Vector 1 [1, 0, 0] added: %s, Basis size: %d\n",
                   added1 ? "YES" : "NO", num_basis);

    // Second vector (orthogonal to first)
    double* v2 = new double[dim];
    v2[0] = 0.0; v2[1] = 1.0; v2[2] = 0.0;
    bool added2 = OrthoUtil::schmidt_add(v2, basis_list, &num_basis,
                                        max_basis, dim, threshold);
    outfile->Printf("Vector 2 [0, 1, 0] added: %s, Basis size: %d\n",
                   added2 ? "YES" : "NO", num_basis);

    // Third vector (linearly dependent on first two - should be rejected)
    double* v3 = new double[dim];
    v3[0] = 1.0; v3[1] = 1.0; v3[2] = 0.0; // Linear combination of v1 and v2
    bool added3 = OrthoUtil::schmidt_add(v3, basis_list, &num_basis,
                                        max_basis, dim, threshold);
    outfile->Printf("Vector 3 [1, 1, 0] added: %s, Basis size: %d\n",
                   added3 ? "YES" : "NO", num_basis);
    outfile->Printf("  (Expected: NO, because it's in span of basis)\n");

    // Fourth vector (has component orthogonal to span)
    double* v4 = new double[dim];
    v4[0] = 0.0; v4[1] = 0.0; v4[2] = 1.0;
    bool added4 = OrthoUtil::schmidt_add(v4, basis_list, &num_basis,
                                        max_basis, dim, threshold);
    outfile->Printf("Vector 4 [0, 0, 1] added: %s, Basis size: %d\n",
                   added4 ? "YES" : "NO", num_basis);

    outfile->Printf("\nFinal basis:\n");
    print_mat(basis_list, num_basis, dim, "outfile");

    delete[] v1;
    delete[] v2;
    delete[] v3;
    delete[] v4;
    free_block(basis_list);
}

// Example 4: Comparison of Classical vs Modified Gram-Schmidt
void example_compare_classical_vs_modified() {
    outfile->Printf("\n==> Example 4: Classical vs Modified GS <==\n\n");

    int n_vecs = 3;
    int dim = 3;

    // Create two identical matrices
    double** A_classical = block_matrix(n_vecs, dim);
    double** A_modified = block_matrix(n_vecs, dim);

    // Fill with the same values
    for (int i = 0; i < n_vecs; i++) {
        for (int j = 0; j < dim; j++) {
            double val = (i == j) ? 1.0 : 0.1;  // Nearly diagonal
            A_classical[i][j] = val;
            A_modified[i][j] = val;
        }
    }

    outfile->Printf("Original matrix:\n");
    print_mat(A_classical, n_vecs, dim, "outfile");

    // Apply Classical GS
    OrthoUtil::classical_gram_schmidt(A_classical, n_vecs, dim);
    outfile->Printf("\nAfter Classical Gram-Schmidt:\n");
    print_mat(A_classical, n_vecs, dim, "outfile");
    bool classical_orthonormal = OrthoUtil::check_orthonormal(A_classical, n_vecs, dim);

    // Apply Modified GS
    OrthoUtil::modified_gram_schmidt(A_modified, n_vecs, dim);
    outfile->Printf("\nAfter Modified Gram-Schmidt:\n");
    print_mat(A_modified, n_vecs, dim, "outfile");
    bool modified_orthonormal = OrthoUtil::check_orthonormal(A_modified, n_vecs, dim);

    outfile->Printf("\nClassical GS orthonormality: %s\n",
                   classical_orthonormal ? "PASSED" : "FAILED");
    outfile->Printf("Modified GS orthonormality: %s\n",
                   modified_orthonormal ? "PASSED" : "FAILED");

    free_block(A_classical);
    free_block(A_modified);
}

// Main example driver (not actually called in normal Psi4 execution)
void run_ortho_util_examples() {
    outfile->Printf("\n");
    outfile->Printf("  =======================================================\n");
    outfile->Printf("  Orthogonalization Utility Library Examples\n");
    outfile->Printf("  =======================================================\n");

    example_modified_gram_schmidt();
    example_vector_orthogonalization();
    example_schmidt_add();
    example_compare_classical_vs_modified();

    outfile->Printf("\n");
    outfile->Printf("  =======================================================\n");
    outfile->Printf("  Examples completed successfully!\n");
    outfile->Printf("  =======================================================\n");
    outfile->Printf("\n");
}

}  // namespace psi
