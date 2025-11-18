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
** \file ortho_util.cc
** \brief Implementation of common orthogonalization utilities
** \ingroup QT
*/

#include "ortho_util.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include <cmath>
#include <cstring>

namespace psi {

void OrthoUtil::classical_gram_schmidt(double** A, int rows, int cols, std::string out_fname) {
    // Use the existing schmidt function
    schmidt(A, rows, cols, out_fname);
}

void OrthoUtil::modified_gram_schmidt(double** A, int rows, int cols) {
    double norm;

    // Modified Gram-Schmidt: orthogonalize rows
    // This algorithm is more numerically stable than Classical GS
    for (int k = 0; k < rows; k++) {
        // Normalize current row
        norm = 0.0;
        for (int i = 0; i < cols; i++) {
            norm += A[k][i] * A[k][i];
        }
        norm = sqrt(norm);

        if (norm > 1.0e-14) {  // Avoid division by zero
            for (int i = 0; i < cols; i++) {
                A[k][i] /= norm;
            }
        }

        // Orthogonalize subsequent rows against current row
        for (int j = k + 1; j < rows; j++) {
            double overlap = 0.0;
            for (int i = 0; i < cols; i++) {
                overlap += A[k][i] * A[j][i];
            }

            for (int i = 0; i < cols; i++) {
                A[j][i] -= overlap * A[k][i];
            }
        }
    }
}

void OrthoUtil::orthogonalize_vector(double* v, double** basis, int num_basis, int dim, double* overlaps) {
    // Orthogonalize vector v against all basis vectors
    for (int i = 0; i < num_basis; i++) {
        double overlap = dot_product(v, basis[i], dim);

        if (overlaps != nullptr) {
            overlaps[i] = overlap;
        }

        // v = v - overlap * basis[i]
        C_DAXPY(dim, -overlap, basis[i], 1, v, 1);
    }
}

double OrthoUtil::normalize_vector(double* v, int dim) {
    double norm = C_DNRM2(dim, v, 1);

    if (norm > 1.0e-14) {  // Avoid division by zero
        double inv_norm = 1.0 / norm;
        C_DSCAL(dim, inv_norm, v, 1);
    }

    return norm;
}

double OrthoUtil::dot_product(const double* v1, const double* v2, int dim) {
    return C_DDOT(dim, v1, 1, v2, 1);
}

bool OrthoUtil::schmidt_add(double* v, double** basis_list, int* num_basis, int max_basis, int dim, double threshold) {
    // Check if we have room for another vector
    if (*num_basis >= max_basis) {
        return false;
    }

    // Orthogonalize against existing basis
    orthogonalize_vector(v, basis_list, *num_basis, dim);

    // Normalize and check threshold
    double norm = normalize_vector(v, dim);

    // If norm is below threshold, don't add the vector
    if (norm < threshold) {
        return false;
    }

    // Add the vector to the basis
    memcpy(basis_list[*num_basis], v, dim * sizeof(double));
    (*num_basis)++;

    return true;
}

void OrthoUtil::modified_gram_schmidt_columns(double** A, int rows, int cols) {
    double norm;

    // Modified Gram-Schmidt for columns instead of rows
    for (int k = 0; k < cols; k++) {
        // Normalize current column
        norm = 0.0;
        for (int i = 0; i < rows; i++) {
            norm += A[i][k] * A[i][k];
        }
        norm = sqrt(norm);

        if (norm > 1.0e-14) {  // Avoid division by zero
            for (int i = 0; i < rows; i++) {
                A[i][k] /= norm;
            }
        }

        // Orthogonalize subsequent columns against current column
        for (int j = k + 1; j < cols; j++) {
            double overlap = 0.0;
            for (int i = 0; i < rows; i++) {
                overlap += A[i][k] * A[i][j];
            }

            for (int i = 0; i < rows; i++) {
                A[i][j] -= overlap * A[i][k];
            }
        }
    }
}

bool OrthoUtil::check_orthonormal(double** A, int num_vecs, int dim, double tolerance) {
    for (int i = 0; i < num_vecs; i++) {
        for (int j = i; j < num_vecs; j++) {
            double overlap = dot_product(A[i], A[j], dim);
            double expected = (i == j) ? 1.0 : 0.0;

            if (std::fabs(overlap - expected) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

}  // namespace psi
