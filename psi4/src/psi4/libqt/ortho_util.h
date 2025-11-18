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
** \file ortho_util.h
** \brief Common orthogonalization utility library for Psi4
** \ingroup QT
**
** This library provides a unified interface for various orthogonalization
** algorithms used throughout Psi4, including:
** - Classical Gram-Schmidt orthogonalization
** - Modified Gram-Schmidt orthogonalization
** - Schmidt orthogonalization with vector lists
**
** The utilities support both matrix-based and vector-based orthogonalization,
** making them suitable for use in CC, EOM, DETCI, and other modules.
*/

#pragma once

#include <string>

namespace psi {

// Forward declarations
class PSI_API OrthoUtil {
public:
    /**
     * @brief Classical Gram-Schmidt orthogonalization of matrix rows
     *
     * Orthogonalizes the rows of matrix A using the Classical Gram-Schmidt process.
     * This is a wrapper around the existing schmidt() function.
     *
     * @param A Matrix to orthogonalize (modified in place)
     * @param rows Number of rows in A
     * @param cols Number of columns in A
     * @param out_fname Output file name for debugging (default: "outfile")
     */
    static void classical_gram_schmidt(double** A, int rows, int cols,
                                       std::string out_fname = "outfile");

    /**
     * @brief Modified Gram-Schmidt orthogonalization of matrix rows
     *
     * Orthogonalizes the rows of matrix A using the Modified Gram-Schmidt process.
     * This method is numerically more stable than Classical GS for nearly-dependent vectors.
     *
     * @param A Matrix to orthogonalize (modified in place)
     * @param rows Number of rows in A
     * @param cols Number of columns in A
     */
    static void modified_gram_schmidt(double** A, int rows, int cols);

    /**
     * @brief Orthogonalize a vector against a list of orthonormal vectors
     *
     * Given a set of orthonormal vectors stored as rows in matrix basis,
     * orthogonalize vector v against all of them.
     *
     * @param v Vector to orthogonalize (modified in place)
     * @param basis Matrix containing orthonormal basis vectors as rows
     * @param num_basis Number of basis vectors
     * @param dim Dimension of vectors
     * @param overlaps Optional array to store overlap integrals (can be nullptr)
     */
    static void orthogonalize_vector(double* v, double** basis,
                                     int num_basis, int dim,
                                     double* overlaps = nullptr);

    /**
     * @brief Normalize a vector
     *
     * Normalizes vector v to unit length using the 2-norm.
     *
     * @param v Vector to normalize (modified in place)
     * @param dim Dimension of vector
     * @return Norm of the vector before normalization
     */
    static double normalize_vector(double* v, int dim);

    /**
     * @brief Compute dot product of two vectors
     *
     * @param v1 First vector
     * @param v2 Second vector
     * @param dim Dimension of vectors
     * @return Dot product v1 · v2
     */
    static double dot_product(const double* v1, const double* v2, int dim);

    /**
     * @brief Schmidt orthogonalization with threshold checking
     *
     * Orthogonalize vector v against a list of vectors and add it to the list
     * only if its norm after orthogonalization exceeds the threshold.
     *
     * @param v Vector to orthogonalize and potentially add
     * @param basis_list Matrix to store basis vectors (rows)
     * @param num_basis Current number of basis vectors (updated if vector is added)
     * @param max_basis Maximum number of basis vectors allowed
     * @param dim Dimension of vectors
     * @param threshold Minimum norm threshold for adding vector
     * @return True if vector was added to basis, false otherwise
     */
    static bool schmidt_add(double* v, double** basis_list, int* num_basis,
                           int max_basis, int dim, double threshold);

    /**
     * @brief Orthogonalize matrix columns (instead of rows)
     *
     * Performs Modified Gram-Schmidt orthogonalization on the columns of A.
     *
     * @param A Matrix to orthogonalize (modified in place)
     * @param rows Number of rows in A
     * @param cols Number of columns in A
     */
    static void modified_gram_schmidt_columns(double** A, int rows, int cols);

    /**
     * @brief Check orthonormality of a set of vectors
     *
     * Verify that the rows of matrix A form an orthonormal set.
     *
     * @param A Matrix containing vectors as rows
     * @param num_vecs Number of vectors
     * @param dim Dimension of vectors
     * @param tolerance Tolerance for orthonormality check
     * @return True if vectors are orthonormal within tolerance
     */
    static bool check_orthonormal(double** A, int num_vecs, int dim,
                                  double tolerance = 1.0e-10);
};

}  // namespace psi
