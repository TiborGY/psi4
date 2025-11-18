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
** \file scf_S_inverse_sqrt.cc
** \brief Computation of S^{-1/2} and S^{1/2} matrices for MCSCF
** \ingroup MCSCF
**
** This file implements symmetric orthogonalization for the MCSCF module,
** computing both the inverse square root and square root of the overlap matrix.
**
** \par Algorithm: Symmetric Orthogonalization
**
** The symmetric orthogonalization algorithm computes:
** 1. Diagonalize overlap matrix: S = L λ L^T
** 2. Compute inverse square root: S^{-1/2} = L λ^{-1/2} L^T
** 3. Compute square root: S^{1/2} = L λ^{1/2} L^T
**
** This is mathematically equivalent to the symmetric orthogonalization
** method used in standard SCF calculations for transforming non-orthogonal
** AO basis functions to orthogonal basis functions.
**
** \par MCSCF-Specific Requirements
**
** The MCSCF module requires both matrices for different transformations:
** - **S^{-1/2}**: Transform operators H, F to orthogonal AO basis
** - **S^{1/2}**: Transform DIIS error matrices back from orthogonal AO basis
**
** Both matrices are computed once during initialization and stored for
** repeated use throughout MCSCF iterations.
**
** \par Implementation Notes
**
** This implementation uses MCSCF-specific data structures (SBlockMatrix)
** for integration with the MCSCF module infrastructure. It maintains a
** self-contained implementation rather than using BasisSetOrthogonalization
** due to type incompatibility and additional functional requirements
** (computing both S^{-1/2} and S^{1/2}).
**
** \see psi4::BasisSetOrthogonalization::compute_symmetric_orthog() for the
**      mathematically equivalent standard implementation using SharedMatrix
**      types (psi4/libmints/orthog.h)
** \see psi4/libmints/orthog.h for comprehensive basis set orthogonalization
**      with multiple algorithm options
*/

#include <cmath>
#include "scf.h"

namespace psi {
namespace mcscf {

/*!
 * \brief Construct S^{-1/2} and S^{1/2} matrices via symmetric orthogonalization
 *
 * This function computes the inverse square root and square root of the overlap
 * matrix S using symmetric eigendecomposition. The algorithm is:
 *
 * 1. Diagonalize: S = L λ L^T
 * 2. Compute: S^{-1/2} = L λ^{-1/2} L^T
 * 3. Compute: S^{1/2} = L λ^{1/2} L^T
 *
 * This is mathematically equivalent to symmetric orthogonalization used in
 * standard SCF calculations. The MCSCF module maintains its own implementation
 * using SBlockMatrix types for integration with MCSCF-specific infrastructure.
 *
 * Both matrices are stored (S_sqrt_inv and S_sqrt) for repeated use in:
 * - S^{-1/2}: Transforming operators to orthogonal AO basis (initial guess, iterations)
 * - S^{1/2}: Transforming operators back from orthogonal AO basis (DIIS)
 *
 * \par Mathematical Details
 *
 * Given overlap matrix S, the symmetric orthogonalization matrix X satisfies:
 * \f[
 *     X = S^{-1/2}, \quad X^T S X = I
 * \f]
 *
 * The eigendecomposition approach ensures X is symmetric:
 * \f[
 *     S = L \lambda L^T, \quad X = L \lambda^{-1/2} L^T
 * \f]
 *
 * This preserves maximum similarity between transformed and original basis
 * functions compared to canonical orthogonalization.
 *
 * \note This implementation does not check for linear dependencies. All eigenvalues
 *       are assumed to be non-zero. For ill-conditioned overlap matrices, consider
 *       using canonical orthogonalization with threshold filtering.
 *
 * \see psi4::BasisSetOrthogonalization for general basis orthogonalization
 *      using SharedMatrix types (psi4/libmints/orthog.h)
 * \see psi4::BasisSetOrthogonalization::compute_symmetric_orthog() for the
 *      standard implementation of this algorithm
 * \see psi4::BasisSetOrthogonalization::compute_canonical_orthog() for an
 *      alternative that handles linear dependencies
 */
void SCF::construct_S_inverse_sqrt() {
    SBlockVector lambda("lambda", nirreps, sopi);
    SBlockMatrix L("L", nirreps, sopi, sopi);
    SBlockMatrix Lambda("Lambda", nirreps, sopi, sopi);

    // Step 1: Diagonalize overlap matrix S = L λ L^T
    S.diagonalize(L, lambda);

    //   lambda->print();
    //   L->print();

    // Step 2: Compute S^{-1/2} = L λ^{-1/2} L^T
    for (int h = 0; h < nirreps; ++h) {
        for (int i = 0; i < sopi[h]; ++i) {
            Lambda->set(h, i, i, 1.0 / sqrt(lambda->get(h, i)));
        }
    }

    T.multiply(false, true, Lambda, L);  // T = λ^{-1/2} L^T
    S_sqrt_inv.multiply(false, false, L, T);  // S^{-1/2} = L T

    // Step 3: Compute S^{1/2} = L λ^{1/2} L^T
    for (int h = 0; h < nirreps; ++h) {
        for (int i = 0; i < sopi[h]; ++i) {
            Lambda->set(h, i, i, sqrt(lambda->get(h, i)));
        }
    }

    T.multiply(false, true, Lambda, L);  // T = λ^{1/2} L^T
    S_sqrt.multiply(false, false, L, T);  // S^{1/2} = L T
}

}  // namespace mcscf
}  // namespace psi
