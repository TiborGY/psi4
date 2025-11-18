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

/*! \file
    \ingroup CCLAMBDA
    \brief Orthogonalization of Lambda vectors against R vectors

    This file implements orthogonalization of Lambda (left eigenvector) amplitudes
    against previously converged R (right eigenvector) vectors in coupled cluster
    calculations. The algorithm follows the Modified Gram-Schmidt pattern similar to
    the general OrthoUtil library (psi4/libqt/ortho_util.h) and EOM-CC utilities
    (psi4/cc/cceom/eom_ortho_util.h).

    Key differences from standard orthogonalization:
    - Orthogonalizes Lambda vectors against R vectors (not L against L)
    - Uses special scaling factor: -overlap / (1.0 - R0*R0) for biorthogonal states
    - Works with pre-computed spin-adapted R vectors for efficiency

    See also: psi4/libqt/ortho_util.h for general orthogonalization utilities
              psi4/cc/cceom/eom_ortho_util.h for EOM-CC specific utilities
*/
#include <cstdio>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"  // For general OrthoUtil reference
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cclambda {

/*!
** \brief Compute dot product of Lambda and R vectors (RHF)
**
** Computes <L|R> for RHF Lambda and R vectors using spin-adapted R.
** This is conceptually similar to OrthoUtil::dot_product() but for
** biorthogonal CC vectors.
**
** For RHF: <L|R> = 2*<L1|R1> + <L2|2*R2-R2^T>
**
** \param IRR Irreducible representation
** \param R_index Index of R vector
** \return Overlap <L|R>
*/
double LRi_dot(int IRR, int R_index);

/*!
** \brief Orthogonalize Lambda vector against R vector
**
** Performs L = L - (overlap/(1-R0^2)) * R, removing the component of L
** along R. The special scaling factor accounts for the biorthogonal nature
** of left and right eigenvectors in EOM-CC theory.
**
** This is conceptually similar to OrthoUtil::orthogonalize_vector() but
** adapted for biorthogonal Lambda/R pairs.
**
** \param IRR Irreducible representation
** \param R_index Index of R vector to orthogonalize against
** \param overlap Pre-computed <L|R> overlap
** \param R0 R0 component of the R vector
*/
void LRi_minus(int IRR, int R_index, double overlap, double R0);

/*!
** \brief Orthogonalize Lambda vector against all R vectors
**
** This function implements Modified Gram-Schmidt orthogonalization of the
** current Lambda vector against all previously converged R vectors. This
** is essential for computing multiple excited-state Lambda vectors.
**
** Algorithm (similar to OrthoUtil::schmidt_add without the "add" part):
**   for each R vector in same irrep:
**       overlap = <L|R>
**       L = L - overlap * R  (with biorthogonal scaling)
**
** Currently supports RHF only (ROHF/UHF to be added).
**
** \param pL_params Array of Lambda parameters for all states
** \param current_L Index of current Lambda vector being computed
*/
void ortho_Rs(struct L_Params *pL_params, int current_L) {
    int L_root, L_irr;
    int R_root, R_irr;
    double overlap;
    int R;

    // Currently only RHF is supported
    if (params.ref != 0) return;

    L_irr = pL_params[current_L].irrep;
    L_root = pL_params[current_L].root;

    // Modified Gram-Schmidt: orthogonalize L against all R vectors in same irrep
    // Algorithm: for each R vector, compute overlap and subtract projection
    for (R = 1; R < params.nstates; ++R) {
        // Skip self (though L and R indices typically differ)
        if (R == current_L) continue;

        R_irr = pL_params[R].irrep;
        R_root = pL_params[R].root;

        // Only orthogonalize against R vectors in the same irrep
        if (L_irr != R_irr) continue;

        // Compute overlap <L|R> (analogous to OrthoUtil::dot_product)
        if (params.ref == 0) overlap = LRi_dot(L_irr, R_root);

        // Add R0 component if full_matrix calculation
        if (L_root == -1) overlap += pL_params[R].R0;

        // outfile->Printf("Overlap with R[%d][%d]: %15.10lf\n", R_irr, R_root, overlap);

        // Orthogonalize: L = L - (overlap/(1-R0^2)) * R
        // (analogous to OrthoUtil::orthogonalize_vector but with biorthogonal scaling)
        LRi_minus(L_irr, R_root, overlap, pL_params[R].R0);
    }
    return;
}

double LRi_dot(int IRR, int R_index) {
    dpdfile2 R1, L1;
    dpdbuf4 R2, L2;
    double overlap;
    char R1A_lbl[32], lbl[32];

    // Compute singles contribution: 2*<L1|R1>
    // Factor of 2 accounts for RHF spin adaptation (alpha and beta are equal)
    sprintf(R1A_lbl, "RIA %d %d", IRR, R_index);
    global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, IRR, 0, 1, R1A_lbl);
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, IRR, 0, 1, "New LIA");
    overlap = 2.0 * global_dpd_->file2_dot(&L1, &R1);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&L1);

    // Compute doubles contribution: <L2|2*R2 - R2^T>
    // Uses pre-computed spin-adapted R2 for efficiency
    sprintf(lbl, "2RIjAb - RIjbA %d %d", IRR, R_index);
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, IRR, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, IRR, 0, 5, 0, 5, 0, "New LIjAb");
    overlap += global_dpd_->buf4_dot(&L2, &R2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&R2);

    return overlap;
}

void LRi_minus(int IRR, int R_index, double overlap, double R0) {
    dpdfile2 R1, L1;
    dpdbuf4 R2, L2;
    char R1A_lbl[32], lbl[32];

    // Compute scaling factor for biorthogonal Lambda/R pairs
    // For biorthogonal eigenvectors: factor = -overlap / (1 - R0*R0)
    // This differs from standard orthogonalization which uses -overlap
    double scale_factor = -overlap / (1.0 - R0 * R0);

    // Orthogonalize singles: L1 = L1 + scale_factor * R1
    sprintf(R1A_lbl, "RIA %d %d", IRR, R_index);
    global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, IRR, 0, 1, R1A_lbl);
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, IRR, 0, 1, "New LIA");
    global_dpd_->file2_axpy(&R1, &L1, scale_factor, 0);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&L1);

    // Orthogonalize doubles: L2 = L2 + scale_factor * R2
    sprintf(lbl, "RIjAb %d %d", IRR, R_index);
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, IRR, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, IRR, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_axpy(&R2, &L2, scale_factor);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&R2);

    // For RHF, copy alpha Lambda to beta Lambda (spin symmetry)
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, IRR, 0, 1, "New LIA");
    global_dpd_->file2_copy(&L1, PSIF_CC_LAMBDA, "New Lia");
    global_dpd_->file2_close(&L1);

    return;
}

}  // namespace cclambda
}  // namespace psi
