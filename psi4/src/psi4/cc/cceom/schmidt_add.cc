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
    \ingroup CCEOM
    \brief Schmidt orthogonalization and vector addition for EOM-CC

    This file implements Schmidt orthogonalization for EOM-CC vectors using
    the Modified Gram-Schmidt algorithm. The implementation follows the same
    algorithmic pattern as the general OrthoUtil library (psi4/libqt/ortho_util.h)
    but is adapted for DPD (Distributed Paired Data) objects used in coupled
    cluster calculations.

    See also: psi4/libqt/ortho_util.h for the general-purpose orthogonalization
              utilities and eom_ortho_util.h for DPD-aware wrappers.
*/
#include <cstdio>
#include <cmath>
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"  // For general OrthoUtil reference
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#include "eom_ortho_util.h"  // DPD-aware orthogonalization utilities
#define EXTERN
#include "globals.h"

namespace psi {
namespace cceom {

/*!
** \brief Orthogonalize and add EOM vector to basis (ROHF/UHF)
**
** This function implements the Modified Gram-Schmidt algorithm for EOM-CC vectors:
** 1. Orthogonalize residual R against all existing basis vectors C[i]
** 2. Normalize the orthogonalized vector
** 3. Add to basis if norm exceeds threshold
**
** This is the DPD-object version of OrthoUtil::schmidt_add() for coupled cluster
** calculations. The algorithm follows the same pattern:
**   for i in range(num_basis):
**       overlap = <R|C[i]>
**       R = R - overlap * C[i]    # Remove component along C[i]
**   norm = ||R||
**   if norm > threshold:
**       R = R / norm               # Normalize
**       C[num_basis] = R          # Add to basis
**       num_basis++
**
** \param RIA Alpha singles amplitudes (modified in place)
** \param Ria Beta singles amplitudes (modified in place)
** \param RIJAB Alpha-alpha doubles amplitudes (modified in place)
** \param Rijab Beta-beta doubles amplitudes (modified in place)
** \param RIjAb Alpha-beta doubles amplitudes (modified in place)
** \param numCs Number of basis vectors (updated if vector is added)
** \param irrep Irreducible representation
*/

extern double norm_C(dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);

extern void scm_C(dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a);

/* use for ROHF and UHF */
void schmidt_add(dpdfile2 *RIA, dpdfile2 *Ria, dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb, int *numCs, int irrep) {
    double overlap;  // Overlap with basis vector
    double norm;
    int i;
    dpdfile2 Cme, CME;
    dpdbuf4 CMNEF, Cmnef, CMnEf;
    char CME_lbl[32], Cme_lbl[32], CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32];

    // Modified Gram-Schmidt: Orthogonalize R against all existing basis vectors
    // Algorithm: for each basis vector C[i], compute overlap and subtract projection
    for (i = 0; i < *numCs; i++) {
        sprintf(CME_lbl, "%s %d", "CME", i);
        sprintf(Cme_lbl, "%s %d", "Cme", i);
        sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
        sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
        sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);

        // Load basis vector C[i] from disk
        global_dpd_->file2_init(&CME, PSIF_EOM_CME, irrep, 0, 1, CME_lbl);
        global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, irrep, 2, 7, 2, 7, 0, CMNEF_lbl);
        if (params.eom_ref == 1) {
            global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, irrep, 0, 1, Cme_lbl);
            global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, irrep, 2, 7, 2, 7, 0, Cmnef_lbl);
            global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, irrep, 0, 5, 0, 5, 0, CMnEf_lbl);
        } else if (params.eom_ref == 2) {
            global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, irrep, 2, 3, Cme_lbl);
            global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, irrep, 12, 17, 12, 17, 0, Cmnef_lbl);
            global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, irrep, 22, 28, 22, 28, 0, CMnEf_lbl);
        }

        // Orthogonalize: Compute overlap and subtract projection
        // Uses wrapper from eom_ortho_util.h (analogous to OrthoUtil::orthogonalize_vector)
        orthogonalize_C(RIA, Ria, RIJAB, Rijab, RIjAb, &CME, &Cme, &CMNEF, &Cmnef, &CMnEf, &overlap);

        // outfile->Printf( "Overlap with basis vector %d = %20.14f\n", i, overlap);

        global_dpd_->file2_close(&CME);
        global_dpd_->file2_close(&Cme);
        global_dpd_->buf4_close(&CMNEF);
        global_dpd_->buf4_close(&Cmnef);
        global_dpd_->buf4_close(&CMnEf);
    }

    // Compute norm of orthogonalized vector (analogous to OrthoUtil::normalize_vector)
    norm = norm_C(RIA, Ria, RIJAB, Rijab, RIjAb);
    // outfile->Printf( "Norm of residual (TDC) = %20.14f\n", norm);

    // Threshold check: only add vector if norm exceeds tolerance
    // This prevents adding linearly dependent vectors to the basis
    if (norm < eom_params.schmidt_add_residual_tol) {
        return;  // Vector is (nearly) in span of existing basis, don't add
    } else {
        // Normalize the vector: R = R / ||R||
        scm_C(RIA, Ria, RIJAB, Rijab, RIjAb, 1.0 / norm);

        // Add normalized vector to basis as C[numCs]
        sprintf(CME_lbl, "%s %d", "CME", *numCs);
        sprintf(Cme_lbl, "%s %d", "Cme", *numCs);
        sprintf(CMNEF_lbl, "%s %d", "CMNEF", *numCs);
        sprintf(Cmnef_lbl, "%s %d", "Cmnef", *numCs);
        sprintf(CMnEf_lbl, "%s %d", "CMnEf", *numCs);

        global_dpd_->file2_copy(RIA, PSIF_EOM_CME, CME_lbl);
        global_dpd_->file2_copy(Ria, PSIF_EOM_Cme, Cme_lbl);
        global_dpd_->buf4_copy(RIJAB, PSIF_EOM_CMNEF, CMNEF_lbl);
        global_dpd_->buf4_copy(Rijab, PSIF_EOM_Cmnef, Cmnef_lbl);
        global_dpd_->buf4_copy(RIjAb, PSIF_EOM_CMnEf, CMnEf_lbl);

        ++(*numCs);  // Increment basis size
    }
    return;
}

/*!
** \brief Orthogonalize and add EOM vector to basis (RHF)
**
** RHF version of schmidt_add with proper spin adaptation.
** For RHF, the singles have a factor of 2 and doubles use 2*RIjAb - RIjbA.
**
** This is the DPD-object RHF version of OrthoUtil::schmidt_add().
**
** \param RIA Singles amplitudes (modified in place)
** \param RIjAb Doubles amplitudes (modified in place)
** \param numCs Number of basis vectors (updated if vector is added)
** \param irrep Irreducible representation
*/
void schmidt_add_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, int *numCs, int irrep) {
    double overlap, norm, R0, C0;
    int i;
    dpdfile2 CME;
    dpdbuf4 CMnEf, CAB1, CAB2;
    dpdfile2 R1;
    dpdbuf4 R2a, R2b;
    char CME_lbl[32], Cme_lbl[32], CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32], C0_lbl[32];

    if (params.full_matrix) psio_read_entry(PSIF_EOM_R, "R0", (char *)&R0, sizeof(double));

    // Modified Gram-Schmidt for RHF: orthogonalize against all basis vectors
    for (i = 0; i < *numCs; i++) {
        sprintf(CME_lbl, "%s %d", "CME", i);
        sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);

        // Load basis vector
        global_dpd_->file2_init(&CME, PSIF_EOM_CME, irrep, 0, 1, CME_lbl);
        global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, irrep, 0, 5, 0, 5, 0, CMnEf_lbl);

        // Orthogonalize using RHF-specific wrapper (handles spin adaptation)
        // Uses wrapper from eom_ortho_util.h (analogous to OrthoUtil::orthogonalize_vector)
        orthogonalize_C_RHF(RIA, RIjAb, &CME, &CMnEf, &overlap, irrep);

        // Handle full_matrix case (includes R0 component)
        if (params.full_matrix) {
            sprintf(C0_lbl, "%s %d", "C0", i);
            psio_read_entry(PSIF_EOM_CME, C0_lbl, (char *)&C0, sizeof(double));
            // Add R0*C0 contribution to overlap
            overlap += C0 * R0;
            // Orthogonalize R0 component
            R0 = R0 - 1.0 * overlap * C0;
        }

        // outfile->Printf( "Overlap with basis vector %d = %20.14f\n", i, overlap);

        global_dpd_->file2_close(&CME);
        global_dpd_->buf4_close(&CMnEf);
    }

    // Compute RHF norm with proper spin adaptation
    // For RHF: ||v||^2 = 2*<singles|singles> + 2*<RIjAb|RIjAb> - <RIjAb|RIjbA>
    global_dpd_->buf4_sort(RIjAb, PSIF_EOM_TMP, pqsr, 0, 5, "RIjbA");
    global_dpd_->buf4_init(&R2b, PSIF_EOM_TMP, irrep, 0, 5, 0, 5, 0, "RIjbA");

    norm = 2.0 * global_dpd_->file2_dot_self(RIA);
    norm += 2.0 * global_dpd_->buf4_dot_self(RIjAb);
    norm -= global_dpd_->buf4_dot(RIjAb, &R2b);
    if (params.full_matrix) norm += R0 * R0;
    norm = sqrt(norm);

    global_dpd_->buf4_close(&R2b);

    // outfile->Printf( "Norm of residual (TDC) = %20.14f\n", norm);

    // Threshold check: only add vector if norm exceeds tolerance
    if (norm < eom_params.schmidt_add_residual_tol) {
        return;  // Vector is (nearly) in span of existing basis, don't add
    } else {
        // Normalize the vector (analogous to OrthoUtil::normalize_vector)
        if (params.full_matrix) R0 *= 1.0 / norm;
        global_dpd_->file2_scm(RIA, 1.0 / norm);
        global_dpd_->buf4_scm(RIjAb, 1.0 / norm);

#ifdef EOM_DEBUG
        // Verify normalization
        global_dpd_->buf4_sort(RIjAb, PSIF_EOM_TMP, pqsr, 0, 5, "RIjbA");
        global_dpd_->buf4_init(&R2b, PSIF_EOM_TMP, irrep, 0, 5, 0, 5, 0, "RIjbA");
        norm = 2.0 * global_dpd_->file2_dot_self(RIA);
        norm += 2.0 * global_dpd_->buf4_dot_self(RIjAb);
        norm -= global_dpd_->buf4_dot(RIjAb, &R2b);
        if (params.full_matrix) norm += R0 * R0;
        norm = sqrt(norm);
        outfile->Printf("Norm of final new C in schmidt_add(): %20.15lf\n", norm);
        global_dpd_->buf4_close(&R2b);
#endif

        // Add normalized vector to basis as C[numCs]
        sprintf(CME_lbl, "%s %d", "CME", *numCs);
        sprintf(CMnEf_lbl, "%s %d", "CMnEf", *numCs);

        global_dpd_->file2_copy(RIA, PSIF_EOM_CME, CME_lbl);
        global_dpd_->buf4_copy(RIjAb, PSIF_EOM_CMnEf, CMnEf_lbl);

        // Generate AA and BB C2 vectors from AB vector for RHF
        // C(IJ,AB) = C(ij,ab) = C(Ij,Ab) - C(Ij,bA)
        global_dpd_->buf4_copy(RIjAb, PSIF_EOM_TMP, "CMnEf");
        global_dpd_->buf4_sort(RIjAb, PSIF_EOM_TMP, pqsr, 0, 5, "CMnfE");

        global_dpd_->buf4_init(&CAB1, PSIF_EOM_TMP, irrep, 0, 5, 0, 5, 0, "CMnEf");
        global_dpd_->buf4_init(&CAB2, PSIF_EOM_TMP, irrep, 0, 5, 0, 5, 0, "CMnfE");
        global_dpd_->buf4_axpy(&CAB2, &CAB1, -1.0);
        global_dpd_->buf4_close(&CAB2);
        global_dpd_->buf4_close(&CAB1);

        // Store R0 component if using full_matrix
        if (params.full_matrix) {
            sprintf(C0_lbl, "%s %d", "C0", *numCs);
            psio_write_entry(PSIF_EOM_CME, C0_lbl, (char *)&R0, sizeof(double));
        }
        ++(*numCs);  // Increment basis size
    }
    return;
}

}  // namespace cceom
}  // namespace psi
