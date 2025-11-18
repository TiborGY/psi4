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
** \file eom_ortho_util.cc
** \brief Implementation of DPD-aware orthogonalization utilities for EOM-CC
** \ingroup CCEOM
*/

#include "eom_ortho_util.h"
#include "psi4/libdpd/dpd.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cceom {

double dot_product_C(dpdfile2 *RIA, dpdfile2 *Ria, dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb,
                     dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf) {
    double dotval = 0.0;

    // Compute dot product for singles amplitudes (both spins)
    dotval += global_dpd_->file2_dot(RIA, CME);
    dotval += global_dpd_->file2_dot(Ria, Cme);

    // Compute dot product for doubles amplitudes (all spin cases)
    dotval += global_dpd_->buf4_dot(RIJAB, CMNEF);
    dotval += global_dpd_->buf4_dot(Rijab, Cmnef);
    dotval += global_dpd_->buf4_dot(RIjAb, CMnEf);

    return dotval;
}

void orthogonalize_C(dpdfile2 *RIA, dpdfile2 *Ria, dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb,
                     dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf,
                     double *overlap) {
    // Compute overlap <R|C>
    *overlap = dot_product_C(RIA, Ria, RIJAB, Rijab, RIjAb, CME, Cme, CMNEF, Cmnef, CMnEf);

    // Perform R = R - overlap * C (orthogonalization step)
    // This implements the Modified Gram-Schmidt algorithm: remove component of R along C
    global_dpd_->file2_axpy(CME, RIA, -1.0 * (*overlap), 0);
    global_dpd_->file2_axpy(Cme, Ria, -1.0 * (*overlap), 0);
    global_dpd_->buf4_axpy(CMNEF, RIJAB, -1.0 * (*overlap));
    global_dpd_->buf4_axpy(Cmnef, Rijab, -1.0 * (*overlap));
    global_dpd_->buf4_axpy(CMnEf, RIjAb, -1.0 * (*overlap));
}

double dot_product_C_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, dpdfile2 *CME, dpdbuf4 *CMnEf, int irrep) {
    double dotval = 0.0;
    dpdbuf4 R2a, R2b;

    // RHF singles: factor of 2 for spin adaptation
    dotval = 2.0 * global_dpd_->file2_dot(RIA, CME);

    // RHF doubles: spin-adapt the residual
    // Compute <2*RIjAb - RIjbA | CMnEf>
    global_dpd_->buf4_copy(RIjAb, PSIF_EOM_TMP, "RIjAb_tmp");
    global_dpd_->buf4_sort(RIjAb, PSIF_EOM_TMP, pqsr, 0, 5, "RIjbA_tmp");

    global_dpd_->buf4_init(&R2a, PSIF_EOM_TMP, irrep, 0, 5, 0, 5, 0, "RIjAb_tmp");
    global_dpd_->buf4_init(&R2b, PSIF_EOM_TMP, irrep, 0, 5, 0, 5, 0, "RIjbA_tmp");

    // Form 2*RIjAb - RIjbA
    global_dpd_->buf4_scm(&R2a, 2.0);
    global_dpd_->buf4_axpy(&R2b, &R2a, -1.0);

    // Dot product with basis vector
    dotval += global_dpd_->buf4_dot(&R2a, CMnEf);

    global_dpd_->buf4_close(&R2a);
    global_dpd_->buf4_close(&R2b);

    return dotval;
}

void orthogonalize_C_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, dpdfile2 *CME, dpdbuf4 *CMnEf, double *overlap, int irrep) {
    // Compute overlap with proper RHF spin adaptation
    *overlap = dot_product_C_RHF(RIA, RIjAb, CME, CMnEf, irrep);

    // Orthogonalize: R = R - overlap * C
    global_dpd_->file2_axpy(CME, RIA, -1.0 * (*overlap), 0);
    global_dpd_->buf4_axpy(CMnEf, RIjAb, -1.0 * (*overlap));
}

}  // namespace cceom
}  // namespace psi
