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
    \ingroup CCENERGY
    \brief UHF DIIS implementation using libdiis/DIISManager (PROOF OF CONCEPT)

    This extends the POC to UHF references, handling separate alpha and beta
    spin components.

    UHF amplitudes:
    - T1: alpha and beta singles (T1a, T1b)
    - T2: alpha-alpha, beta-beta, and alpha-beta doubles (T2aa, T2bb, T2ab)

    Code reduction: ~365 lines → ~110 lines (70% reduction)

    Author: Claude (Anthropic AI)
    Date: 2025-11-18
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/libdiis/diismanager.h"
#include "MOInfo.h"
#include "Params.h"
#include "psi4/cc/ccwave.h"

namespace psi {
namespace ccenergy {

/*
** diis_UHF_libdiis: DIIS extrapolation for UHF CCSD amplitudes using libdiis
**
** Replaces ~365 lines of custom DIIS code with ~110 lines using libdiis.
**
** UHF requires handling 5 amplitude components (T1a, T1b, T2aa, T2bb, T2ab)
** libdiis variadic templates handle this seamlessly.
**
** Algorithm:
**   1. Compute error vectors for all 5 amplitude components
**   2. Add all error and amplitude vectors to DIIS subspace
**   3. Extrapolate if enough vectors available
*/

void CCEnergyWavefunction::diis_UHF_libdiis(int iter) {
    // Need at least 2 iterations for DIIS
    if (iter < 2) {
        return;
    }

    auto nirreps = moinfo_.nirreps;

    // DPD file/buffer objects for amplitudes and error vectors
    dpdfile2 T1a_new, T1a_old, R1a;
    dpdfile2 T1b_new, T1b_old, R1b;
    dpdbuf4 T2aa_new, T2aa_old, R2aa;
    dpdbuf4 T2bb_new, T2bb_old, R2bb;
    dpdbuf4 T2ab_new, T2ab_old, R2ab;

    //
    // Step 1: Compute error vectors using DPD operations
    //

    // T1a amplitudes (alpha singles)
    // R1a = T1a_new - T1a_old
    global_dpd_->file2_init(&T1a_new, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_init(&T1a_old, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_copy(&T1a_new, PSIF_CC_OEI, "R1a_IA");
    global_dpd_->file2_init(&R1a, PSIF_CC_OEI, 0, 0, 1, "R1a_IA");
    global_dpd_->file2_axpy(&T1a_old, &R1a, -1.0, 0);
    global_dpd_->file2_close(&T1a_old);

    // T1b amplitudes (beta singles)
    // R1b = T1b_new - T1b_old
    global_dpd_->file2_init(&T1b_new, PSIF_CC_OEI, 0, 2, 3, "New tia");
    global_dpd_->file2_init(&T1b_old, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_copy(&T1b_new, PSIF_CC_OEI, "R1b_ia");
    global_dpd_->file2_init(&R1b, PSIF_CC_OEI, 0, 2, 3, "R1b_ia");
    global_dpd_->file2_axpy(&T1b_old, &R1b, -1.0, 0);
    global_dpd_->file2_close(&T1b_old);

    // T2aa amplitudes (alpha-alpha)
    // R2aa = T2aa_new - T2aa_old
    global_dpd_->buf4_init(&T2aa_new, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_init(&T2aa_old, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_copy(&T2aa_new, PSIF_CC_TAMPS, "R2aa_IJAB");
    global_dpd_->buf4_init(&R2aa, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "R2aa_IJAB");
    global_dpd_->buf4_axpy(&T2aa_old, &R2aa, -1.0);
    global_dpd_->buf4_close(&T2aa_old);

    // T2bb amplitudes (beta-beta)
    // R2bb = T2bb_new - T2bb_old
    global_dpd_->buf4_init(&T2bb_new, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
    global_dpd_->buf4_init(&T2bb_old, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
    global_dpd_->buf4_copy(&T2bb_new, PSIF_CC_TAMPS, "R2bb_ijab");
    global_dpd_->buf4_init(&R2bb, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "R2bb_ijab");
    global_dpd_->buf4_axpy(&T2bb_old, &R2bb, -1.0);
    global_dpd_->buf4_close(&T2bb_old);

    // T2ab amplitudes (alpha-beta)
    // R2ab = T2ab_new - T2ab_old
    global_dpd_->buf4_init(&T2ab_new, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    global_dpd_->buf4_init(&T2ab_old, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_copy(&T2ab_new, PSIF_CC_TAMPS, "R2ab_IjAb");
    global_dpd_->buf4_init(&R2ab, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "R2ab_IjAb");
    global_dpd_->buf4_axpy(&T2ab_old, &R2ab, -1.0);
    global_dpd_->buf4_close(&T2ab_old);

    //
    // Step 2: Add error and amplitude vectors to DIIS subspace
    //
    // libdiis variadic template handles all 5 components
    //
    bool added = ccsd_diis_manager_->add_entry(
        &R1a, &R1b, &R2aa, &R2bb, &R2ab,              // Error vectors
        &T1a_new, &T1b_new, &T2aa_new, &T2bb_new, &T2ab_new  // Amplitude vectors
    );

    //
    // Step 3: Extrapolate if we have enough vectors in the subspace
    //
    int subspace_size = ccsd_diis_manager_->subspace_size();
    if (subspace_size >= 2) {
        // Extrapolate: updates all amplitude vectors in-place
        ccsd_diis_manager_->extrapolate(&T1a_new, &T1b_new, &T2aa_new, &T2bb_new, &T2ab_new);

        outfile->Printf("  DIIS: extrapolated with %d vectors\n", subspace_size);
    }

    //
    // Step 4: Cleanup - close all DPD objects
    //
    global_dpd_->file2_close(&R1a);
    global_dpd_->file2_close(&R1b);
    global_dpd_->buf4_close(&R2aa);
    global_dpd_->buf4_close(&R2bb);
    global_dpd_->buf4_close(&R2ab);
    global_dpd_->file2_close(&T1a_new);
    global_dpd_->file2_close(&T1b_new);
    global_dpd_->buf4_close(&T2aa_new);
    global_dpd_->buf4_close(&T2bb_new);
    global_dpd_->buf4_close(&T2ab_new);
}

}  // namespace ccenergy
}  // namespace psi
