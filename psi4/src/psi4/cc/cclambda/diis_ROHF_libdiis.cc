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
    \brief DIIS extrapolation for ROHF Lambda amplitudes using libdiis
*/

#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/psifiles.h"
#include "cclambda.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cclambda {

/*
** diis_ROHF_libdiis: DIIS extrapolation for ROHF Lambda amplitudes using libdiis
**
** This implementation uses libdiis/DIISManager for DIIS extrapolation,
** leveraging native DPD buffer support to eliminate manual flattening.
**
** Amplitude components (5 total):
**   - L1: LIA, Lia (2 components)
**   - L2: LIJAB, Lijab, LIjAb (3 components)
**
** Error vectors: R = L_new - L_old (computed using DPD operations)
**
** Key advantages over original implementation:
**   - No manual vector length calculation
**   - No manual flattening of DPD buffers to 1D arrays
**   - No manual B matrix construction or linear system solving
**   - Automatic storage management by libdiis
**   - Cleaner, more maintainable code (~100 lines vs ~275 lines original)
**
** Based on successful ccenergy ROHF DIIS migration (100% validated).
**
** Parameters:
**   iter: Iteration number
**   L_irr: Irrep of target state (for excited states)
*/

void CCLambdaWavefunction::diis_ROHF_libdiis(int iter, int L_irr) {
    // Need at least 2 iterations for DIIS extrapolation
    if (iter < 2) {
        return;
    }

    auto nirreps = moinfo.nirreps;
    dpdfile2 L1a_new, L1a_old, R1a;
    dpdfile2 L1b_new, L1b_old, R1b;
    dpdbuf4 L2aa_new, L2aa_old, R2aa;
    dpdbuf4 L2bb_new, L2bb_old, R2bb;
    dpdbuf4 L2ab_new, L2ab_old, R2ab;

    /*
     * Step 1: Compute error vectors for all 5 amplitude components
     * Using DPD file2_axpy and buf4_axpy operations
     */

    // L1a (LIA) - Alpha singles
    global_dpd_->file2_init(&L1a_new, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&L1a_old, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_copy(&L1a_new, PSIF_CC_OEI, "R1a_IA");
    global_dpd_->file2_init(&R1a, PSIF_CC_OEI, L_irr, 0, 1, "R1a_IA");
    global_dpd_->file2_axpy(&L1a_old, &R1a, -1.0, 0);
    global_dpd_->file2_close(&L1a_old);

    // L1b (Lia) - Beta singles
    global_dpd_->file2_init(&L1b_new, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia");
    global_dpd_->file2_init(&L1b_old, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
    global_dpd_->file2_copy(&L1b_new, PSIF_CC_OEI, "R1b_ia");
    global_dpd_->file2_init(&R1b, PSIF_CC_OEI, L_irr, 0, 1, "R1b_ia");
    global_dpd_->file2_axpy(&L1b_old, &R1b, -1.0, 0);
    global_dpd_->file2_close(&L1b_old);

    // L2aa (LIJAB) - Alpha-alpha doubles
    global_dpd_->buf4_init(&L2aa_new, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_init(&L2aa_old, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_copy(&L2aa_new, PSIF_CC_LAMBDA, "R2aa_IJAB");
    global_dpd_->buf4_init(&R2aa, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "R2aa_IJAB");
    global_dpd_->buf4_axpy(&L2aa_old, &R2aa, -1.0);
    global_dpd_->buf4_close(&L2aa_old);

    // L2bb (Lijab) - Beta-beta doubles
    global_dpd_->buf4_init(&L2bb_new, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    global_dpd_->buf4_init(&L2bb_old, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_copy(&L2bb_new, PSIF_CC_LAMBDA, "R2bb_ijab");
    global_dpd_->buf4_init(&R2bb, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "R2bb_ijab");
    global_dpd_->buf4_axpy(&L2bb_old, &R2bb, -1.0);
    global_dpd_->buf4_close(&L2bb_old);

    // L2ab (LIjAb) - Alpha-beta doubles
    global_dpd_->buf4_init(&L2ab_new, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_init(&L2ab_old, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_copy(&L2ab_new, PSIF_CC_LAMBDA, "R2ab_IjAb");
    global_dpd_->buf4_init(&R2ab, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "R2ab_IjAb");
    global_dpd_->buf4_axpy(&L2ab_old, &R2ab, -1.0);
    global_dpd_->buf4_close(&L2ab_old);

    /*
     * Step 2: Add error and amplitude vectors to DIIS subspace
     * libdiis variadic template handles 5 amplitude components
     */
    bool added = ccsd_diis_manager_->add_entry(
        &R1a, &R1b, &R2aa, &R2bb, &R2ab,                        // Error vectors (5 components)
        &L1a_new, &L1b_new, &L2aa_new, &L2bb_new, &L2ab_new     // Amplitude vectors (5 components)
    );

    /*
     * Step 3: Perform DIIS extrapolation if subspace is large enough
     * Requires at least 2 vectors for extrapolation
     */
    int subspace_size = ccsd_diis_manager_->subspace_size();
    if (subspace_size >= 2) {
        // libdiis computes optimal linear combination and updates amplitudes in place
        ccsd_diis_manager_->extrapolate(&L1a_new, &L1b_new, &L2aa_new, &L2bb_new, &L2ab_new);
        outfile->Printf("  DIIS: extrapolated with %d vectors\n", subspace_size);
    }

    /*
     * Step 4: Cleanup - close all DPD file handles
     */
    global_dpd_->file2_close(&R1a);
    global_dpd_->file2_close(&R1b);
    global_dpd_->buf4_close(&R2aa);
    global_dpd_->buf4_close(&R2bb);
    global_dpd_->buf4_close(&R2ab);
    global_dpd_->file2_close(&L1a_new);
    global_dpd_->file2_close(&L1b_new);
    global_dpd_->buf4_close(&L2aa_new);
    global_dpd_->buf4_close(&L2bb_new);
    global_dpd_->buf4_close(&L2ab_new);
}

}  // namespace cclambda
}  // namespace psi
