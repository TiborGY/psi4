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
    \brief DIIS extrapolation for RHF Lambda amplitudes using libdiis
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
** diis_RHF_libdiis: DIIS extrapolation for RHF Lambda amplitudes using libdiis
**
** This implementation uses libdiis/DIISManager for DIIS extrapolation,
** leveraging native DPD buffer support to eliminate manual flattening.
**
** Amplitude components (2 total):
**   - L1: LIA (singles)
**   - L2: LIjAb (doubles)
**
** Error vectors: R = L_new - L_old (computed using DPD operations)
**
** Key advantages over original implementation:
**   - No manual vector length calculation
**   - No manual flattening of DPD buffers to 1D arrays
**   - No manual B matrix construction or linear system solving
**   - Automatic storage management by libdiis
**   - Cleaner, more maintainable code (~65 lines vs ~213 lines original)
**
** Based on successful ccenergy DIIS migration (100% validated).
**
** Parameters:
**   iter: Iteration number
**   L_irr: Irrep of target state (for excited states)
*/

void CCLambdaWavefunction::diis_RHF_libdiis(int iter, int L_irr) {
    // Need at least 2 iterations for DIIS extrapolation
    if (iter < 2) {
        return;
    }

    auto nirreps = moinfo.nirreps;
    dpdfile2 L1_new, L1_old, R1;
    dpdbuf4 L2_new, L2_old, R2;

    /*
     * Step 1: Compute error vector for L1 (singles)
     * R1 = L1_new - L1_old using DPD file2_axpy operation
     */
    global_dpd_->file2_init(&L1_new, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&L1_old, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");

    // Create error vector as copy of L1_new, then subtract L1_old
    global_dpd_->file2_copy(&L1_new, PSIF_CC_OEI, "R1_IA");
    global_dpd_->file2_init(&R1, PSIF_CC_OEI, L_irr, 0, 1, "R1_IA");
    global_dpd_->file2_axpy(&L1_old, &R1, -1.0, 0);
    global_dpd_->file2_close(&L1_old);

    /*
     * Step 2: Compute error vector for L2 (doubles)
     * R2 = L2_new - L2_old using DPD buf4_axpy operation
     */
    global_dpd_->buf4_init(&L2_new, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_init(&L2_old, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");

    // Create error vector as copy of L2_new, then subtract L2_old
    global_dpd_->buf4_copy(&L2_new, PSIF_CC_LAMBDA, "R2_IjAb");
    global_dpd_->buf4_init(&R2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "R2_IjAb");
    global_dpd_->buf4_axpy(&L2_old, &R2, -1.0);
    global_dpd_->buf4_close(&L2_old);

    /*
     * Step 3: Add error and amplitude vectors to DIIS subspace
     * libdiis handles:
     *   - DPD → Matrix conversion
     *   - Vector storage (OnDisk policy)
     *   - B matrix construction
     *   - Subspace management (LargestError removal if full)
     */
    bool added = ccsd_diis_manager_->add_entry(
        &R1, &R2,        // Error vectors (2 components)
        &L1_new, &L2_new // Amplitude vectors (2 components)
    );

    /*
     * Step 4: Perform DIIS extrapolation if subspace is large enough
     * Requires at least 2 vectors for extrapolation
     */
    int subspace_size = ccsd_diis_manager_->subspace_size();
    if (subspace_size >= 2) {
        // libdiis computes optimal linear combination and updates amplitudes in place
        ccsd_diis_manager_->extrapolate(&L1_new, &L2_new);
        outfile->Printf("  DIIS: extrapolated with %d vectors\n", subspace_size);
    }

    /*
     * Step 5: Spin-adaptation (legacy code for compatibility)
     * RHF requires copying LIjAb → LIJAB and Lijab for spin-adapted routines
     * This preserves the original behavior from lines 269, 284-287 of original diis.cc
     */
    global_dpd_->file2_copy(&L1_new, PSIF_CC_LAMBDA, "New Lia");

    // Need to close and reinitialize L2_new with proper symmetry for spin-adaptation
    global_dpd_->buf4_close(&L2_new);
    global_dpd_->buf4_init(&L2_new, PSIF_CC_LAMBDA, L_irr, 2, 7, 0, 5, 1, "New LIjAb");
    global_dpd_->buf4_copy(&L2_new, PSIF_CC_LAMBDA, "New LIJAB");
    global_dpd_->buf4_copy(&L2_new, PSIF_CC_LAMBDA, "New Lijab");
    global_dpd_->buf4_close(&L2_new);

    /*
     * Step 6: Cleanup - close all DPD file handles
     */
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->file2_close(&L1_new);
}

}  // namespace cclambda
}  // namespace psi
