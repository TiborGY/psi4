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
    \brief RHF DIIS implementation using libdiis/DIISManager (PROOF OF CONCEPT)

    This is a proof-of-concept implementation to validate that libdiis can
    successfully replace the custom DIIS implementation for RHF CCSD amplitudes.

    Key advantages over custom implementation:
    - Eliminates ~200 lines of manual buffer management code
    - libdiis handles DPD buffer conversion automatically
    - Automatic B matrix construction and linear system solving
    - Consistent DIIS behavior with other Psi4 modules (occ, dfocc, dct)

    Technical approach:
    - Compute error vectors (R1, R2) using DPD operations
    - Pass DPD buffers directly to DIISManager
    - libdiis converts DPD → Matrix internally for storage/extrapolation
    - libdiis extrapolates and updates T1_new, T2_new in-place

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
#include "psi4/libdiis/diismanager.h"  // NEW: libdiis interface
#include "MOInfo.h"
#include "Params.h"
#include "psi4/cc/ccwave.h"

namespace psi {
namespace ccenergy {

/*
** diis_RHF_libdiis: DIIS extrapolation for RHF CCSD amplitudes using libdiis
**
** This function replaces ~250 lines of custom DIIS code with ~60 lines that
** leverage the centralized libdiis/DIISManager infrastructure.
**
** Parameters:
**   iter - Current CC iteration number
**
** Algorithm:
**   1. Compute error vectors: R1 = T1_new - T1_old, R2 = T2_new - T2_old
**   2. Add error and amplitude vectors to DIIS subspace
**   3. Extrapolate new amplitudes (if enough vectors available)
**
** Note: DIISManager was initialized in CCEnergyWavefunction constructor
*/
void CCEnergyWavefunction::diis_RHF_libdiis(int iter) {
    // Early exit if not enough iterations for DIIS
    if (iter < 2) {
        return;  // Need at least 2 iterations
    }

    auto nirreps = moinfo_.nirreps;
    dpdfile2 T1_new, T1_old, R1;
    dpdbuf4 T2_new, T2_old, R2;

    // ==========================================
    // Step 1: Compute error vector for T1
    // ==========================================
    //
    // R1 = T1_new - T1_old
    //
    // Using DPD operations instead of manual flattening:
    // - Initialize R1 as copy of T1_new
    // - Subtract T1_old using axpy operation

    global_dpd_->file2_init(&T1_new, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_init(&T1_old, PSIF_CC_TAMPS, 0, 0, 1, "tIA");

    // Copy T1_new to R1
    global_dpd_->file2_copy(&T1_new, PSIF_CC_OEI, "R1_IA");
    global_dpd_->file2_init(&R1, PSIF_CC_OEI, 0, 0, 1, "R1_IA");

    // R1 = R1 - T1_old (axpy: Y = Y + alpha*X, with alpha=-1.0)
    global_dpd_->file2_axpy(&T1_old, &R1, -1.0, 0);

    global_dpd_->file2_close(&T1_old);
    // Keep R1 and T1_new open for add_entry() call

    // ==========================================
    // Step 2: Compute error vector for T2
    // ==========================================
    //
    // R2 = T2_new - T2_old
    //
    // Same approach as T1, using DPD operations

    global_dpd_->buf4_init(&T2_new, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_init(&T2_old, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");

    // Copy T2_new to R2
    global_dpd_->buf4_copy(&T2_new, PSIF_CC_TAMPS, "R2_IjAb");
    global_dpd_->buf4_init(&R2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "R2_IjAb");

    // R2 = R2 - T2_old
    global_dpd_->buf4_axpy(&T2_old, &R2, -1.0);

    global_dpd_->buf4_close(&T2_old);
    // Keep R2 and T2_new open for add_entry() call

    // ==========================================
    // Step 3: Add vectors to DIIS subspace
    // ==========================================
    //
    // libdiis expects: error vectors first, then amplitude vectors
    // It will automatically:
    // - Convert DPD buffers to Matrix format
    // - Store to PSIF_LIBDIIS file
    // - Update B matrix
    // - Manage subspace size (remove oldest or largest error)

    bool added = ccsd_diis_manager_->add_entry(&R1, &R2, &T1_new, &T2_new);

    if (!added) {
        outfile->Printf("  DIIS: Warning - add_entry failed at iteration %d\n", iter);
    }

    // ==========================================
    // Step 4: Extrapolate if enough vectors
    // ==========================================
    //
    // Once we have at least 2 vectors (mindiis), perform extrapolation
    // libdiis will:
    // - Build B matrix from stored error vectors
    // - Solve DIIS linear system with conditioning
    // - Compute extrapolated amplitudes
    // - Update T1_new and T2_new in-place

    int subspace_size = ccsd_diis_manager_->subspace_size();

    if (subspace_size >= 2) {
        // Perform extrapolation
        ccsd_diis_manager_->extrapolate(&T1_new, &T2_new);

        outfile->Printf("  DIIS: extrapolated with %d vectors\n", subspace_size);
    } else {
        outfile->Printf("  DIIS: stored vector (subspace = %d, need >= 2)\n", subspace_size);
    }

    // ==========================================
    // Step 5: Cleanup
    // ==========================================
    //
    // Close all DPD files
    // Note: T1_new and T2_new contain extrapolated values if DIIS was performed

    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->file2_close(&T1_new);
    global_dpd_->buf4_close(&T2_new);
}

}  // namespace ccenergy
}  // namespace psi
