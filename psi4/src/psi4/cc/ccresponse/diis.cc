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
    \ingroup ccresponse
    \brief DIIS extrapolation for response amplitudes using libdiis
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <map>
#include <memory>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {

// Static map to store DIISManager instances per perturbation/frequency combination
// Each perturbation (e.g., "Mu", "P", "L") at each frequency needs its own DIIS history
static std::map<std::string, std::shared_ptr<DIISManager>> diis_managers_;

/*
** DIIS: Direct inversion in the iterative subspace routine to
** accelerate convergence of the response amplitude equations.
**
** This implementation uses libdiis/DIISManager for DIIS extrapolation,
** leveraging native DPD buffer support to eliminate manual flattening.
**
** Response amplitude components (2 total for RHF):
**   - X1: X_pert_IA (singles, e.g., "X_Mu_IA (0.000)")
**   - X2: X_pert_IjAb (doubles, e.g., "X_Mu_IjAb (0.000)")
**
** Error vectors: R = X_new - X_old (computed using DPD operations)
**
** Key features:
**   - Per-perturbation DIIS: Each perturbation/frequency has its own DIISManager
**   - Frequency-dependent: Labels include omega (static: 0.000, dynamic: e.g., 0.072)
**   - No manual flattening: Uses DPD file2_axpy and buf4_axpy operations
**   - Automatic B matrix construction and linear system solving via libdiis
**
** Original custom implementation replaced 2025-11-18:
** - Eliminates ~193 lines of manual DIIS code
** - Uses native DPD operations (file2_axpy, buf4_axpy)
** - Automatic storage management by libdiis
** - Cleaner, more maintainable code (~80 lines vs ~193 lines original)
**
** Based on successful ccenergy and cclambda DIIS migrations (100% validated).
**
** Parameters:
**   iter: Iteration number
**   pert: Perturbation name (e.g., "Mu", "P", "L")
**   irrep: Irrep of perturbation
**   omega: Frequency (0.0 for static properties, non-zero for dynamic)
*/

void diis(int iter, const char *pert, int irrep, double omega) {
    // Need at least 2 iterations for DIIS extrapolation
    if (iter < 2) {
        return;
    }

    auto nirreps = moinfo.nirreps;

    // Create unique key for this perturbation/frequency combination
    // Format: "Mu_0.000000" or "P_0.072000"
    char omega_str[32];
    sprintf(omega_str, "%.6f", omega);
    std::string key = std::string(pert) + "_" + omega_str;

    // Initialize DIISManager on first use for this perturbation/frequency
    if (diis_managers_.find(key) == diis_managers_.end()) {
        std::string label = std::string("Response DIIS ") + pert;
        if (omega != 0.0) {
            label += " ω=" + std::string(omega_str);
        }

        diis_managers_[key] = std::make_shared<DIISManager>(
            8,                                          // max 8 vectors
            label,                                      // label with perturbation/frequency
            DIISManager::RemovalPolicy::LargestError,  // removal policy
            DIISManager::StoragePolicy::OnDisk         // storage policy (same as original)
        );
    }

    auto& manager = diis_managers_[key];

    /*
     * Compute error vectors for response amplitudes
     * Response equations use X vectors instead of T or Lambda
     */
    char lbl[32];
    dpdfile2 X1_new, X1_old, R1;
    dpdbuf4 X2_new, X2_old, R2;

    /*
     * Step 1: Compute error vector for X1 (singles)
     * R1 = X1_new - X1_old using DPD file2_axpy operation
     */
    sprintf(lbl, "New X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1_new, PSIF_CC_OEI, irrep, 0, 1, lbl);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1_old, PSIF_CC_OEI, irrep, 0, 1, lbl);

    // Create error vector as copy of X1_new, then subtract X1_old
    global_dpd_->file2_copy(&X1_new, PSIF_CC_OEI, "R1_IA_tmp");
    global_dpd_->file2_init(&R1, PSIF_CC_OEI, irrep, 0, 1, "R1_IA_tmp");
    global_dpd_->file2_axpy(&X1_old, &R1, -1.0, 0);
    global_dpd_->file2_close(&X1_old);

    /*
     * Step 2: Compute error vector for X2 (doubles)
     * R2 = X2_new - X2_old using DPD buf4_axpy operation
     */
    sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2_new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2_old, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

    // Create error vector as copy of X2_new, then subtract X2_old
    global_dpd_->buf4_copy(&X2_new, PSIF_CC_LR, "R2_IjAb_tmp");
    global_dpd_->buf4_init(&R2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "R2_IjAb_tmp");
    global_dpd_->buf4_axpy(&X2_old, &R2, -1.0);
    global_dpd_->buf4_close(&X2_old);

    /*
     * Step 3: Add error and amplitude vectors to DIIS subspace
     * libdiis handles:
     *   - DPD → Matrix conversion
     *   - Vector storage (OnDisk policy)
     *   - B matrix construction
     *   - Subspace management (LargestError removal if full)
     */
    bool added = manager->add_entry(
        &R1, &R2,        // Error vectors (2 components)
        &X1_new, &X2_new // Amplitude vectors (2 components)
    );

    /*
     * Step 4: Perform DIIS extrapolation if subspace is large enough
     * Requires at least 2 vectors for extrapolation
     */
    int subspace_size = manager->subspace_size();
    if (subspace_size >= 2) {
        // libdiis computes optimal linear combination and updates amplitudes in place
        manager->extrapolate(&X1_new, &X2_new);
        outfile->Printf("  DIIS: extrapolated with %d vectors\n", subspace_size);
    }

    /*
     * Step 5: Cleanup - close all DPD file handles
     */
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&R2);
    global_dpd_->file2_close(&X1_new);
    global_dpd_->buf4_close(&X2_new);
}

}  // namespace ccresponse
}  // namespace psi
