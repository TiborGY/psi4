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
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "psi4/cc/ccwave.h"

namespace psi {
namespace ccenergy {

int CCEnergyWavefunction::converged(double ediff) {
    double rms = 0.0;
    dpdfile2 T1, T1old;
    dpdbuf4 T2, T2old;

    if (params_.ref == 0) { /** RHF **/
        // T1 amplitudes
        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
        global_dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 0, 1, "tIA");
        rms += global_dpd_->file2_sq_diff(&T1, &T1old);
        global_dpd_->file2_close(&T1);
        global_dpd_->file2_close(&T1old);

        // T2 amplitudes
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
        global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        rms += global_dpd_->buf4_sq_diff(&T2, &T2old);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_close(&T2old);

    } else if (params_.ref == 1) { /** ROHF **/
        // T1 (IA) amplitudes
        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
        global_dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 0, 1, "tIA");
        rms += global_dpd_->file2_sq_diff(&T1, &T1old);
        global_dpd_->file2_close(&T1);
        global_dpd_->file2_close(&T1old);

        // T1 (ia) amplitudes
        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "New tia");
        global_dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 0, 1, "tia");
        rms += global_dpd_->file2_sq_diff(&T1, &T1old);
        global_dpd_->file2_close(&T1);
        global_dpd_->file2_close(&T1old);

        // T2 (IJAB) amplitudes
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
        global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
        rms += global_dpd_->buf4_sq_diff(&T2, &T2old);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_close(&T2old);

        // T2 (ijab) amplitudes
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
        global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
        rms += global_dpd_->buf4_sq_diff(&T2, &T2old);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_close(&T2old);

        // T2 (IjAb) amplitudes
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
        global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        rms += global_dpd_->buf4_sq_diff(&T2, &T2old);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_close(&T2old);

    } else if (params_.ref == 2) { /** UHF **/
        // T1 (IA) amplitudes
        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
        global_dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 0, 1, "tIA");
        rms += global_dpd_->file2_sq_diff(&T1, &T1old);
        global_dpd_->file2_close(&T1);
        global_dpd_->file2_close(&T1old);

        // T1 (ia) amplitudes
        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "New tia");
        global_dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 2, 3, "tia");
        rms += global_dpd_->file2_sq_diff(&T1, &T1old);
        global_dpd_->file2_close(&T1);
        global_dpd_->file2_close(&T1old);

        // T2 (IJAB) amplitudes
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
        global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
        rms += global_dpd_->buf4_sq_diff(&T2, &T2old);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_close(&T2old);

        // T2 (ijab) amplitudes
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
        global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
        rms += global_dpd_->buf4_sq_diff(&T2, &T2old);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_close(&T2old);

        // T2 (IjAb) amplitudes
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
        global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
        rms += global_dpd_->buf4_sq_diff(&T2, &T2old);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_close(&T2old);
    }

    rms = sqrt(rms);
    moinfo_.conv = rms;

    if ((rms < params_.convergence) && (std::fabs(ediff) < params_.e_convergence))
        return 1;
    else
        return 0;
}
}  // namespace ccenergy
}  // namespace psi
