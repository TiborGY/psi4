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
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "cclambda.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cclambda {

int CCLambdaWavefunction::converged(int L_irr) {
    double rms = 0.0;
    dpdfile2 L1, L1old;
    dpdbuf4 L2, L2old;

    // L1 (IA) - all reference types
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&L1old, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    rms += global_dpd_->file2_sq_diff(&L1, &L1old);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&L1old);

    if (params.ref == 0) {
        rms *= 2.0;  // RHF: account for spin symmetry
    }

    // L1 (ia) - ROHF/UHF only
    if (params.ref == 1) {
        global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia");
        global_dpd_->file2_init(&L1old, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
        rms += global_dpd_->file2_sq_diff(&L1, &L1old);
        global_dpd_->file2_close(&L1);
        global_dpd_->file2_close(&L1old);
    } else if (params.ref == 2) {
        global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 2, 3, "New Lia");
        global_dpd_->file2_init(&L1old, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");
        rms += global_dpd_->file2_sq_diff(&L1, &L1old);
        global_dpd_->file2_close(&L1);
        global_dpd_->file2_close(&L1old);
    }

    // L2 (IJAB) - ROHF/UHF only
    if (params.ref == 1 || params.ref == 2) {
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
        global_dpd_->buf4_init(&L2old, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
        rms += global_dpd_->buf4_sq_diff(&L2, &L2old);
        global_dpd_->buf4_close(&L2);
        global_dpd_->buf4_close(&L2old);
    }

    // L2 (ijab) - ROHF/UHF only
    if (params.ref == 1) {
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
        global_dpd_->buf4_init(&L2old, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
        rms += global_dpd_->buf4_sq_diff(&L2, &L2old);
        global_dpd_->buf4_close(&L2);
        global_dpd_->buf4_close(&L2old);
    } else if (params.ref == 2) {
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
        global_dpd_->buf4_init(&L2old, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
        rms += global_dpd_->buf4_sq_diff(&L2, &L2old);
        global_dpd_->buf4_close(&L2);
        global_dpd_->buf4_close(&L2old);
    }

    // L2 (IjAb) - all reference types
    if (params.ref == 0 || params.ref == 1) {
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
        global_dpd_->buf4_init(&L2old, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    } else if (params.ref == 2) {
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
        global_dpd_->buf4_init(&L2old, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    }
    rms += global_dpd_->buf4_sq_diff(&L2, &L2old);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&L2old);

    rms = sqrt(rms);
    moinfo.conv = rms;

    if (rms < params.convergence)
        return 1;
    else
        return 0;
}

}  // namespace cclambda
}  // namespace psi
