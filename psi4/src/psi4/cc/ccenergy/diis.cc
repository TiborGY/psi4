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
#include "psi4/libmints/matrix.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "Params.h"
#include "psi4/cc/ccwave.h"

namespace psi {
namespace ccenergy {

/*
** DIIS: Direct inversion in the iterative subspace routine to
** accelerate convergence of the CCSD amplitude equations.
**
** Implementation uses libdiis/DIISManager for all reference types,
** leveraging centralized DIIS infrastructure and native DPD buffer support.
**
** Original custom implementations replaced 2025-11-18:
** - Eliminates ~980 lines of duplicate DIIS code
** - Uses native DPD operations (file2_axpy, buf4_axpy)
** - Automatic B matrix construction and linear system solving
** - Consistent DIIS behavior with occ/dfocc modules
*/

void CCEnergyWavefunction::diis(int iter) {
    // Use libdiis implementation for all reference types
    if (params_.ref == 0)
        diis_RHF_libdiis(iter);
    else if (params_.ref == 1)
        diis_ROHF_libdiis(iter);
    else if (params_.ref == 2)
        diis_UHF_libdiis(iter);
}

void CCEnergyWavefunction::diis_invert_B(double** B, double* C, int dimension, double tolerance) {
    auto B2 = std::make_shared<Matrix>("B2", dimension, dimension);
    double** Bp = B2->pointer();
    ::memcpy((void*)Bp[0], B[0], sizeof(double) * dimension * dimension);

    auto* Sp = new double[dimension];
    auto* Tp = new double[dimension];

    bool is_zero = false;
    for (int i = 0; i < dimension - 1; i++) {
        if (Bp[i][i] <= 0.0) is_zero = true;
    }

    if (is_zero) {
        for (int i = 0; i < dimension; i++) {
            Sp[i] = 1.0;
        }
    } else {
        for (int i = 0; i < dimension - 1; i++) {
            Sp[i] = pow(Bp[i][i], -1.0 / 2.0);
        }
        Sp[dimension - 1] = 1.0;
    }

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            Bp[i][j] *= Sp[i] * Sp[j];
        }
    }

    B2->power(-1.0, tolerance);

    C_DGEMV('N', dimension, dimension, 1.0, Bp[0], dimension, C, 1, 0.0, Tp, 1);

    for (int i = 0; i < dimension; i++) {
        C[i] = Sp[i] * Tp[i];
    }

    delete[] Sp;
    delete[] Tp;
}

}  // namespace ccenergy
}  // namespace psi
