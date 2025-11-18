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
    \brief DIIS extrapolation dispatcher for Lambda amplitudes
*/

#include "cclambda.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cclambda {

/*
** DIIS: Direct inversion in the iterative subspace routine to
** accelerate convergence of the Lambda amplitude equations.
**
** Implementation uses libdiis/DIISManager for all reference types,
** leveraging centralized DIIS infrastructure and native DPD buffer support.
**
** Original custom implementations replaced 2025-11-18:
** - Eliminates ~768 lines of duplicate DIIS code
** - Uses native DPD operations (file2_axpy, buf4_axpy)
** - Automatic B matrix construction and linear system solving
** - Consistent DIIS behavior with ccenergy, occ, dfocc modules
**
** Migration based on proven ccenergy DIIS consolidation (100% successful).
*/

void CCLambdaWavefunction::diis(int iter, int L_irr) {
    // Use libdiis implementation for all reference types
    if (params.ref == 0)
        diis_RHF_libdiis(iter, L_irr);
    else if (params.ref == 1)
        diis_ROHF_libdiis(iter, L_irr);
    else if (params.ref == 2)
        diis_UHF_libdiis(iter, L_irr);
}

}  // namespace cclambda
}  // namespace psi
