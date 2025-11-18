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
    \ingroup DPD
    \brief Computes the sum of squared differences between two dpdbuf4 objects.
*/
#include <cstdio>
#include "dpd.h"

namespace psi {

/*!
** buf4_rms_diff(): Computes the sum of squared differences (RMS^2) between
** two dpdbuf4 objects. This is useful for convergence checking in iterative
** CC methods.
**
** \param Buf1 = pointer to the first dpdbuf4 (typically "new" amplitudes)
** \param Buf2 = pointer to the second dpdbuf4 (typically "old" amplitudes)
**
** Returns: sum of (Buf1[i] - Buf2[i])^2 over all elements
**
** Note: The square root is NOT taken, allowing the caller to accumulate
** contributions from multiple buffers before taking sqrt once.
*/
double DPD::buf4_rms_diff(dpdbuf4 *Buf1, dpdbuf4 *Buf2) {
    int h, nirreps, my_irrep;
    int row, col;
    double rms = 0.0, diff;

    nirreps = Buf1->params->nirreps;
    my_irrep = Buf1->file.my_irrep;

    // Process each irrep block
    for (h = 0; h < nirreps; h++) {
        buf4_mat_irrep_init(Buf1, h);
        buf4_mat_irrep_rd(Buf1, h);
        buf4_mat_irrep_init(Buf2, h);
        buf4_mat_irrep_rd(Buf2, h);

        // Compute sum of squared differences for this irrep
        for (row = 0; row < Buf1->params->rowtot[h]; row++)
            for (col = 0; col < Buf1->params->coltot[h ^ my_irrep]; col++) {
                diff = Buf1->matrix[h][row][col] - Buf2->matrix[h][row][col];
                rms += diff * diff;
            }

        buf4_mat_irrep_close(Buf1, h);
        buf4_mat_irrep_close(Buf2, h);
    }

    return rms;
}

}  // namespace psi
