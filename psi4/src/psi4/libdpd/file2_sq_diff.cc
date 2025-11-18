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
    \brief Computes the sum of squared differences between two dpdfile2 objects.
*/
#include <cstdio>
#include "dpd.h"

namespace psi {

/*!
** file2_sq_diff(): Computes the sum of squared differences between
** two dpdfile2 objects. This is useful for convergence checking in iterative
** CC methods.
**
** \param File1 = pointer to the first dpdfile2 (typically "new" amplitudes)
** \param File2 = pointer to the second dpdfile2 (typically "old" amplitudes)
**
** Returns: sum of (File1[i] - File2[i])^2 over all elements
**
** Note: The square root is NOT taken, allowing the caller to accumulate
** contributions from multiple files before taking sqrt once at the end.
*/
double DPD::file2_sq_diff(dpdfile2 *File1, dpdfile2 *File2) {
    int h, nirreps, my_irrep;
    int row, col;
    double sum_sq_diff = 0.0, diff;

    nirreps = File1->params->nirreps;
    my_irrep = File1->my_irrep;

    // Initialize and read both files
    file2_mat_init(File1);
    file2_mat_rd(File1);
    file2_mat_init(File2);
    file2_mat_rd(File2);

    // Compute sum of squared differences
    for (h = 0; h < nirreps; h++)
        for (row = 0; row < File1->params->rowtot[h]; row++)
            for (col = 0; col < File1->params->coltot[h ^ my_irrep]; col++) {
                diff = File1->matrix[h][row][col] - File2->matrix[h][row][col];
                sum_sq_diff += diff * diff;
            }

    // Clean up
    file2_mat_close(File1);
    file2_mat_close(File2);

    return sum_sq_diff;
}

}  // namespace psi
