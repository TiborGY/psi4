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
    \ingroup CC
    \brief Shared utilities for coupled cluster methods
*/

#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "ccwave.h"

namespace psi {
namespace ccenergy {

/**
 * halftrans(): Transform the last two indices of a dpdbuf4 between MO and SO bases
 *
 * This function is used in AO-basis algorithms for CCSD(T) and Lambda equations.
 * It enables workflows like:
 *   (ij|ab) MO → (ij|pq) AO → contract with (pq|rs) → (ij|rs) AO → (ij|cd) MO
 *
 * This is DIFFERENT from libtrans half-transformations:
 * - libtrans: (AO|AO) → (MO|AO) → (MO|MO) [full transform in two steps]
 * - halftrans: (MO|MO) ↔ (MO|AO) [partial transform for AO-basis algorithms]
 *
 * Do NOT attempt to replace this with libtrans. The use cases are distinct:
 * - libtrans transforms AO integrals → MO integrals (forward only)
 * - halftrans transforms MO quantities ↔ mixed MO/AO format (bidirectional)
 *
 * Used by:
 * - ccenergy/BT2_AO.cc: AO-basis algorithm for (abcd) contribution to T2
 * - cclambda/BL2_AO.cc: AO-basis algorithm for Lambda equations
 *
 * Activated when user sets AO_BASIS = DISK or DIRECT (not default).
 *
 * Algorithm:
 * For type = 0 (MO → SO):
 *   1. Read MO buffer: (ij|cd) where c,d are MO indices
 *   2. Transform d: (ij|c,MO) × C2^T → (ij|c,SO) [intermediate X]
 *   3. Transform c: C1 × (ij|c,SO) → (ij|SO,SO) = (ij|pq)
 *   4. Write result to SO buffer with scaling: alpha*source + beta*target
 *
 * For type = 1 (SO → MO):
 *   1. Read SO buffer: (ij|pq) where p,q are SO indices
 *   2. Transform q: (ij|p,SO) × C2 → (ij|p,MO) [intermediate X]
 *   3. Transform p: C1^T × (ij|p,MO) → (ij|MO,MO) = (ij|cd)
 *   4. Write result to MO buffer with scaling
 *
 * @param Buf1 MO-basis dpdbuf4 (already initialized)
 * @param dpdnum1 DPD instance number for Buf1
 * @param Buf2 SO-basis dpdbuf4 (already initialized)
 * @param dpdnum2 DPD instance number for Buf2
 * @param C1 Left transformation matrix (SO x MO, symmetry blocked)
 * @param C2 Right transformation matrix (SO x MO, symmetry blocked)
 * @param nirreps Number of irreps
 * @param mo_row MO index offsets for each (h, Gc) pair
 * @param so_row SO index offsets for each (h, Gc) pair
 * @param mospi_left Number of MOs per irrep for left index (c in examples)
 * @param mospi_right Number of MOs per irrep for right index (d in examples)
 * @param sospi Number of SOs per irrep
 * @param type 0 = MO → SO; 1 = SO → MO
 * @param alpha Scaling factor for source buffer
 * @param beta Scaling factor for target buffer (allows accumulation)
 */
void CCEnergyWavefunction::halftrans(dpdbuf4 *Buf1, int dpdnum1, dpdbuf4 *Buf2, int dpdnum2, double ***C1, double ***C2,
                                     int nirreps, int **mo_row, int **so_row, int *mospi_left, int *mospi_right,
                                     int *sospi, int type, double alpha, double beta) {
    int Gd, cd, pq;
    double **X;

    for (int h = 0; h < nirreps; h++) {
        dpd_set_default(dpdnum1);
        global_dpd_->buf4_mat_irrep_init(Buf1, h);

        dpd_set_default(dpdnum2);
        global_dpd_->buf4_mat_irrep_init(Buf2, h);

        if (type == 0) { /* alpha * Buf1 --> beta * Buf2 */
            if (alpha != 0.0) {
                dpd_set_default(dpdnum1);
                global_dpd_->buf4_mat_irrep_rd(Buf1, h);
            }
            if (beta != 0.0) {
                dpd_set_default(dpdnum2);
                global_dpd_->buf4_mat_irrep_rd(Buf2, h);
            }
        }
        if (type == 1) { /* alpha * Buf2 --> beta * Buf1 */
            if (alpha != 0.0) {
                dpd_set_default(dpdnum2);
                global_dpd_->buf4_mat_irrep_rd(Buf2, h);
            }
            if (beta != 0.0) {
                dpd_set_default(dpdnum1);
                global_dpd_->buf4_mat_irrep_rd(Buf1, h);
            }
        }

        for (int Gc = 0; Gc < nirreps; Gc++) {
            Gd = h ^ Gc;

            cd = mo_row[h][Gc];
            pq = so_row[h][Gc];

            if (mospi_left[Gc] && mospi_right[Gd] && sospi[Gc] && sospi[Gd]) {
                if (type == 0) {  // MO → SO transformation
                    X = block_matrix(mospi_left[Gc], sospi[Gd]);

                    for (int ij = 0; ij < Buf1->params->rowtot[h]; ij++) {
                        // Step 1: Transform right index (d) from MO to SO
                        // X(c,q) = (ij|c,d) × C2(q,d)^T
                        C_DGEMM('n', 't', mospi_left[Gc], sospi[Gd], mospi_right[Gd], 1.0, &(Buf1->matrix[h][ij][cd]),
                                mospi_right[Gd], &(C2[Gd][0][0]), mospi_right[Gd], 0.0, &(X[0][0]), sospi[Gd]);

                        // Step 2: Transform left index (c) from MO to SO
                        // (ij|p,q) = C1(p,c) × X(c,q)
                        C_DGEMM('n', 'n', sospi[Gc], sospi[Gd], mospi_left[Gc], alpha, &(C1[Gc][0][0]), mospi_left[Gc],
                                &(X[0][0]), sospi[Gd], beta, &(Buf2->matrix[h][ij][pq]), sospi[Gd]);
                    }
                } else {  // SO → MO transformation (backtransform)
                    X = block_matrix(sospi[Gc], mospi_right[Gd]);

                    for (int ij = 0; ij < Buf1->params->rowtot[h]; ij++) {
                        // Step 1: Transform right index (q) from SO to MO
                        // X(p,d) = (ij|p,q) × C2(q,d)
                        C_DGEMM('n', 'n', sospi[Gc], mospi_right[Gd], sospi[Gd], 1.0, &(Buf2->matrix[h][ij][pq]),
                                sospi[Gd], &(C2[Gd][0][0]), mospi_right[Gd], 0.0, &(X[0][0]), mospi_right[Gd]);

                        // Step 2: Transform left index (p) from SO to MO
                        // (ij|c,d) = C1(p,c)^T × X(p,d)
                        C_DGEMM('t', 'n', mospi_left[Gc], mospi_right[Gd], sospi[Gc], alpha, &(C1[Gc][0][0]),
                                mospi_left[Gc], &(X[0][0]), mospi_right[Gd], beta, &(Buf1->matrix[h][ij][cd]),
                                mospi_right[Gd]);
                    }
                }

                free_block(X);
            }
        }

        dpd_set_default(dpdnum1);
        if (type == 1) global_dpd_->buf4_mat_irrep_wrt(Buf1, h);
        global_dpd_->buf4_mat_irrep_close(Buf1, h);

        dpd_set_default(dpdnum2);
        if (type == 0) global_dpd_->buf4_mat_irrep_wrt(Buf2, h);
        global_dpd_->buf4_mat_irrep_close(Buf2, h);
    }
}

}  // namespace ccenergy
}  // namespace psi
