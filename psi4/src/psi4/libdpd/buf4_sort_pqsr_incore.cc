/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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
    \brief Incore pqsr -> pqrs sort, with and without adding the result to a target dpdbuf4 that already exists.
    The if constexpr results in a compile-time if, auto-generating the two different version of this function.
*/

#include "dpd.h"

namespace psi {
template <bool axpy>
void DPD::buf4_sort_pqsr_incore(dpdbuf4 &OutBuf, const dpdbuf4 &InBuf, const int32_t nirreps, const int32_t my_irrep,
                                const double alpha /*=1.0*/) {
    for (int32_t h = 0; h < nirreps; h++) {
        const auto r_irrep = h ^ my_irrep;

        for (int64_t pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
            const auto p = OutBuf.params->roworb[h][pq][0];
            const auto q = OutBuf.params->roworb[h][pq][1];
            const auto row = InBuf.params->rowidx[p][q];

            for (int64_t rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                const auto r = OutBuf.params->colorb[r_irrep][rs][0];
                const auto s = OutBuf.params->colorb[r_irrep][rs][1];
                const auto sr = InBuf.params->colidx[s][r];
                if constexpr (axpy) {
                    OutBuf.matrix[h][pq][rs] += alpha * InBuf.matrix[h][row][sr];
                } else {
                    OutBuf.matrix[h][pq][rs] = InBuf.matrix[h][row][sr];
                }
            }
        }
    }
}
}  // namespace psi
