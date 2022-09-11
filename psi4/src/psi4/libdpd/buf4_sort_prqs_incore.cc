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
    \brief Incore pqrs <- prqs sort, with and without adding the result to a target dpdbuf4 that already exists.
    The if constexpr results in a compile-time if, auto-generating the two different version of this function.
*/

#include "dpd.h"

namespace psi {
template <bool axpy>
void DPD::buf4_sort_prqs_incore(dpdbuf4 &OutBuf, const dpdbuf4 &InBuf, const int32_t nirreps, const int32_t my_irrep,
                                const double alpha /*=1.0*/) {
    for (int32_t h = 0; h < nirreps; h++) {
        const auto r_irrep = h ^ my_irrep;

        for (int32_t Gp = 0; Gp < nirreps; Gp++) {
            const auto Gq = Gp ^ h;

            for (int32_t Gr = 0; Gr < nirreps; Gr++) {
                const auto Gs = Gr ^ r_irrep;
                // Irreps on the source
                const auto Gpr = Gp ^ Gr;
                const auto Gqs = Gq ^ Gs;

                for (int64_t p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                    const auto P = OutBuf.params->poff[Gp] + p;

                    for (int64_t q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                        const auto Q = OutBuf.params->qoff[Gq] + q;
                        const auto pq = OutBuf.params->rowidx[P][Q];

                        for (int64_t r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                            const auto R = OutBuf.params->roff[Gr] + r;
                            const auto pr = InBuf.params->rowidx[P][R];

                            for (int64_t s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                const auto S = OutBuf.params->soff[Gs] + s;
                                const auto rs = OutBuf.params->colidx[R][S];
                                const auto qs = InBuf.params->colidx[Q][S];
                                if constexpr (axpy) {
                                    OutBuf.matrix[h][pq][rs] += alpha * InBuf.matrix[Gpr][pr][qs];
                                } else {
                                    OutBuf.matrix[h][pq][rs] = InBuf.matrix[Gpr][pr][qs];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
}  // namespace psi
