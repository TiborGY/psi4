/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <sstream>
#include "psi4/libciomr/libciomr.h"
#include "dpd.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"

namespace psi {

/* Invariant check for DPD pairorb int-matrices, run at teardown.
 *
 * After DPD::init() completes, every non-null pairorb[pair][h] (an
 * int** of shape pairtot[pair][h] x 2) contains only non-negative
 * absolute orbital indices: the -1 sentinel values written during
 * allocation (init.cc:363-367 and 516-523) are unconditionally
 * overwritten by the fill loops (init.cc:382-487 and 536-554), which
 * visit each (p, q) pair exactly once via the monotonically
 * incremented count[h] counter. No consumer of libdpd writes to
 * pairorb (or to roworb/colorb, which alias pairorb) after init,
 * so the non-negativity invariant holds from init to teardown.
 *
 * If, at teardown, any entry is still negative or has become
 * garbage, one of two things has happened:
 *   (a) DPD::init() never finished filling every slot (an internal
 *       libdpd bug), or
 *   (b) a consumer of libdpd wrote out of bounds and overwrote a
 *       region of memory belonging to a pairorb int-matrix. The
 *       stability4 rhf.cc OOB crash was exactly this case: a
 *       negative column index in a Matrix::set() call corrupted
 *       the heap metadata of an adjacent chunk, surfacing as a bare
 *       glibc "free(): invalid next size" abort in free_int_matrix
 *       (close.cc) with no indication of the real culprit. A
 *       different allocator layout, or a larger OOB, would have
 *       clobbered DPD's own pairorb data directly, and this check
 *       catches that variant.
 *
 * Running this scan always-on (rather than under #ifdef DPD_DEBUG)
 * honors the request that libdpd ergonomics turn this class of
 * bug into a thrown PSIEXCEPTION instead of a bare abort in user
 * builds. The cost is O(sum pairtot[pair][h] * 2) integer reads,
 * executed once per dpd_close() -- negligible against the dozens
 * of free() calls it accompanies and against any real SCF/CC run.
 *
 * Throwing from a destructor is generally discouraged. We guard
 * with std::uncaught_exceptions(): if a stack-unwind is already
 * in flight (e.g. dpd_close was called from ~IntegralTransform
 * during exception propagation of some other psi4 error), throwing
 * here would call std::terminate, losing both diagnostics. In that
 * case we print the corruption message to outfile and stderr and
 * fall through to the existing free path (which glibc will then
 * abort with its own message). The non-unwind case throws cleanly
 * so that psi4's normal Python-level error reporting surfaces it
 * as a PsiException.
 *
 * Note: we deliberately do NOT call dpd_error() here -- we are
 * already inside ~DPD(), and dpd_error() calls dpd_close() which
 * calls delete on this again, which would recurse. A plain throw
 * is the safe choice.
 */
static void dpd_pairorb_integrity_check(int ****pairorb, int **pairtot, int num_pairs, int nirreps) {
    for (int pair = 0; pair < num_pairs; ++pair) {
        if (!pairorb[pair]) continue;
        for (int h = 0; h < nirreps; ++h) {
            int rows = pairtot[pair][h];
            if (rows <= 0 || !pairorb[pair][h]) continue;
            for (int r = 0; r < rows; ++r) {
                int p = pairorb[pair][h][r][0];
                int q = pairorb[pair][h][r][1];
                if (p < 0 || q < 0) {
                    std::ostringstream oss;
                    oss << "DPD pairorb integrity check failed at teardown: pair=" << pair << " irrep=" << h
                        << " row=" << r << " orb=(" << p << "," << q << "). A consumer wrote out of bounds into a "
                        << "libdpd pairorb int-matrix (or DPD::init never finished). Run under AddressSanitizer to "
                        << "locate the offending write.";
                    const std::string msg = oss.str();
                    if (std::uncaught_exceptions() == 0) {
                        throw PSIEXCEPTION(msg);
                    } else {
                        /* Already unwinding -- throwing here would call
                         * std::terminate. Print loudly and let the existing
                         * free path surface glibc's abort, which is still
                         * (slightly) more legible than a silent terminate. */
                        if (outfile) outfile->Printf("%s\n", msg.c_str());
                        std::fprintf(stderr, "%s\n", msg.c_str());
                        return;
                    }
                }
            }
        }
    }
}

DPD::~DPD() {
    int h, i, j, k, cnt;
    /*  dpd_file2_cache_print(stdout); */
    file2_cache_close();
    /*  dpd_file4_cache_print(stdout);*/
    file4_cache_close();

    if (params4)
        for (i = 0; i < num_pairs; i++)
            for (j = 0; j < num_pairs; j++) free_int_matrix(params4[i][j].start13);

    if (orboff) {
        for (i = 0; i < num_subspaces; i++) free(orboff[i]);
        free(orboff);
    }

    if (pairidx && pairorb) {
        /* Verify pairorb int-matrices were not corrupted between
         * dpd_init and teardown. See comment above
         * dpd_pairorb_integrity_check(). */
        dpd_pairorb_integrity_check(pairorb, pairtot, num_pairs, nirreps);
        for (i = 0; i < num_subspaces; i++) {
            for (j = 0; j < 5; j++) {
                free_int_matrix(pairidx[5 * i + j]);
                for (k = 0; k < nirreps; k++)
                    if (pairtot[5 * i + j][k]) free_int_matrix(pairorb[5 * i + j][k]);
                free(pairorb[5 * i + j]);
            }
        }
        for (i = 0, cnt = 5 * num_subspaces; i < num_subspaces; i++) {
            for (j = i + 1; j < num_subspaces; j++, cnt += 2) {
                free_int_matrix(pairidx[cnt]);
                free_int_matrix(pairidx[cnt + 1]);
                for (k = 0; k < nirreps; k++) {
                    if (pairtot[cnt][k]) free_int_matrix(pairorb[cnt][k]);
                    if (pairtot[cnt + 1][k]) free_int_matrix(pairorb[cnt + 1][k]);
                }
                free(pairorb[cnt]);
                free(pairorb[cnt + 1]);
            }
        }
        free(pairidx);
        free(pairorb);
    }

    if (orbs2 && orbidx2) {
        for (i = 0; i < num_subspaces; i++) {
            free(orbidx2[i]);
            for (j = 0; j < nirreps; j++) {
                if (orbspi[i][j]) free(orbs2[i][j]);
            }
            free(orbs2[i]);
        }
        free(orbidx2);
        free(orbs2);
    }

    if (orbspi && orbsym) {
        for (i = 0; i < num_subspaces; i++) {
            free(orbspi[i]);
            free(orbsym[i]);
        }
        free(orbspi);
        free(orbsym);
    }

    if (pairtot) free_int_matrix(pairtot);

    if (numorbs) free(numorbs);

    if (params4) {
        for (i = 0; i < num_pairs; i++) free(params4[i]);
        free(params4);
    }
    if (params2) {
        for (i = 0; i < num_subspaces; i++) free(params2[i]);
        free(params2);
    }

    /*
    printf("memory = %d; memfree = %d\n",
    dpd_main.memory, dpd_main.memfree);
  */
}

}  // namespace psi
