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

#include "sapt2.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
namespace sapt {

void SAPT2::elst12() {
    double e_elst120 = elst120(wBAA_, wBRR_, CHFA_, PSIF_SAPT_AMPS, "pAA Density Matrix", "pRR Density Matrix",
                               "Y2 AR Amplitudes", foccA_, noccA_, nvirA_);

    if (debug_) {
        outfile->Printf("    Elst120,r           = %18.12lf [Eh]\n", e_elst120);
    }

    double e_elst102 = elst120(wABB_, wASS_, CHFB_, PSIF_SAPT_AMPS, "pBB Density Matrix", "pSS Density Matrix",
                               "Y2 BS Amplitudes", foccB_, noccB_, nvirB_);

    if (debug_) {
        outfile->Printf("    Elst102,r           = %18.12lf [Eh]\n\n", e_elst102);
    }

    e_elst12_ = e_elst120 + e_elst102;

    if (print_) {
        outfile->Printf("    Elst12,r            = %18.12lf [Eh]\n", e_elst12_);
    }
}

double SAPT2::elst120(double **wBAA, double **wBRR, double **CHFA, int ampfile, const char *pAAlabel,
                      const char *pRRlabel, const char *Ylabel, size_t foccA, size_t noccA, size_t nvirA) {
    int aoccA = noccA - foccA;

    auto pAA = std::make_shared<Matrix>("pAA", aoccA, aoccA);
    psio_->read_entry(ampfile, pAAlabel, (char *)pAA->get_pointer(), sizeof(double) * aoccA * aoccA);

    auto pRR = std::make_shared<Matrix>("pRR", nvirA, nvirA);
    psio_->read_entry(ampfile, pRRlabel, (char *)pRR->get_pointer(), sizeof(double) * nvirA * nvirA);

    auto yAR = std::make_shared<Matrix>("yAR", aoccA, nvirA);
    psio_->read_entry(ampfile, Ylabel, (char *)yAR->get_pointer(), sizeof(double) * aoccA * nvirA);

    double e1 = 0.0, e2 = 0.0, e3 = 0.0;

    double **pAAp = pAA->pointer();
    for (int a = 0; a < aoccA; a++) {
        e1 -= 2.0 * C_DDOT(aoccA, pAAp[a], 1, &(wBAA[a + foccA][foccA]), 1);
    }

    e2 += 2.0 * C_DDOT(nvirA * nvirA, pRR->get_pointer(), 1, wBRR[0], 1);
    e3 += 4.0 * C_DDOT(aoccA * nvirA, yAR->get_pointer(), 1, CHFA[foccA], 1);

    if (debug_) {
        outfile->Printf("\n    Elst12_1            = %18.12lf [Eh]\n", e1);
        outfile->Printf("    Elst12_2            = %18.12lf [Eh]\n", e2);
        outfile->Printf("    Elst12_3            = %18.12lf [Eh]\n", e3);
    }

    return (e1 + e2 + e3);
}
}  // namespace sapt
}  // namespace psi
