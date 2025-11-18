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

#include "sapt2p3.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
namespace sapt {

void SAPT2p3::exch_disp30() {
    auto tARBS = std::make_shared<Matrix>("tARBS", aoccA_ * nvirA_, aoccB_ * nvirB_);
    psio_->read_entry(PSIF_SAPT_AMPS, "Disp30 uARBS Amplitudes", (char *)tARBS->get_pointer(),
                      sizeof(double) * aoccA_ * nvirA_ * aoccB_ * nvirB_);

    auto vARBS = std::make_shared<Matrix>("vARBS", noccA_ * nvirA_, noccB_ * nvirB_);
    psio_->read_entry(PSIF_SAPT_AMPS, "Exch-Disp V_ARBS", (char *)vARBS->get_pointer(),
                      sizeof(double) * noccA_ * nvirA_ * noccB_ * nvirB_);
    double **vARBSp = vARBS->pointer();
    double **tARBSp = tARBS->pointer();

    double ex_1 = 0.0;

    for (int a = 0, ar = 0; a < aoccA_; a++) {
        for (int r = 0; r < nvirA_; r++, ar++) {
            int aarr = (a + foccA_) * nvirA_ + r;
            ex_1 -= 2.0 * C_DDOT(aoccB_ * nvirB_, &(vARBSp[aarr][foccB_ * nvirB_]), 1, tARBSp[ar], 1);
        }
    }

    // Automatic cleanup via shared_ptr

    double ex_2 = exch_disp30_20();
    double ex_3 = exch_disp30_02();
    double ex_4 = exch_disp30_22();

    e_exch_disp30_ = ex_1 + ex_2 + ex_3 + ex_4;

    if (debug_) {
        outfile->Printf("\n    Exch-Disp_1         = %18.12lf [Eh]\n", ex_1);
        outfile->Printf("    Exch-Disp_2         = %18.12lf [Eh]\n", ex_2);
        outfile->Printf("    Exch-Disp_3         = %18.12lf [Eh]\n", ex_3);
        outfile->Printf("    Exch-Disp_4         = %18.12lf [Eh]\n", ex_4);
    }
    if (print_) {
        outfile->Printf("    Exch-Disp30         = %18.12lf [Eh]\n", e_exch_disp30_);
    }
}

double SAPT2p3::exch_disp30_20() {
    double energy = 0.0;

    auto uARAR = std::make_shared<Matrix>("uARAR", aoccA_ * nvirA_, aoccA_ * nvirA_);
    double **uARARp = uARAR->pointer();
    double **B_p_AR = get_AR_ints(1, foccA_);
    auto T_p_AR = std::make_shared<Matrix>("T_p_AR", aoccA_ * nvirA_, ndf_ + 3);
    psio_->read_entry(PSIF_SAPT_AMPS, "T AR Intermediates", (char *)T_p_AR->get_pointer(),
                      sizeof(double) * aoccA_ * nvirA_ * (ndf_ + 3));

    C_DGEMM('N', 'T', aoccA_ * nvirA_, aoccA_ * nvirA_, ndf_ + 3, 1.0, &(B_p_AR[0][0]), ndf_ + 3, T_p_AR->get_pointer(),
            ndf_ + 3, 0.0, uARAR->get_pointer(), aoccA_ * nvirA_);

    // T_p_AR automatically cleaned up

    for (int ar = 0; ar < aoccA_ * nvirA_; ar++) {
        for (int a1r1 = 0; a1r1 < ar; a1r1++) {
            double tval = uARARp[ar][a1r1] + uARARp[a1r1][ar];
            uARARp[a1r1][ar] = tval;
            uARARp[ar][a1r1] = tval;
        }
    }

    C_DSCAL(aoccA_ * nvirA_, 2.0, uARAR->get_pointer(), aoccA_ * nvirA_ + 1);

    for (int a = 0, ar = 0; a < aoccA_; a++) {
        for (int r = 0; r < nvirA_; r++, ar++) {
            for (int aa = 0, aarr = 0; aa < aoccA_; aa++) {
                for (int rr = 0; rr < nvirA_; rr++, aarr++) {
                    double denom =
                        evalsA_[a + foccA_] + evalsA_[aa + foccA_] - evalsA_[r + noccA_] - evalsA_[rr + noccA_];
                    uARARp[ar][aarr] /= denom;
                }
            }
        }
    }

    auto U_p_AR = std::make_shared<Matrix>("U_p_AR", aoccA_ * nvirA_, ndf_ + 3);

    C_DGEMM('N', 'N', aoccA_ * nvirA_, ndf_ + 3, aoccA_ * nvirA_, 1.0, uARAR->get_pointer(), aoccA_ * nvirA_, &(B_p_AR[0][0]),
            ndf_ + 3, 0.0, U_p_AR->get_pointer(), ndf_ + 3);

    std::vector<double> X(nvirA_);

    for (int a = 0; a < aoccA_; a++) {
        for (int a1 = 0; a1 <= a; a1++) {
            for (int r = 0; r < nvirA_; r++) {
                int ar = a * nvirA_ + r;
                int a1r = a1 * nvirA_ + r;
                C_DCOPY(nvirA_, &(uARARp[ar][a1 * nvirA_]), 1, X.data(), 1);
                C_DCOPY(nvirA_, &(uARARp[a1r][a * nvirA_]), 1, &(uARARp[ar][a1 * nvirA_]), 1);
                C_DCOPY(nvirA_, X.data(), 1, &(uARARp[a1r][a * nvirA_]), 1);
            }
        }
    }

    // X automatically cleaned up

    auto U_p_ApR = std::make_shared<Matrix>("U_p_ApR", aoccA_ * nvirA_, ndf_ + 3);

    C_DGEMM('N', 'N', aoccA_ * nvirA_, ndf_ + 3, aoccA_ * nvirA_, 1.0, uARAR->get_pointer(), aoccA_ * nvirA_, &(B_p_AR[0][0]),
            ndf_ + 3, 0.0, U_p_ApR->get_pointer(), ndf_ + 3);

    free_block(B_p_AR);
    // uARAR automatically cleaned up

    double **B_p_RB = get_RB_ints(1);
    auto X_p_AR = std::make_shared<Matrix>("X_p_AR", aoccA_ * nvirA_, ndf_ + 3);
    double **X_p_ARp = X_p_AR->pointer();

    for (int r = 0; r < nvirA_; r++) {
        C_DGEMM('N', 'N', aoccA_, ndf_ + 3, noccB_, 1.0, &(sAB_[foccA_][0]), nmoB_, &(B_p_RB[r * noccB_][0]), ndf_ + 3,
                0.0, &(X_p_ARp[r][0]), nvirA_ * (ndf_ + 3));
    }

    energy = C_DDOT(aoccA_ * nvirA_ * (ndf_ + 3), U_p_ApR->get_pointer(), 1, X_p_AR->get_pointer(), 1);

    energy -= 2.0 * C_DDOT(aoccA_ * nvirA_ * (ndf_ + 3), U_p_AR->get_pointer(), 1, X_p_AR->get_pointer(), 1);

    free_block(B_p_RB);
    // X_p_AR automatically cleaned up

    auto xAR = std::make_shared<Matrix>("xAR", aoccA_, nvirA_);
    auto yAR = std::make_shared<Matrix>("yAR", aoccA_, nvirA_);

    C_DGEMM('N', 'T', aoccA_, nvirA_, noccB_, 1.0, &(sAB_[foccA_][0]), nmoB_, &(sAB_[noccA_][0]), nmoB_, 0.0,
            xAR->get_pointer(), nvirA_);

    C_DGEMV('n', aoccA_ * nvirA_, ndf_ + 3, 1.0, U_p_ApR->get_pointer(), ndf_ + 3, diagBB_, 1, 0.0, yAR->get_pointer(), 1);

    energy += 2.0 * C_DDOT(aoccA_ * nvirA_, xAR->get_pointer(), 1, yAR->get_pointer(), 1);

    C_DGEMV('n', aoccA_ * nvirA_, ndf_ + 3, 1.0, U_p_AR->get_pointer(), (ndf_ + 3), diagBB_, 1, 0.0, yAR->get_pointer(), 1);

    energy -= 4.0 * C_DDOT(aoccA_ * nvirA_, xAR->get_pointer(), 1, yAR->get_pointer(), 1);

    // xAR, yAR automatically cleaned up

    auto A_p_AB = std::make_shared<Matrix>("A_p_AB", aoccA_ * noccB_, ndf_ + 3);
    double **A_p_ABp = A_p_AB->pointer();
    auto A_p_BB = std::make_shared<Matrix>("A_p_BB", noccB_ * noccB_, ndf_ + 3);
    double **U_p_ApRp = U_p_ApR->pointer();
    double **U_p_ARp = U_p_AR->pointer();

    for (int a = 0; a < aoccA_; a++) {
        C_DGEMM('T', 'N', noccB_, ndf_ + 3, nvirA_, 1.0, &(sAB_[noccA_][0]), nmoB_, &(U_p_ApRp[a * nvirA_][0]), ndf_ + 3,
                0.0, &(A_p_ABp[a * noccB_][0]), ndf_ + 3);
    }

    C_DGEMM('T', 'N', noccB_, noccB_ * (ndf_ + 3), aoccA_, -1.0, &(sAB_[foccA_][0]), nmoB_, A_p_AB->get_pointer(),
            noccB_ * (ndf_ + 3), 0.0, A_p_BB->get_pointer(), noccB_ * (ndf_ + 3));

    for (int a = 0; a < aoccA_; a++) {
        C_DGEMM('T', 'N', noccB_, ndf_ + 3, nvirA_, 1.0, &(sAB_[noccA_][0]), nmoB_, &(U_p_ARp[a * nvirA_][0]), ndf_ + 3,
                0.0, &(A_p_ABp[a * noccB_][0]), ndf_ + 3);
    }

    C_DGEMM('T', 'N', noccB_, noccB_ * (ndf_ + 3), aoccA_, 2.0, &(sAB_[foccA_][0]), nmoB_, A_p_AB->get_pointer(),
            noccB_ * (ndf_ + 3), 1.0, A_p_BB->get_pointer(), noccB_ * (ndf_ + 3));

    double **B_p_BB = get_BB_ints(1);

    energy += C_DDOT(noccB_ * noccB_ * (ndf_ + 3), A_p_BB->get_pointer(), 1, &(B_p_BB[0][0]), 1);

    // A_p_AB, A_p_BB, U_p_AR, U_p_ApR automatically cleaned up
    free_block(B_p_BB);

    return (4.0 * energy);
}

double SAPT2p3::exch_disp30_02() {
    double energy = 0.0;

    auto uBSBS = std::make_shared<Matrix>("uBSBS", aoccB_ * nvirB_, aoccB_ * nvirB_);
    double **uBSBSp = uBSBS->pointer();
    double **B_p_BS = get_BS_ints(1, foccB_);
    auto T_p_BS = std::make_shared<Matrix>("T_p_BS", aoccB_ * nvirB_, ndf_ + 3);
    psio_->read_entry(PSIF_SAPT_AMPS, "T BS Intermediates", (char *)T_p_BS->get_pointer(),
                      sizeof(double) * aoccB_ * nvirB_ * (ndf_ + 3));

    C_DGEMM('N', 'T', aoccB_ * nvirB_, aoccB_ * nvirB_, ndf_ + 3, 1.0, &(B_p_BS[0][0]), ndf_ + 3, T_p_BS->get_pointer(),
            ndf_ + 3, 0.0, uBSBS->get_pointer(), aoccB_ * nvirB_);

    // T_p_BS automatically cleaned up

    for (int bs = 0; bs < aoccB_ * nvirB_; bs++) {
        for (int b1s1 = 0; b1s1 < bs; b1s1++) {
            double tval = uBSBSp[bs][b1s1] + uBSBSp[b1s1][bs];
            uBSBSp[b1s1][bs] = tval;
            uBSBSp[bs][b1s1] = tval;
        }
    }

    C_DSCAL(aoccB_ * nvirB_, 2.0, uBSBS->get_pointer(), aoccB_ * nvirB_ + 1);

    for (int b = 0, bs = 0; b < aoccB_; b++) {
        for (int s = 0; s < nvirB_; s++, bs++) {
            for (int bb = 0, bbss = 0; bb < aoccB_; bb++) {
                for (int ss = 0; ss < nvirB_; ss++, bbss++) {
                    double denom =
                        evalsB_[b + foccB_] + evalsB_[bb + foccB_] - evalsB_[s + noccB_] - evalsB_[ss + noccB_];
                    uBSBSp[bs][bbss] /= denom;
                }
            }
        }
    }

    auto U_p_BS = std::make_shared<Matrix>("U_p_BS", aoccB_ * nvirB_, ndf_ + 3);

    C_DGEMM('N', 'N', aoccB_ * nvirB_, (ndf_ + 3), aoccB_ * nvirB_, 1.0, uBSBS->get_pointer(), aoccB_ * nvirB_,
            &(B_p_BS[0][0]), ndf_ + 3, 0.0, U_p_BS->get_pointer(), ndf_ + 3);

    std::vector<double> X(nvirB_);

    for (int b = 0; b < aoccB_; b++) {
        for (int b1 = 0; b1 <= b; b1++) {
            for (int s = 0; s < nvirB_; s++) {
                int bs = b * nvirB_ + s;
                int b1s = b1 * nvirB_ + s;
                C_DCOPY(nvirB_, &(uBSBSp[bs][b1 * nvirB_]), 1, X.data(), 1);
                C_DCOPY(nvirB_, &(uBSBSp[b1s][b * nvirB_]), 1, &(uBSBSp[bs][b1 * nvirB_]), 1);
                C_DCOPY(nvirB_, X.data(), 1, &(uBSBSp[b1s][b * nvirB_]), 1);
            }
        }
    }

    // X automatically cleaned up

    auto U_p_BpS = std::make_shared<Matrix>("U_p_BpS", aoccB_ * nvirB_, ndf_ + 3);

    C_DGEMM('N', 'N', aoccB_ * nvirB_, ndf_ + 3, aoccB_ * nvirB_, 1.0, uBSBS->get_pointer(), aoccB_ * nvirB_, &(B_p_BS[0][0]),
            ndf_ + 3, 0.0, U_p_BpS->get_pointer(), ndf_ + 3);

    free_block(B_p_BS);
    // uBSBS automatically cleaned up

    double **B_p_AS = get_AS_ints(1);
    auto X_p_BS = std::make_shared<Matrix>("X_p_BS", aoccB_ * nvirB_, ndf_ + 3);

    C_DGEMM('T', 'N', aoccB_, nvirB_ * (ndf_ + 3), noccA_, 1.0, &(sAB_[0][foccB_]), nmoB_, &(B_p_AS[0][0]),
            nvirB_ * (ndf_ + 3), 0.0, X_p_BS->get_pointer(), nvirB_ * (ndf_ + 3));

    energy = C_DDOT(aoccB_ * nvirB_ * (ndf_ + 3), U_p_BpS->get_pointer(), 1, X_p_BS->get_pointer(), 1);

    energy -= 2.0 * C_DDOT(aoccB_ * nvirB_ * (ndf_ + 3), U_p_BS->get_pointer(), 1, X_p_BS->get_pointer(), 1);

    free_block(B_p_AS);
    // X_p_BS automatically cleaned up

    auto xBS = std::make_shared<Matrix>("xBS", aoccB_, nvirB_);
    auto yBS = std::make_shared<Matrix>("yBS", aoccB_, nvirB_);

    C_DGEMM('T', 'N', aoccB_, nvirB_, noccA_, 1.0, &(sAB_[0][foccB_]), nmoB_, &(sAB_[0][noccB_]), nmoB_, 0.0,
            xBS->get_pointer(), nvirB_);

    C_DGEMV('n', aoccB_ * nvirB_, ndf_ + 3, 1.0, U_p_BpS->get_pointer(), ndf_ + 3, diagAA_, 1, 0.0, yBS->get_pointer(), 1);

    energy += 2.0 * C_DDOT(aoccB_ * nvirB_, xBS->get_pointer(), 1, yBS->get_pointer(), 1);

    C_DGEMV('n', aoccB_ * nvirB_, ndf_ + 3, 1.0, U_p_BS->get_pointer(), ndf_ + 3, diagAA_, 1, 0.0, yBS->get_pointer(), 1);

    energy -= 4.0 * C_DDOT(aoccB_ * nvirB_, xBS->get_pointer(), 1, yBS->get_pointer(), 1);

    // xBS, yBS automatically cleaned up

    auto A_p_BA = std::make_shared<Matrix>("A_p_BA", aoccB_ * noccA_, ndf_ + 3);
    double **A_p_BAp = A_p_BA->pointer();
    auto A_p_AA = std::make_shared<Matrix>("A_p_AA", noccA_ * noccA_, ndf_ + 3);
    double **U_p_BpSp = U_p_BpS->pointer();
    double **U_p_BSp = U_p_BS->pointer();

    for (int b = 0; b < aoccB_; b++) {
        C_DGEMM('N', 'N', noccA_, (ndf_ + 3), nvirB_, 1.0, &(sAB_[0][noccB_]), nmoB_, &(U_p_BpSp[b * nvirB_][0]),
                (ndf_ + 3), 0.0, &(A_p_BAp[b * noccA_][0]), (ndf_ + 3));
    }

    C_DGEMM('N', 'N', noccA_, noccA_ * (ndf_ + 3), aoccB_, -1.0, &(sAB_[0][foccB_]), nmoB_, A_p_BA->get_pointer(),
            noccA_ * (ndf_ + 3), 0.0, A_p_AA->get_pointer(), noccA_ * (ndf_ + 3));

    for (int b = 0; b < aoccB_; b++) {
        C_DGEMM('N', 'N', noccA_, (ndf_ + 3), nvirB_, 1.0, &(sAB_[0][noccB_]), nmoB_, &(U_p_BSp[b * nvirB_][0]),
                (ndf_ + 3), 0.0, &(A_p_BAp[b * noccA_][0]), (ndf_ + 3));
    }

    C_DGEMM('N', 'N', noccA_, noccA_ * (ndf_ + 3), aoccB_, 2.0, &(sAB_[0][foccB_]), nmoB_, A_p_BA->get_pointer(),
            noccA_ * (ndf_ + 3), 1.0, A_p_AA->get_pointer(), noccA_ * (ndf_ + 3));

    double **B_p_AA = get_AA_ints(1);

    energy += C_DDOT(noccA_ * noccA_ * (ndf_ + 3), A_p_AA->get_pointer(), 1, &(B_p_AA[0][0]), 1);

    // A_p_BA, A_p_AA, U_p_BS, U_p_BpS automatically cleaned up
    free_block(B_p_AA);

    return (4.0 * energy);
}

double SAPT2p3::exch_disp30_22() {
    double energy = 0.0;

    auto tARBS = std::make_shared<Matrix>("tARBS", aoccA_ * nvirA_, aoccB_ * nvirB_);
    psio_->read_entry(PSIF_SAPT_AMPS, "tARBS Amplitudes", (char *)tARBS->get_pointer(),
                      sizeof(double) * aoccA_ * nvirA_ * aoccB_ * nvirB_);
    double **tARBSp = tARBS->pointer();

    auto tAS_RB = std::make_shared<Matrix>("tAS_RB", nvirA_, aoccB_);
    double **tAS_RBp = tAS_RB->pointer();
    auto tRB_AS = std::make_shared<Matrix>("tRB_AS", aoccA_, nvirB_);
    double **tRB_ASp = tRB_AS->pointer();

    for (int a = 0, ar = 0; a < aoccA_; a++) {
        for (int r = 0; r < nvirA_; r++, ar++) {
            for (int b = 0, bs = 0; b < aoccB_; b++) {
                for (int s = 0; s < nvirB_; s++, bs++) {
                    tAS_RBp[r][b] += tARBSp[ar][bs] * sAB_[a + foccA_][s + noccB_];
                    tRB_ASp[a][s] += tARBSp[ar][bs] * sAB_[r + noccA_][b + foccB_];
                }
            }
        }
    }

    auto T_p_AR = std::make_shared<Matrix>("T_p_AR", aoccA_ * nvirA_, ndf_ + 3);
    psio_->read_entry(PSIF_SAPT_AMPS, "T AR Intermediates", (char *)T_p_AR->get_pointer(),
                      sizeof(double) * aoccA_ * nvirA_ * (ndf_ + 3));
    double **T_p_ARp = T_p_AR->pointer();

    auto T_p_BS = std::make_shared<Matrix>("T_p_BS", aoccB_ * nvirB_, ndf_ + 3);
    psio_->read_entry(PSIF_SAPT_AMPS, "T BS Intermediates", (char *)T_p_BS->get_pointer(),
                      sizeof(double) * aoccB_ * nvirB_ * (ndf_ + 3));
    double **T_p_BSp = T_p_BS->pointer();

    double **B_p_AR = get_AR_ints(0, foccA_);
    double **B_p_BS = get_BS_ints(0, foccB_);

    auto X_p_BS = std::make_shared<Matrix>("X_p_BS", aoccB_ * nvirB_, ndf_ + 3);
    double **X_p_BSp = X_p_BS->pointer();

    auto xBB = std::make_shared<Matrix>("xBB", aoccB_, aoccB_);
    auto xSS = std::make_shared<Matrix>("xSS", nvirB_, nvirB_);

    C_DGEMM('T', 'N', aoccB_, aoccB_, nvirA_, 1.0, &(sAB_[noccA_][foccB_]), nmoB_, tAS_RBp[0], aoccB_, 0.0,
            xBB->get_pointer(), aoccB_);

    C_DGEMM('N', 'N', aoccB_, nvirB_ * (ndf_ + 3), aoccB_, 2.0, xBB->get_pointer(), aoccB_, &(B_p_BS[0][0]),
            nvirB_ * (ndf_ + 3), 0.0, X_p_BSp[0], nvirB_ * (ndf_ + 3));

    C_DGEMM('T', 'N', nvirB_, nvirB_, aoccA_, 1.0, &(sAB_[foccA_][noccB_]), nmoB_, tRB_ASp[0], nvirB_, 0.0,
            xSS->get_pointer(), nvirB_);

    for (int b = 0; b < aoccB_; b++) {
        C_DGEMM('N', 'N', nvirB_, ndf_ + 3, nvirB_, 2.0, xSS->get_pointer(), nvirB_, &(B_p_BS[b * nvirB_][0]), ndf_ + 3, 1.0,
                &(X_p_BSp[b * nvirB_][0]), ndf_ + 3);
    }

    energy += C_DDOT(aoccB_ * nvirB_ * (ndf_ + 3), T_p_BSp[0], 1, X_p_BSp[0], 1);

    // xBB, xSS, X_p_BS automatically cleaned up

    auto X_p_AR = std::make_shared<Matrix>("X_p_AR", aoccA_ * nvirA_, ndf_ + 3);
    double **X_p_ARp = X_p_AR->pointer();

    auto xAA = std::make_shared<Matrix>("xAA", aoccA_, aoccA_);
    auto xRR = std::make_shared<Matrix>("xRR", nvirA_, nvirA_);

    C_DGEMM('N', 'T', aoccA_, aoccA_, nvirB_, 1.0, &(sAB_[foccA_][noccB_]), nmoB_, tRB_ASp[0], nvirB_, 0.0,
            xAA->get_pointer(), aoccA_);

    C_DGEMM('N', 'N', aoccA_, nvirA_ * (ndf_ + 3), aoccA_, 2.0, xAA->get_pointer(), aoccA_, &(B_p_AR[0][0]),
            nvirA_ * (ndf_ + 3), 0.0, X_p_ARp[0], nvirA_ * (ndf_ + 3));

    C_DGEMM('N', 'T', nvirA_, nvirA_, aoccB_, 1.0, &(sAB_[noccA_][foccB_]), nmoB_, tAS_RBp[0], aoccB_, 0.0,
            xRR->get_pointer(), nvirA_);

    for (int a = 0; a < aoccA_; a++) {
        C_DGEMM('N', 'N', nvirA_, ndf_ + 3, nvirA_, 2.0, xRR->get_pointer(), nvirA_, &(B_p_AR[a * nvirA_][0]), ndf_ + 3, 1.0,
                &(X_p_ARp[a * nvirA_][0]), ndf_ + 3);
    }

    energy += C_DDOT(aoccA_ * nvirA_ * (ndf_ + 3), T_p_ARp[0], 1, X_p_ARp[0], 1);

    // xAA, xRR, X_p_AR automatically cleaned up

    auto A_p_AB = std::make_shared<Matrix>("A_p_AB", aoccA_ * aoccB_, ndf_ + 3);
    double **A_p_ABp = A_p_AB->pointer();
    auto B_p_AB = std::make_shared<Matrix>("B_p_AB", aoccA_ * aoccB_, ndf_ + 3);
    double **B_p_ABp = B_p_AB->pointer();

    for (int a = 0; a < aoccA_; a++) {
        C_DGEMM('T', 'N', aoccB_, ndf_ + 3, nvirA_, 1.0, &(sAB_[noccA_][foccB_]), nmoB_, &(T_p_ARp[a * nvirA_][0]),
                ndf_ + 3, 0.0, &(A_p_ABp[a * aoccB_][0]), ndf_ + 3);
    }

    for (int b = 0; b < aoccB_; b++) {
        C_DGEMM('N', 'N', aoccA_, ndf_ + 3, nvirB_, 1.0, &(sAB_[foccA_][noccB_]), nmoB_, &(T_p_BSp[b * nvirB_][0]),
                ndf_ + 3, 0.0, &(B_p_ABp[b][0]), aoccB_ * (ndf_ + 3));
    }

    energy -= 4.0 * C_DDOT(aoccA_ * aoccB_ * (ndf_ + 3), A_p_ABp[0], 1, B_p_ABp[0], 1);

    for (int a = 0; a < aoccA_; a++) {
        C_DGEMM('T', 'N', aoccB_, ndf_ + 3, nvirA_, 1.0, tAS_RBp[0], aoccB_, &(B_p_AR[a * nvirA_][0]), ndf_ + 3,
                0.0, &(A_p_ABp[a * aoccB_][0]), ndf_ + 3);
    }

    for (int b = 0; b < aoccB_; b++) {
        C_DGEMM('N', 'N', aoccA_, ndf_ + 3, nvirB_, 1.0, tRB_ASp[0], nvirB_, &(B_p_BS[b * nvirB_][0]), ndf_ + 3,
                0.0, &(B_p_ABp[b][0]), aoccB_ * (ndf_ + 3));
    }

    energy -= C_DDOT(aoccA_ * aoccB_ * (ndf_ + 3), A_p_ABp[0], 1, B_p_ABp[0], 1);

    // A_p_AB, B_p_AB, T_p_AR, T_p_BS automatically cleaned up

    auto tABRS = std::make_shared<Matrix>("tABRS", aoccA_ * aoccB_, nvirA_ * nvirB_);
    double **tABRSp = tABRS->pointer();

    for (int a = 0, ar = 0; a < aoccA_; a++) {
        for (int r = 0; r < nvirA_; r++, ar++) {
            for (int b = 0; b < aoccB_; b++) {
                int ab = a * aoccB_ + b;
                C_DCOPY(nvirB_, &(tARBSp[ar][b * nvirB_]), 1, &(tABRSp[ab][r * nvirB_]), 1);
            }
        }
    }

    // tARBS automatically cleaned up

    auto vABRS = std::make_shared<Matrix>("vABRS", aoccA_ * aoccB_, nvirA_ * nvirB_);
    double **vABRSp = vABRS->pointer();

    for (int a = 0, ab = 0; a < aoccA_; a++) {
        for (int b = 0; b < aoccB_; b++, ab++) {
            C_DGEMM('N', 'T', nvirA_, nvirB_, ndf_ + 3, 1.0, &(B_p_AR[a * nvirA_][0]), ndf_ + 3,
                    &(B_p_BS[b * nvirB_][0]), ndf_ + 3, 0.0, &(vABRSp[ab][0]), nvirB_);
        }
    }

    free_block(B_p_AR);
    free_block(B_p_BS);

    auto xAR = std::make_shared<Matrix>("xAR", aoccA_, nvirA_);
    auto ABAB = std::make_shared<Matrix>("ABAB", aoccA_ * aoccB_, aoccA_ * aoccB_);

    for (int a = 0, ab = 0; a < aoccA_; a++) {
        for (int b = 0; b < aoccB_; b++, ab++) {
            C_DGEMM('N', 'T', aoccA_, nvirA_, nvirB_, 1.0, &(sAB_[foccA_][noccB_]), nmoB_, &(tABRSp[ab][0]), nvirB_, 0.0,
                    xAR->get_pointer(), nvirA_);
            C_DGEMM('N', 'N', aoccA_, aoccB_, nvirA_, 1.0, xAR->get_pointer(), nvirA_, &(sAB_[noccA_][foccB_]), nmoB_, 0.0,
                    ABAB->get_pointer(ab), aoccB_);
        }
    }

    // xAR automatically cleaned up

    auto xABRS = std::make_shared<Matrix>("xABRS", aoccA_ * aoccB_, nvirA_ * nvirB_);

    C_DGEMM('T', 'N', aoccA_ * aoccB_, nvirA_ * nvirB_, aoccA_ * aoccB_, 1.0, ABAB->get_pointer(), aoccA_ * aoccB_,
            vABRSp[0], nvirA_ * nvirB_, 0.0, xABRS->get_pointer(), nvirA_ * nvirB_);

    energy -= C_DDOT((long int)aoccA_ * aoccB_ * nvirA_ * nvirB_, tABRSp[0], 1, xABRS->get_pointer(), 1);

    // tABRS, ABAB, vABRS, xABRS, tAS_RB, tRB_AS automatically cleaned up

    return (2.0 * energy);
}
}  // namespace sapt
}  // namespace psi
