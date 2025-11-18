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
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
namespace sapt {

void SAPT2p3::disp30() {
    if (third_order_) {
        double **B_p_AR = get_DF_ints(PSIF_SAPT_AA_DF_INTS, "AR RI Integrals", foccA_, noccA_, 0, nvirA_);
        double **B_p_BS = get_DF_ints(PSIF_SAPT_BB_DF_INTS, "BS RI Integrals", foccB_, noccB_, 0, nvirB_);

        auto tARBS = std::make_shared<Matrix>("tARBS", aoccA_ * nvirA_, aoccB_ * nvirB_);
        auto vARBS = std::make_shared<Matrix>("vARBS", aoccA_ * nvirA_, aoccB_ * nvirB_);
        psio_->read_entry(PSIF_SAPT_AMPS, "Disp30 uARBS Amplitudes", (char *)tARBS->get_pointer(),
                          sizeof(double) * aoccA_ * nvirA_ * aoccB_ * nvirB_);

        C_DGEMM('N', 'T', aoccA_ * nvirA_, aoccB_ * nvirB_, ndf_ + 3, 1.0, B_p_AR[0], ndf_ + 3, B_p_BS[0], ndf_ + 3,
                0.0, vARBS->get_pointer(), aoccB_ * nvirB_);

        e_disp30_ = 4.0 * C_DDOT((long int)aoccA_ * nvirA_ * aoccB_ * nvirB_, vARBS->get_pointer(), 1, tARBS->get_pointer(), 1);

        free_block(B_p_AR);
        free_block(B_p_BS);
    } else {
        double e1 = disp30_1(PSIF_SAPT_AMPS, "tARBS Amplitudes", PSIF_SAPT_AA_DF_INTS, "RR RI Integrals",
                             PSIF_SAPT_BB_DF_INTS, "SS RI Integrals", foccA_, noccA_, nvirA_, foccB_, noccB_, nvirB_);
        double e2 = disp30_2(PSIF_SAPT_AMPS, "tARBS Amplitudes", PSIF_SAPT_AA_DF_INTS, "AA RI Integrals",
                             "RR RI Integrals", PSIF_SAPT_BB_DF_INTS, "BB RI Integrals", "SS RI Integrals", foccA_,
                             noccA_, nvirA_, foccB_, noccB_, nvirB_);

        e_disp30_ = e1 + e2;
    }

    if (print_) {
        outfile->Printf("    Disp30              = %18.12lf [Eh]\n", e_disp30_);
    }
}

double SAPT2p3::disp30_1(int ampfile, const char *amplabel, int AAintfile, const char *RRlabel, int BBintfile,
                         const char *SSlabel, size_t foccA, size_t noccA, size_t nvirA, size_t foccB, size_t noccB, size_t nvirB) {
    size_t aoccA = noccA - foccA;
    size_t aoccB = noccB - foccB;

    auto tARBS = std::make_shared<Matrix>("tARBS", aoccA_ * nvirA_, aoccB_ * nvirB_);
    psio_->read_entry(PSIF_SAPT_AMPS, "tARBS Amplitudes", (char *)tARBS->get_pointer(),
                      sizeof(double) * aoccA_ * nvirA_ * aoccB_ * nvirB_);
    double **tARBSp = tARBS->pointer();

    auto tRSAB = std::make_shared<Matrix>("tRSAB", nvirA * nvirB, aoccA * aoccB);
    double **tRSABp = tRSAB->pointer();

    for (int a = 0, ar = 0; a < aoccA; a++) {
        for (int r = 0; r < nvirA; r++, ar++) {
            for (int b = 0, bs = 0; b < aoccB; b++) {
                for (int s = 0; s < nvirB; s++, bs++) {
                    int ab = a * aoccB + b;
                    int rs = r * nvirB + s;
                    int sr = s * nvirA + r;
                    tRSABp[rs][ab] = tARBSp[ar][bs];
                }
            }
        }
    }

    double energy = 0.0;

    psio_address next_DF_RR = PSIO_ZERO;
    psio_address next_DF_SS = PSIO_ZERO;

    auto B_p_RR = std::make_shared<Matrix>("B_p_RR", nvirA * (nvirA + 1) / 2, ndf_ + 3);
    double **B_p_RRp = B_p_RR->pointer();
    auto B_p_SS = std::make_shared<Matrix>("B_p_SS", nvirB * (nvirB + 1) / 2, ndf_ + 3);
    double **B_p_SSp = B_p_SS->pointer();

    for (int r1 = 0, r1r2 = 0; r1 < nvirA; r1++) {
        for (int r2 = 0; r2 <= r1; r2++, r1r2++) {
            next_DF_RR = psio_get_address(PSIO_ZERO, sizeof(double) * (r1 * nvirA + r2) * (ndf_ + 3));
            psio_->read(AAintfile, RRlabel, (char *)&(B_p_RRp[r1r2][0]), sizeof(double) * (ndf_ + 3), next_DF_RR,
                        &next_DF_RR);
            if (r1 != r2) C_DSCAL(ndf_ + 3, 2.0, B_p_RRp[r1r2], 1);
        }
    }

    for (int s1 = 0, s1s2 = 0; s1 < nvirB; s1++) {
        for (int s2 = 0; s2 <= s1; s2++, s1s2++) {
            next_DF_SS = psio_get_address(PSIO_ZERO, sizeof(double) * (s1 * nvirB + s2) * (ndf_ + 3));
            psio_->read(BBintfile, SSlabel, (char *)&(B_p_SSp[s1s2][0]), sizeof(double) * (ndf_ + 3), next_DF_SS,
                        &next_DF_SS);
            if (s1 != s2) C_DSCAL(ndf_ + 3, 2.0, B_p_SSp[s1s2], 1);
        }
    }

    auto xRS = std::make_shared<Matrix>("xRS", nvirA, nvirB * nvirB);
    double **xRSp = xRS->pointer();
    auto yRS = std::make_shared<Matrix>("yRS", nvirA, nvirB * (nvirB + 1) / 2);
    double **yRSp = yRS->pointer();
    double *zSS = init_array(nvirB * (nvirB + 1) / 2);

    for (int r1 = 0; r1 < nvirA; r1++) {
        C_DGEMM('N', 'T', (r1 + 1) * nvirB, nvirB, aoccA * aoccB, 1.0, tRSAB->get_pointer(), aoccA * aoccB, tRSABp[r1 * nvirB],
                aoccA * aoccB, 0.0, xRS->get_pointer(), nvirB);
        C_DGEMM('N', 'T', (r1 + 1), nvirB * (nvirB + 1) / 2, ndf_ + 3, 1.0, B_p_RRp[ioff_[r1]], ndf_ + 3, B_p_SS->get_pointer(),
                ndf_ + 3, 0.0, yRS->get_pointer(), nvirB * (nvirB + 1) / 2);
        for (int r2 = 0; r2 <= r1; r2++) {
            for (int s1 = 0, s1s2 = 0; s1 < nvirB; s1++) {
                for (int s2 = 0; s2 <= s1; s2++, s1s2++) {
                    zSS[s1s2] = xRSp[r2][s1 * nvirB + s2];
                    zSS[s1s2] += xRSp[r2][s2 * nvirB + s1];
                }
            }
            energy += 2.0 * C_DDOT(nvirB * (nvirB + 1) / 2, zSS, 1, yRSp[r2], 1);
        }
    }

    free(zSS);

    return (energy);
}

double SAPT2p3::disp30_2(int ampfile, const char *amplabel, int AAintfile, const char *AAlabel, const char *RRlabel,
                         int BBintfile, const char *BBlabel, const char *SSlabel, size_t foccA, size_t noccA, size_t nvirA,
                         size_t foccB, size_t noccB, size_t nvirB) {
    size_t aoccA = noccA - foccA;
    size_t aoccB = noccB - foccB;

    auto tARBS = std::make_shared<Matrix>("tARBS", aoccA_ * nvirA_, aoccB_ * nvirB_);
    psio_->read_entry(PSIF_SAPT_AMPS, "tARBS Amplitudes", (char *)tARBS->get_pointer(),
                      sizeof(double) * aoccA_ * nvirA_ * aoccB_ * nvirB_);
    double **tARBSp = tARBS->pointer();

    auto tABRS = std::make_shared<Matrix>("tABRS", aoccA * aoccB, nvirA * nvirB);
    double **tABRSp = tABRS->pointer();

    for (int a = 0, ar = 0; a < aoccA; a++) {
        for (int r = 0; r < nvirA; r++, ar++) {
            for (int b = 0, bs = 0; b < aoccB; b++) {
                for (int s = 0; s < nvirB; s++, bs++) {
                    int ab = a * aoccB + b;
                    int rs = r * nvirB + s;
                    tABRSp[ab][rs] = tARBSp[ar][bs];
                }
            }
        }
    }

    auto t2ABRS = std::make_shared<Matrix>("t2ABRS", aoccA * aoccB, nvirA * nvirB);

    double **B_p_AA = get_DF_ints(AAintfile, AAlabel, foccA, noccA, foccA, noccA);
    double **B_p_BB = get_DF_ints(BBintfile, BBlabel, foccB, noccB, foccB, noccB);

    auto ABAB = std::make_shared<Matrix>("ABAB", aoccA * aoccB, aoccA * aoccB);
    double **ABABp = ABAB->pointer();

    for (int a = 0, ab = 0; a < aoccA; a++) {
        for (int b = 0; b < aoccB; b++, ab++) {
            C_DGEMM('N', 'T', aoccA, aoccB, ndf_ + 3, 1.0, &(B_p_AA[a * aoccA][0]), ndf_ + 3, &(B_p_BB[b * aoccB][0]),
                    ndf_ + 3, 0.0, &(ABABp[ab][0]), aoccB);
        }
    }

    free_block(B_p_AA);
    free_block(B_p_BB);

    C_DGEMM('N', 'N', aoccA * aoccB, nvirA * nvirB, aoccA * aoccB, 1.0, ABAB->get_pointer(), aoccA * aoccB, tABRS->get_pointer(),
            nvirA * nvirB, 1.0, t2ABRS->get_pointer(), nvirA * nvirB);

    auto tBRAS = std::make_shared<Matrix>("tBRAS", aoccB * nvirA, aoccA * nvirB);
    double **tBRASp = tBRAS->pointer();

    for (int a = 0, ab = 0; a < aoccA; a++) {
        for (int b = 0; b < aoccB; b++, ab++) {
            for (int r = 0, rs = 0; r < nvirA; r++) {
                for (int s = 0; s < nvirB; s++, rs++) {
                    int br = b * nvirA + r;
                    int as = a * nvirB + s;
                    tBRASp[br][as] = tABRSp[ab][rs];
                }
            }
        }
    }

    auto t2BRAS = std::make_shared<Matrix>("t2BRAS", aoccB * nvirA, aoccA * nvirB);
    double **t2BRASp = t2BRAS->pointer();
    double **t2ABRSp = t2ABRS->pointer();

    for (int a = 0, ab = 0; a < aoccA; a++) {
        for (int b = 0; b < aoccB; b++, ab++) {
            for (int r = 0, rs = 0; r < nvirA; r++) {
                for (int s = 0; s < nvirB; s++, rs++) {
                    int br = b * nvirA + r;
                    int as = a * nvirB + s;
                    t2BRASp[br][as] = t2ABRSp[ab][rs];
                }
            }
        }
    }

    B_p_BB = get_DF_ints(BBintfile, BBlabel, foccB, noccB, foccB, noccB);
    double **B_p_RR = get_DF_ints(AAintfile, RRlabel, 0, nvirA, 0, nvirA);

    auto BRBR = std::make_shared<Matrix>("BRBR", aoccB * nvirA, aoccB * nvirA);
    double **BRBRp = BRBR->pointer();

    for (int b = 0, br = 0; b < aoccB; b++) {
        for (int r = 0; r < nvirA; r++, br++) {
            C_DGEMM('N', 'T', aoccB, nvirA, ndf_ + 3, 1.0, &(B_p_BB[b * aoccB][0]), ndf_ + 3, &(B_p_RR[r * nvirA][0]),
                    ndf_ + 3, 0.0, &(BRBRp[br][0]), nvirA);
        }
    }

    free_block(B_p_BB);
    free_block(B_p_RR);

    C_DGEMM('N', 'N', aoccB * nvirA, aoccA * nvirB, aoccB * nvirA, -1.0, BRBR->get_pointer(), aoccB * nvirA, tBRAS->get_pointer(),
            aoccA * nvirB, 1.0, t2BRAS->get_pointer(), aoccA * nvirB);

    B_p_AA = get_DF_ints(AAintfile, AAlabel, foccA, noccA, foccA, noccA);
    double **B_p_SS = get_DF_ints(BBintfile, SSlabel, 0, nvirB, 0, nvirB);

    auto ASAS = std::make_shared<Matrix>("ASAS", aoccA * nvirB, aoccA * nvirB);
    double **ASASp = ASAS->pointer();

    for (int a = 0, as = 0; a < aoccA; a++) {
        for (int s = 0; s < nvirB; s++, as++) {
            C_DGEMM('N', 'T', aoccA, nvirB, ndf_ + 3, 1.0, &(B_p_AA[a * aoccA][0]), ndf_ + 3, &(B_p_SS[s * nvirB][0]),
                    ndf_ + 3, 0.0, &(ASASp[as][0]), nvirB);
        }
    }

    free_block(B_p_AA);
    free_block(B_p_SS);

    C_DGEMM('N', 'N', aoccB * nvirA, aoccA * nvirB, aoccA * nvirB, -1.0, tBRAS->get_pointer(), aoccA * nvirB, ASAS->get_pointer(),
            aoccA * nvirB, 1.0, t2BRAS->get_pointer(), aoccA * nvirB);

    double energy = 4.0 * C_DDOT((long int)aoccA * aoccB * nvirA * nvirB, tBRAS->get_pointer(), 1, t2BRAS->get_pointer(), 1);

    return (energy);
}
}  // namespace sapt
}  // namespace psi
