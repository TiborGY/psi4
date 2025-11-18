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
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
namespace sapt {

void SAPT2::ind22() {
    double e_ind220 = ind220();

    if (debug_) {
        outfile->Printf("    Ind220              = %18.12lf [Eh]\n", e_ind220);
    }

    double e_ind202 = ind202();

    if (debug_) {
        outfile->Printf("    Ind202              = %18.12lf [Eh]\n\n", e_ind202);
    }

    e_ind22_ = e_ind220 + e_ind202;
    e_exch_ind22_ = e_ind22_ * (e_exch_ind20_ / e_ind20_);

    if (print_) {
        outfile->Printf("    Ind22               = %18.12lf [Eh]\n", e_ind22_);
    }
}

double SAPT2::ind220() {
    auto iAR = std::make_shared<Matrix>("iAR", aoccA_, nvirA_);
    double **iARp = iAR->pointer();

    for (int a = 0; a < aoccA_; a++) {
        for (int r = 0; r < nvirA_; r++) {
            iARp[a][r] = wBAR_[a + foccA_][r] / (evalsA_[a + foccA_] - evalsA_[r + noccA_]);
        }
    }

    auto iBS = std::make_shared<Matrix>("iBS", aoccB_, nvirB_);
    double **iBSp = iBS->pointer();

    for (int b = 0; b < aoccB_; b++) {
        for (int s = 0; s < nvirB_; s++) {
            iBSp[b][s] = wABS_[b + foccB_][s] / (evalsB_[b + foccB_] - evalsB_[s + noccB_]);
        }
    }

    double energy = 0.0;

    energy += ind220_1(PSIF_SAPT_AA_DF_INTS, "AA RI Integrals", "AR RI Integrals", "RR RI Integrals", PSIF_SAPT_AMPS,
                       "tARAR Amplitudes", iARp, wBAA_, wBRR_, foccA_, noccA_, nvirA_, evalsA_);

    energy += ind220_2(PSIF_SAPT_AMPS, "T2 AR Amplitudes", iARp, wBAA_, wBRR_, foccA_, noccA_, nvirA_);

    energy += ind220_3(PSIF_SAPT_AMPS, "pAA Density Matrix", "pRR Density Matrix", iARp, wBAR_, foccA_, noccA_, nvirA_);

    energy += ind220_4(PSIF_SAPT_AMPS, "Theta AR Intermediates", PSIF_SAPT_AA_DF_INTS, "AR RI Integrals", iARp, foccA_,
                       noccA_, nvirA_);

    energy += ind220_5(PSIF_SAPT_AMPS, "t2ARAR Amplitudes", iARp, foccA_, noccA_, nvirA_, evalsA_);

    energy += ind220_6(PSIF_SAPT_AA_DF_INTS, "AA RI Integrals", "AR RI Integrals", "RR RI Integrals", PSIF_SAPT_AMPS,
                       "tARAR Amplitudes", iARp, foccA_, noccA_, nvirA_);

    energy += ind220_7(PSIF_SAPT_AA_DF_INTS, "AA RI Integrals", "AR RI Integrals", "RR RI Integrals",
                       PSIF_SAPT_BB_DF_INTS, "BS RI Integrals", PSIF_SAPT_AMPS, "T2 AR Amplitudes",
                       "pAA Density Matrix", "pRR Density Matrix", iBSp, foccA_, noccA_, nvirA_, foccB_, noccB_, nvirB_);

    return (energy);
}

double SAPT2::ind202() {
    auto iAR = std::make_shared<Matrix>("iAR", aoccA_, nvirA_);
    double **iARp = iAR->pointer();

    for (int a = 0; a < aoccA_; a++) {
        for (int r = 0; r < nvirA_; r++) {
            iARp[a][r] = wBAR_[a + foccA_][r] / (evalsA_[a + foccA_] - evalsA_[r + noccA_]);
        }
    }

    auto iBS = std::make_shared<Matrix>("iBS", aoccB_, nvirB_);
    double **iBSp = iBS->pointer();

    for (int b = 0; b < aoccB_; b++) {
        for (int s = 0; s < nvirB_; s++) {
            iBSp[b][s] = wABS_[b + foccB_][s] / (evalsB_[b + foccB_] - evalsB_[s + noccB_]);
        }
    }

    double energy = 0.0;

    energy += ind220_1(PSIF_SAPT_BB_DF_INTS, "BB RI Integrals", "BS RI Integrals", "SS RI Integrals", PSIF_SAPT_AMPS,
                       "tBSBS Amplitudes", iBSp, wABB_, wASS_, foccB_, noccB_, nvirB_, evalsB_);

    energy += ind220_2(PSIF_SAPT_AMPS, "T2 BS Amplitudes", iBSp, wABB_, wASS_, foccB_, noccB_, nvirB_);

    energy += ind220_3(PSIF_SAPT_AMPS, "pBB Density Matrix", "pSS Density Matrix", iBSp, wABS_, foccB_, noccB_, nvirB_);

    energy += ind220_4(PSIF_SAPT_AMPS, "Theta BS Intermediates", PSIF_SAPT_BB_DF_INTS, "BS RI Integrals", iBSp, foccB_,
                       noccB_, nvirB_);

    energy += ind220_5(PSIF_SAPT_AMPS, "t2BSBS Amplitudes", iBSp, foccB_, noccB_, nvirB_, evalsB_);

    energy += ind220_6(PSIF_SAPT_BB_DF_INTS, "BB RI Integrals", "BS RI Integrals", "SS RI Integrals", PSIF_SAPT_AMPS,
                       "tBSBS Amplitudes", iBSp, foccB_, noccB_, nvirB_);

    energy += ind220_7(PSIF_SAPT_BB_DF_INTS, "BB RI Integrals", "BS RI Integrals", "SS RI Integrals",
                       PSIF_SAPT_AA_DF_INTS, "AR RI Integrals", PSIF_SAPT_AMPS, "T2 BS Amplitudes",
                       "pBB Density Matrix", "pSS Density Matrix", iARp, foccB_, noccB_, nvirB_, foccA_, noccA_, nvirA_);

    return (energy);
}

double SAPT2::ind220_1(int intfile, const char *AAlabel, const char *ARlabel, const char *RRlabel, int ampfile,
                       const char *tlabel, double **iAR, double **wBAA, double **wBRR, size_t foccA, size_t noccA, size_t nvirA,
                       double *evalsA) {
    int aoccA = noccA - foccA;

    auto C_p_AR = std::make_shared<Matrix>("C_p_AR", aoccA * nvirA, ndf_ + 3);
    double **B_p_RR = get_DF_ints(intfile, RRlabel, 0, nvirA, 0, nvirA);

    C_DGEMM('N', 'N', aoccA, nvirA * (ndf_ + 3), nvirA, 1.0, iAR[0], nvirA, B_p_RR[0], nvirA * (ndf_ + 3), 0.0,
            C_p_AR->get_pointer(), nvirA * (ndf_ + 3));

    free_block(B_p_RR);

    double **B_p_AA = get_DF_ints(intfile, AAlabel, foccA, noccA, foccA, noccA);
    double **C_p_ARp = C_p_AR->pointer();

    for (int a = 0; a < aoccA; a++) {
        C_DGEMM('T', 'N', nvirA, ndf_ + 3, aoccA, -1.0, iAR[0], nvirA, B_p_AA[a * aoccA], ndf_ + 3, 1.0,
                C_p_ARp[a * nvirA], ndf_ + 3);
    }

    free_block(B_p_AA);

    auto xARAR = std::make_shared<Matrix>("xARAR", aoccA * nvirA, aoccA * nvirA);
    double **xARARp = xARAR->pointer();
    double **B_p_AR = get_DF_ints(intfile, ARlabel, foccA, noccA, 0, nvirA);

    C_DGEMM('N', 'T', aoccA * nvirA, aoccA * nvirA, ndf_ + 3, 1.0, C_p_AR->get_pointer(), ndf_ + 3, B_p_AR[0], ndf_ + 3, 0.0,
            xARAR->get_pointer(), aoccA * nvirA);

    free_block(B_p_AR);

    auto tARAR = std::make_shared<Matrix>("tARAR", aoccA * nvirA, aoccA * nvirA);
    psio_->read_entry(ampfile, tlabel, (char *)tARAR->get_pointer(), sizeof(double) * aoccA * nvirA * aoccA * nvirA);

    C_DGEMM('N', 'N', aoccA, nvirA * aoccA * nvirA, aoccA, -1.0, &(wBAA[foccA][foccA]), noccA, tARAR->get_pointer(),
            nvirA * aoccA * nvirA, 1.0, xARAR->get_pointer(), nvirA * aoccA * nvirA);

    C_DGEMM('N', 'T', aoccA * nvirA * aoccA, nvirA, nvirA, 1.0, tARAR->get_pointer(), nvirA, wBRR[0], nvirA, 1.0, xARAR->get_pointer(), nvirA);

    symmetrize(xARAR->get_pointer(), aoccA, nvirA);

    auto yARAR = std::make_shared<Matrix>("yARAR", aoccA * nvirA, aoccA * nvirA);
    double **yARARp = yARAR->pointer();
    C_DCOPY((long int)aoccA * nvirA * aoccA * nvirA, xARAR->get_pointer(), 1, yARAR->get_pointer(), 1);
    antisym(yARARp, aoccA, nvirA);

    for (int a = 0, ar = 0; a < aoccA; a++) {
        for (int r = 0; r < nvirA; r++, ar++) {
            for (int aa = 0, aarr = 0; aa < aoccA; aa++) {
                for (int rr = 0; rr < nvirA; rr++, aarr++) {
                    xARARp[ar][aarr] /= evalsA[a + foccA] + evalsA[aa + foccA] - evalsA[r + noccA] - evalsA[rr + noccA];
                }
            }
        }
    }

    double energy = C_DDOT((long int)aoccA * nvirA * aoccA * nvirA, xARAR->get_pointer(), 1, yARAR->get_pointer(), 1);

    if (debug_) {
        outfile->Printf("\n    Ind22_1             = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}

double SAPT2::ind220_2(int ampfile, const char *tlabel, double **iAR, double **wBAA, double **wBRR, size_t foccA,
                       size_t noccA, size_t nvirA) {
    size_t aoccA = noccA - foccA;

    auto tAR = std::make_shared<Matrix>("tAR", aoccA, nvirA);
    psio_->read_entry(ampfile, tlabel, (char *)tAR->get_pointer(), sizeof(double) * aoccA * nvirA);

    auto zAR = std::make_shared<Matrix>("zAR", aoccA, nvirA);

    C_DGEMM('N', 'T', aoccA, nvirA, nvirA, 1.0, iAR[0], nvirA, wBRR[0], nvirA, 0.0, zAR->get_pointer(), nvirA);

    C_DGEMM('N', 'N', aoccA, nvirA, aoccA, -1.0, &(wBAA[foccA][foccA]), noccA, iAR[0], nvirA, 1.0, zAR->get_pointer(), nvirA);

    double energy = 4.0 * C_DDOT((long int)aoccA * nvirA, tAR->get_pointer(), 1, zAR->get_pointer(), 1);

    if (debug_) {
        outfile->Printf("    Ind22_2             = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}

double SAPT2::ind220_3(int ampfile, const char *AAlabel, const char *RRlabel, double **iAR, double **wBAR, size_t foccA,
                       size_t noccA, size_t nvirA) {
    size_t aoccA = noccA - foccA;

    auto pAA = std::make_shared<Matrix>("pAA", aoccA, aoccA);
    auto pRR = std::make_shared<Matrix>("pRR", nvirA, nvirA);

    psio_->read_entry(ampfile, AAlabel, (char *)pAA->get_pointer(), sizeof(double) * aoccA * aoccA);
    psio_->read_entry(ampfile, RRlabel, (char *)pRR->get_pointer(), sizeof(double) * nvirA * nvirA);

    auto xAA = std::make_shared<Matrix>("xAA", aoccA, aoccA);
    auto xRR = std::make_shared<Matrix>("xRR", nvirA, nvirA);

    C_DGEMM('N', 'T', aoccA, aoccA, nvirA, 1.0, iAR[0], nvirA, wBAR[foccA], nvirA, 0.0, xAA->get_pointer(), aoccA);
    C_DGEMM('T', 'N', nvirA, nvirA, aoccA, 1.0, iAR[0], nvirA, wBAR[foccA], nvirA, 0.0, xRR->get_pointer(), nvirA);

    double energy = 0.0;

    energy -= 2.0 * C_DDOT(aoccA * aoccA, pAA->get_pointer(), 1, xAA->get_pointer(), 1);
    energy -= 2.0 * C_DDOT(nvirA * nvirA, pRR->get_pointer(), 1, xRR->get_pointer(), 1);

    if (debug_) {
        outfile->Printf("    Ind22_3             = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}

double SAPT2::ind220_4(int ampfile, const char *thetalabel, int intfile, const char *ARlabel, double **iAR, size_t foccA,
                       size_t noccA, size_t nvirA) {
    size_t aoccA = noccA - foccA;

    auto xAA = std::make_shared<Matrix>("xAA", aoccA, aoccA);
    auto xRR = std::make_shared<Matrix>("xRR", nvirA, nvirA);

    C_DGEMM('N', 'T', aoccA, aoccA, nvirA, 1.0, iAR[0], nvirA, iAR[0], nvirA, 0.0, xAA->get_pointer(), aoccA);
    C_DGEMM('T', 'N', nvirA, nvirA, aoccA, 1.0, iAR[0], nvirA, iAR[0], nvirA, 0.0, xRR->get_pointer(), nvirA);

    double **B_p_AR = get_DF_ints(intfile, ARlabel, foccA, noccA, 0, nvirA);
    auto C_p_AR = std::make_shared<Matrix>("C_p_AR", aoccA * nvirA, ndf_ + 3);

    C_DGEMM('N', 'N', aoccA, nvirA * (ndf_ + 3), aoccA, 1.0, xAA->get_pointer(), aoccA, B_p_AR[0], nvirA * (ndf_ + 3), 0.0,
            C_p_AR->get_pointer(), nvirA * (ndf_ + 3));

    for (int a = 0; a < aoccA; a++) {
        C_DGEMM('N', 'N', nvirA, ndf_ + 3, nvirA, 1.0, xRR->get_pointer(), nvirA, B_p_AR[a * nvirA], ndf_ + 3, 1.0,
                C_p_AR->get_pointer(a * nvirA), ndf_ + 3);
    }

    free_block(B_p_AR);

    auto T_p_AR = std::make_shared<Matrix>("T_p_AR", aoccA * nvirA, ndf_ + 3);
    psio_->read_entry(ampfile, thetalabel, (char *)T_p_AR->get_pointer(), sizeof(double) * aoccA * nvirA * (ndf_ + 3));

    double energy = -2.0 * C_DDOT(aoccA * nvirA * (ndf_ + 3), C_p_AR->get_pointer(), 1, T_p_AR->get_pointer(), 1);

    if (debug_) {
        outfile->Printf("    Ind22_4             = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}

double SAPT2::ind220_5(int ampfile, const char *tlabel, double **iAR, size_t foccA, size_t noccA, size_t nvirA, double *evalsA) {
    size_t aoccA = noccA - foccA;

    auto tARAR = std::make_shared<Matrix>("tARAR", aoccA * nvirA, aoccA * nvirA);
    psio_->read_entry(ampfile, tlabel, (char *)tARAR->get_pointer(), sizeof(double) * aoccA * nvirA * aoccA * nvirA);
    double **tARARp = tARAR->pointer();
    antisym(tARARp, aoccA, nvirA);

    for (int a = 0, ar = 0; a < aoccA; a++) {
        for (int r = 0; r < nvirA; r++, ar++) {
            for (int aa = 0, aarr = 0; aa < aoccA; aa++) {
                for (int rr = 0; rr < nvirA; rr++, aarr++) {
                    tARARp[ar][aarr] *= evalsA[a + foccA] + evalsA[aa + foccA] - evalsA[r + noccA] - evalsA[rr + noccA];
                }
            }
        }
    }

    auto xAR = std::make_shared<Matrix>("xAR", aoccA, nvirA);

    C_DGEMV('n', aoccA * nvirA, aoccA * nvirA, 1.0, tARAR->get_pointer(), aoccA * nvirA, iAR[0], 1, 0.0, xAR->get_pointer(), 1);

    double energy = 2.0 * C_DDOT(aoccA * nvirA, xAR->get_pointer(), 1, iAR[0], 1);

    if (debug_) {
        outfile->Printf("    Ind22_5             = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}

double SAPT2::ind220_6(int intfile, const char *AAlabel, const char *ARlabel, const char *RRlabel, int ampfile,
                       const char *tlabel, double **iAR, size_t foccA, size_t noccA, size_t nvirA) {
    int aoccA = noccA - foccA;

    double **B_p_AR = get_DF_ints(intfile, ARlabel, foccA, noccA, 0, nvirA);
    auto gARAR = std::make_shared<Matrix>("gARAR", aoccA * nvirA, aoccA * nvirA);

    C_DGEMM('N', 'T', aoccA * nvirA, aoccA * nvirA, ndf_ + 3, 2.0, B_p_AR[0], ndf_ + 3, B_p_AR[0], ndf_ + 3, 0.0,
            gARAR->get_pointer(), aoccA * nvirA);

    free_block(B_p_AR);

    double **B_p_AA = get_DF_ints(intfile, AAlabel, foccA, noccA, foccA, noccA);
    double **B_p_RR = get_DF_ints(intfile, RRlabel, 0, nvirA, 0, nvirA);
    double **gARARp = gARAR->pointer();

    for (int a = 0, ar = 0; a < aoccA; a++) {
        for (int r = 0; r < nvirA; r++, ar++) {
            C_DGEMM('N', 'T', aoccA, nvirA, ndf_ + 3, -1.0, B_p_AA[a * aoccA], ndf_ + 3, B_p_RR[r * nvirA], ndf_ + 3,
                    1.0, gARARp[ar], nvirA);
        }
    }

    free_block(B_p_AA);
    free_block(B_p_RR);

    auto xAR = std::make_shared<Matrix>("xAR", aoccA, nvirA);
    auto yAR = std::make_shared<Matrix>("yAR", aoccA, nvirA);

    C_DGEMV('n', aoccA * nvirA, aoccA * nvirA, 1.0, gARAR->get_pointer(), aoccA * nvirA, iAR[0], 1, 0.0, xAR->get_pointer(), 1);

    auto tARAR = std::make_shared<Matrix>("tARAR", aoccA * nvirA, aoccA * nvirA);
    psio_->read_entry(ampfile, tlabel, (char *)tARAR->get_pointer(), sizeof(double) * aoccA * nvirA * aoccA * nvirA);
    double **tARARp = tARAR->pointer();
    antisym(tARARp, aoccA, nvirA);

    C_DGEMV('n', aoccA * nvirA, aoccA * nvirA, 1.0, tARAR->get_pointer(), aoccA * nvirA, iAR[0], 1, 0.0, yAR->get_pointer(), 1);

    double energy = -4.0 * C_DDOT(aoccA * nvirA, xAR->get_pointer(), 1, yAR->get_pointer(), 1);

    if (debug_) {
        outfile->Printf("    Ind22_6             = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}

double SAPT2::ind220_7(int AAfile, const char *AAlabel, const char *ARlabel, const char *RRlabel, int BBfile,
                       const char *BSlabel, int ampfile, const char *tlabel, const char *pAAlabel, const char *pRRlabel,
                       double **iBS, size_t foccA, size_t noccA, size_t nvirA, size_t foccB, size_t noccB, size_t nvirB) {
    size_t aoccA = noccA - foccA;
    size_t aoccB = noccB - foccB;

    auto pAA = std::make_shared<Matrix>("pAA", aoccA, aoccA);
    auto tAR = std::make_shared<Matrix>("tAR", aoccA, nvirA);
    auto pRR = std::make_shared<Matrix>("pRR", nvirA, nvirA);

    psio_->read_entry(ampfile, pAAlabel, (char *)pAA->get_pointer(), sizeof(double) * aoccA * aoccA);
    psio_->read_entry(ampfile, tlabel, (char *)tAR->get_pointer(), sizeof(double) * aoccA * nvirA);
    psio_->read_entry(ampfile, pRRlabel, (char *)pRR->get_pointer(), sizeof(double) * nvirA * nvirA);

    double *W = init_array(ndf_ + 3);
    double *X = init_array(ndf_ + 3);
    double *Y = init_array(ndf_ + 3);
    double *Z = init_array(ndf_ + 3);

    double **B_p_AA = get_DF_ints(AAfile, AAlabel, foccA, noccA, foccA, noccA);

    C_DGEMV('t', aoccA * aoccA, ndf_ + 3, 1.0, B_p_AA[0], ndf_ + 3, pAA->get_pointer(), 1, 0.0, W, 1);

    free_block(B_p_AA);

    double **B_p_RR = get_DF_ints(AAfile, RRlabel, 0, nvirA, 0, nvirA);

    C_DGEMV('t', nvirA * nvirA, ndf_ + 3, 1.0, B_p_RR[0], ndf_ + 3, pRR->get_pointer(), 1, 0.0, X, 1);

    free_block(B_p_RR);

    double **B_p_AR = get_DF_ints(AAfile, ARlabel, foccA, noccA, 0, nvirA);

    C_DGEMV('t', aoccA * nvirA, ndf_ + 3, 1.0, B_p_AR[0], ndf_ + 3, tAR->get_pointer(), 1, 0.0, Y, 1);

    free_block(B_p_AR);

    double **B_p_BS = get_DF_ints(BBfile, BSlabel, foccB, noccB, 0, nvirB);

    C_DGEMV('t', aoccB * nvirB, ndf_ + 3, 1.0, B_p_BS[0], ndf_ + 3, iBS[0], 1, 0.0, Z, 1);

    free_block(B_p_BS);

    double energy = 0.0;

    energy -= 8.0 * C_DDOT(ndf_ + 3, W, 1, Z, 1);
    energy += 8.0 * C_DDOT(ndf_ + 3, X, 1, Z, 1);
    energy += 16.0 * C_DDOT(ndf_ + 3, Y, 1, Z, 1);

    free(W);
    free(X);
    free(Y);
    free(Z);

    if (debug_) {
        outfile->Printf("    Ind22_7             = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}
}  // namespace sapt
}  // namespace psi
