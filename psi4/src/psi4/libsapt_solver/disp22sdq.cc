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

#include "sapt2p.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
namespace sapt {

void SAPT2p::disp22sdq() {
    double e_disp211 = disp211();

    if (debug_) {
        outfile->Printf("    Disp211             = %18.12lf [Eh]\n", e_disp211);
    }

    double e_disp220s = disp220s(PSIF_SAPT_AMPS, "T2 AR Amplitudes", "T AR Intermediates", PSIF_SAPT_AA_DF_INTS,
                                 "AA RI Integrals", "RR RI Integrals", foccA_, noccA_, nvirA_);

    if (debug_) {
        outfile->Printf("    Disp220 (S)         = %18.12lf [Eh]\n", e_disp220s);
    }

    double e_disp202s = disp220s(PSIF_SAPT_AMPS, "T2 BS Amplitudes", "T BS Intermediates", PSIF_SAPT_BB_DF_INTS,
                                 "BB RI Integrals", "SS RI Integrals", foccB_, noccB_, nvirB_);

    if (debug_) {
        outfile->Printf("    Disp202 (S)         = %18.12lf [Eh]\n", e_disp202s);
    }

    double e_disp220d = disp220d_1(PSIF_SAPT_AMPS, "t2ARAR Amplitudes", "T AR Intermediates", PSIF_SAPT_AA_DF_INTS,
                                   "AR RI Integrals", foccA_, noccA_, nvirA_);
    e_disp220d += disp220d_2(PSIF_SAPT_AMPS, "gARAR x tARBS", "Theta AR Intermediates", PSIF_SAPT_BB_DF_INTS,
                             "BS RI Integrals", foccA_, noccA_, nvirA_, foccB_, noccB_, nvirB_, evalsA_, evalsB_, 'N');

    if (debug_) {
        outfile->Printf("    Disp220 (D)         = %18.12lf [Eh]\n", e_disp220d);
    }

    double e_disp202d = disp220d_1(PSIF_SAPT_AMPS, "t2BSBS Amplitudes", "T BS Intermediates", PSIF_SAPT_BB_DF_INTS,
                                   "BS RI Integrals", foccB_, noccB_, nvirB_);
    e_disp202d += disp220d_2(PSIF_SAPT_AMPS, "gBSBS x tARBS", "Theta BS Intermediates", PSIF_SAPT_AA_DF_INTS,
                             "AR RI Integrals", foccB_, noccB_, nvirB_, foccA_, noccA_, nvirA_, evalsB_, evalsA_, 'T');

    if (debug_) {
        outfile->Printf("    Disp202 (D)         = %18.12lf [Eh]\n", e_disp202d);
    }

    double e_disp220q =
        disp220q_1(PSIF_SAPT_AMPS, "tARAR Amplitudes", "T AR Intermediates", "Theta AR Intermediates", aoccA_, nvirA_);
    e_disp220q += disp220q_2(PSIF_SAPT_AMPS, "pAA Density Matrix", "pRR Density Matrix", "T AR Intermediates",
                             PSIF_SAPT_AA_DF_INTS, "AR RI Integrals", foccA_, noccA_, nvirA_);
    e_disp220q += disp220q_3(PSIF_SAPT_AMPS, "tARAR Amplitudes", "tARBS Amplitudes", 'N', PSIF_SAPT_AA_DF_INTS,
                             "AR RI Integrals", foccA_, noccA_, nvirA_, foccB_, noccB_, nvirB_);
    e_disp220q += disp220q_4(PSIF_SAPT_AMPS, "tARAR Amplitudes", "tARBS Amplitudes", 'N', PSIF_SAPT_AA_DF_INTS,
                             "AR RI Integrals", foccA_, noccA_, nvirA_, foccB_, noccB_, nvirB_);

    if (debug_) {
        outfile->Printf("    Disp220 (Q)         = %18.12lf [Eh]\n", e_disp220q);
    }

    double e_disp202q =
        disp220q_1(PSIF_SAPT_AMPS, "tBSBS Amplitudes", "T BS Intermediates", "Theta BS Intermediates", aoccB_, nvirB_);
    e_disp202q += disp220q_2(PSIF_SAPT_AMPS, "pBB Density Matrix", "pSS Density Matrix", "T BS Intermediates",
                             PSIF_SAPT_BB_DF_INTS, "BS RI Integrals", foccB_, noccB_, nvirB_);
    e_disp202q += disp220q_3(PSIF_SAPT_AMPS, "tBSBS Amplitudes", "tARBS Amplitudes", 'T', PSIF_SAPT_BB_DF_INTS,
                             "BS RI Integrals", foccB_, noccB_, nvirB_, foccA_, noccA_, nvirA_);
    e_disp202q += disp220q_4(PSIF_SAPT_AMPS, "tBSBS Amplitudes", "tARBS Amplitudes", 'T', PSIF_SAPT_BB_DF_INTS,
                             "BS RI Integrals", foccB_, noccB_, nvirB_, foccA_, noccA_, nvirA_);

    if (debug_) {
        outfile->Printf("    Disp202 (Q)         = %18.12lf [Eh]\n\n", e_disp202q);
    }

    e_disp22sdq_ = e_disp211 + e_disp220s + e_disp202s + e_disp220d + e_disp202d + e_disp220q + e_disp202q;

    if (print_) {
        outfile->Printf("    Disp22 (SDQ)        = %18.12lf [Eh]\n", e_disp22sdq_);
    }
}

double SAPT2p::disp211() {
    double energy = 0.0;

    auto xARBS = std::make_shared<Matrix>("xARBS", aoccA_ * nvirA_, aoccB_ * nvirB_);
    double **xARBSp = xARBS->pointer();
    auto yARBS = std::make_shared<Matrix>("yARBS", aoccA_ * nvirA_, aoccB_ * nvirB_);

    psio_->read_entry(PSIF_SAPT_AMPS, "gBSBS x tARBS", (char *)xARBS->get_pointer(),
                      sizeof(double) * aoccA_ * nvirA_ * aoccB_ * nvirB_);

    psio_->read_entry(PSIF_SAPT_AMPS, "gARAR x tARBS", (char *)yARBS->get_pointer(),
                      sizeof(double) * aoccA_ * nvirA_ * aoccB_ * nvirB_);

    double **B_p_AR = get_DF_ints(PSIF_SAPT_AA_DF_INTS, "AR RI Integrals", foccA_, noccA_, 0, nvirA_);
    auto T_p_BS = std::make_shared<Matrix>("T_p_BS", aoccB_ * nvirB_, ndf_ + 3);
    psio_->read_entry(PSIF_SAPT_AMPS, "Theta BS Intermediates", (char *)T_p_BS->get_pointer(),
                      sizeof(double) * aoccB_ * nvirB_ * (ndf_ + 3));

    C_DGEMM('N', 'T', aoccA_ * nvirA_, aoccB_ * nvirB_, ndf_ + 3, 1.0, B_p_AR[0], ndf_ + 3, T_p_BS->get_pointer(), ndf_ + 3, 1.0,
            xARBS->get_pointer(), aoccB_ * nvirB_);

    free_block(B_p_AR);

    auto T_p_AR = std::make_shared<Matrix>("T_p_AR", aoccA_ * nvirA_, ndf_ + 3);
    psio_->read_entry(PSIF_SAPT_AMPS, "Theta AR Intermediates", (char *)T_p_AR->get_pointer(),
                      sizeof(double) * aoccA_ * nvirA_ * (ndf_ + 3));
    double **B_p_BS = get_DF_ints(PSIF_SAPT_BB_DF_INTS, "BS RI Integrals", foccB_, noccB_, 0, nvirB_);

    C_DGEMM('N', 'T', aoccA_ * nvirA_, aoccB_ * nvirB_, ndf_ + 3, 1.0, T_p_AR->get_pointer(), ndf_ + 3, B_p_BS[0], ndf_ + 3, 1.0,
            yARBS->get_pointer(), aoccB_ * nvirB_);

    free_block(B_p_BS);

    for (int a = 0, ar = 0; a < aoccA_; a++) {
        for (int r = 0; r < nvirA_; r++, ar++) {
            for (int b = 0, bs = 0; b < aoccB_; b++) {
                for (int s = 0; s < nvirB_; s++, bs++) {
                    xARBSp[ar][bs] /=
                        evalsA_[a + foccA_] + evalsB_[b + foccB_] - evalsA_[r + noccA_] - evalsB_[s + noccB_];
                }
            }
        }
    }

    energy = 8.0 * C_DDOT(aoccA_ * nvirA_ * aoccB_ * nvirB_, xARBS->get_pointer(), 1, yARBS->get_pointer(), 1);

    psio_->read_entry(PSIF_SAPT_AMPS, "tARBS Amplitudes", (char *)xARBS->get_pointer(),
                      sizeof(double) * aoccA_ * nvirA_ * aoccB_ * nvirB_);

    C_DGEMM('N', 'T', aoccA_ * nvirA_, aoccB_ * nvirB_, ndf_ + 3, 1.0, T_p_AR->get_pointer(), ndf_ + 3, T_p_BS->get_pointer(), ndf_ + 3, 0.0,
            yARBS->get_pointer(), aoccB_ * nvirB_);

    energy += 8.0 * C_DDOT(aoccA_ * nvirA_ * aoccB_ * nvirB_, xARBS->get_pointer(), 1, yARBS->get_pointer(), 1);

    return (energy);
}

double SAPT2p::disp220s(int ampfile, const char *tlabel, const char *thetalabel, int intfile, const char *AAlabel,
                        const char *RRlabel, size_t foccA, size_t noccA, size_t nvirA) {
    size_t aoccA = noccA - foccA;

    auto yAR = std::make_shared<Matrix>("yAR", aoccA, nvirA);

    auto T_p_AR = std::make_shared<Matrix>("T_p_AR", aoccA * nvirA, ndf_ + 3);
    psio_->read_entry(ampfile, thetalabel, (char *)T_p_AR->get_pointer(), sizeof(double) * aoccA * nvirA * (ndf_ + 3));

    double **B_p_RR = get_DF_ints(intfile, RRlabel, 0, nvirA, 0, nvirA);

    C_DGEMM('N', 'T', aoccA, nvirA, nvirA * (ndf_ + 3), 1.0, T_p_AR->get_pointer(), nvirA * (ndf_ + 3), B_p_RR[0],
            nvirA * (ndf_ + 3), 0.0, yAR->get_pointer(), nvirA);

    free_block(B_p_RR);

    double **B_p_AA = get_DF_ints(intfile, AAlabel, foccA, noccA, foccA, noccA);

    for (int a = 0; a < aoccA; a++) {
        C_DGEMM('N', 'T', aoccA, nvirA, ndf_ + 3, -1.0, B_p_AA[a * aoccA], ndf_ + 3, T_p_AR->get_pointer(a * nvirA), ndf_ + 3, 1.0,
                yAR->get_pointer(), nvirA);
    }

    free_block(B_p_AA);

    auto tAR = std::make_shared<Matrix>("tAR", aoccA, nvirA);
    psio_->read_entry(ampfile, tlabel, (char *)tAR->get_pointer(), sizeof(double) * aoccA * nvirA);

    double energy = 8.0 * C_DDOT(aoccA * nvirA, tAR->get_pointer(), 1, yAR->get_pointer(), 1);

    return (energy);
}

double SAPT2p::disp220d_1(int ampfile, const char *tlabel, const char *thetalabel, int intfile, const char *ARlabel,
                          size_t foccA, size_t noccA, size_t nvirA) {
    size_t aoccA = noccA - foccA;

    std::vector<double> xARAR((long int)aoccA * nvirA * aoccA * nvirA);
    double **B_p_AR = get_DF_ints(intfile, ARlabel, foccA, noccA, 0, nvirA);
    auto T_p_AR = std::make_shared<Matrix>("T_p_AR", aoccA * nvirA, ndf_ + 3);
    psio_->read_entry(ampfile, thetalabel, (char *)T_p_AR->get_pointer(), sizeof(double) * aoccA * nvirA * (ndf_ + 3));

    C_DGEMM('N', 'T', aoccA * nvirA, aoccA * nvirA, ndf_ + 3, 1.0, T_p_AR->get_pointer(), ndf_ + 3, B_p_AR[0], ndf_ + 3, 0.0, xARAR.data(),
            aoccA * nvirA);

    symmetrize(xARAR.data(), aoccA, nvirA);
    antisym(xARAR.data(), aoccA, nvirA);

    free_block(B_p_AR);

    std::vector<double> t2ARAR((long int)aoccA * nvirA * aoccA * nvirA);
    psio_->read_entry(ampfile, tlabel, (char *)t2ARAR.data(), sizeof(double) * aoccA * nvirA * aoccA * nvirA);

    double energy = 4.0 * C_DDOT((long int)aoccA * nvirA * aoccA * nvirA, xARAR.data(), 1, t2ARAR.data(), 1);

    if (debug_) {
        outfile->Printf("\n    Disp22d_1           = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}

double SAPT2p::disp220d_2(int ampfile, const char *glabel, const char *thetalabel, int intfile, const char *BSlabel,
                          size_t foccA, size_t noccA, size_t nvirA, size_t foccB, size_t noccB, size_t nvirB, double *evalsA,
                          double *evalsB, const char trans) {
    size_t aoccA = noccA - foccA;
    size_t aoccB = noccB - foccB;

    double energy = 0.0;

    auto T_p_AR = std::make_shared<Matrix>("T_p_AR", aoccA * nvirA, ndf_ + 3);
    psio_->read_entry(ampfile, thetalabel, (char *)T_p_AR->get_pointer(), sizeof(double) * aoccA * nvirA * (ndf_ + 3));
    double **B_p_BS = get_DF_ints(intfile, BSlabel, foccB, noccB, 0, nvirB);

    if (trans == 'n' || trans == 'N') {
        auto yARBS = std::make_shared<Matrix>("yARBS", aoccA * nvirA, aoccB * nvirB);
        double **yARBSp = yARBS->pointer();
        psio_->read_entry(ampfile, glabel, (char *)yARBS->get_pointer(), sizeof(double) * aoccA * nvirA * aoccB * nvirB);

        C_DGEMM('N', 'T', aoccA * nvirA, aoccB * nvirB, ndf_ + 3, 1.0, T_p_AR->get_pointer(), ndf_ + 3, B_p_BS[0], ndf_ + 3, 1.0,
                yARBS->get_pointer(), aoccB * nvirB);

        for (int a = 0, ar = 0; a < aoccA; a++) {
            for (int r = 0; r < nvirA; r++, ar++) {
                for (int b = 0, bs = 0; b < aoccB; b++) {
                    for (int s = 0; s < nvirB; s++, bs++) {
                        double tval = yARBSp[ar][bs];
                        tval *= tval;
                        energy += 4.0 * tval /
                                  (evalsA[a + foccA] + evalsB[b + foccB] - evalsA[r + noccA] - evalsB[s + noccB]);
                    }
                }
            }
        }
    } else if (trans == 't' || trans == 'T') {
        auto yBSAR = std::make_shared<Matrix>("yBSAR", aoccB * nvirB, aoccA * nvirA);
        double **yBSARp = yBSAR->pointer();
        psio_->read_entry(ampfile, glabel, (char *)yBSAR->get_pointer(), sizeof(double) * aoccA * nvirA * aoccB * nvirB);

        C_DGEMM('N', 'T', aoccB * nvirB, aoccA * nvirA, ndf_ + 3, 1.0, B_p_BS[0], ndf_ + 3, T_p_AR->get_pointer(), ndf_ + 3, 1.0,
                yBSAR->get_pointer(), aoccA * nvirA);

        for (int b = 0, bs = 0; b < aoccB; b++) {
            for (int s = 0; s < nvirB; s++, bs++) {
                for (int a = 0, ar = 0; a < aoccA; a++) {
                    for (int r = 0; r < nvirA; r++, ar++) {
                        double tval = yBSARp[bs][ar];
                        tval *= tval;
                        energy += 4.0 * tval /
                                  (evalsA[a + foccA] + evalsB[b + foccB] - evalsA[r + noccA] - evalsB[s + noccB]);
                    }
                }
            }
        }
    } else
        throw PsiException("You want me to do what to that matrix?", __FILE__, __LINE__);

    free_block(B_p_BS);

    if (debug_) {
        outfile->Printf("    Disp22d_2           = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}

double SAPT2p::disp220q_1(int ampfile, const char *tlabel, const char *Tlabel, const char *thetalabel, size_t aoccA,
                          size_t nvirA) {
    double energy = 0.0;

    auto thetaARAR = std::make_shared<Matrix>("thetaARAR", aoccA * nvirA, aoccA * nvirA);
    psio_->read_entry(ampfile, tlabel, (char *)thetaARAR->get_pointer(), sizeof(double) * aoccA * nvirA * aoccA * nvirA);
    double **thetaARARp = thetaARAR->pointer();
    antisym(thetaARARp, aoccA, nvirA);

    auto T_p_AR = std::make_shared<Matrix>("T_p_AR", aoccA * nvirA, ndf_ + 3);
    psio_->read_entry(ampfile, Tlabel, (char *)T_p_AR->get_pointer(), sizeof(double) * aoccA * nvirA * (ndf_ + 3));

    auto theta_p_AR = std::make_shared<Matrix>("theta_p_AR", aoccA * nvirA, ndf_ + 3);
    psio_->read_entry(ampfile, thetalabel, (char *)theta_p_AR->get_pointer(), sizeof(double) * aoccA * nvirA * (ndf_ + 3));

    auto xARAR = std::make_shared<Matrix>("xARAR", aoccA * nvirA, aoccA * nvirA);

    C_DGEMM('N', 'T', aoccA * nvirA, aoccA * nvirA, ndf_ + 3, 1.0, T_p_AR->get_pointer(), ndf_ + 3, theta_p_AR->get_pointer(), ndf_ + 3, 0.0,
            xARAR->get_pointer(), aoccA * nvirA);

    energy = 4.0 * C_DDOT((long int)aoccA * nvirA * aoccA * nvirA, xARAR->get_pointer(), 1, thetaARAR->get_pointer(), 1);

    if (debug_) {
        outfile->Printf("\n    Disp22q_1           = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}

double SAPT2p::disp220q_2(int ampfile, const char *pAAlabel, const char *pRRlabel, const char *Tlabel, int intfile,
                          const char *ARlabel, size_t foccA, size_t noccA, size_t nvirA) {
    size_t aoccA = noccA - foccA;

    auto pAA = std::make_shared<Matrix>("pAA", aoccA, aoccA);
    auto pRR = std::make_shared<Matrix>("pRR", nvirA, nvirA);

    psio_->read_entry(ampfile, pAAlabel, (char *)pAA->get_pointer(), sizeof(double) * aoccA * aoccA);
    psio_->read_entry(ampfile, pRRlabel, (char *)pRR->get_pointer(), sizeof(double) * nvirA * nvirA);

    auto qAA = std::make_shared<Matrix>("qAA", aoccA, aoccA);
    auto qRR = std::make_shared<Matrix>("qRR", nvirA, nvirA);

    double **B_p_AR = get_DF_ints(intfile, ARlabel, foccA, noccA, 0, nvirA);
    auto T_p_AR = std::make_shared<Matrix>("T_p_AR", aoccA * nvirA, ndf_ + 3);
    psio_->read_entry(ampfile, Tlabel, (char *)T_p_AR->get_pointer(), sizeof(double) * aoccA * nvirA * (ndf_ + 3));

    C_DGEMM('N', 'T', aoccA, aoccA, nvirA * (ndf_ + 3), 1.0, B_p_AR[0], nvirA * (ndf_ + 3), T_p_AR->get_pointer(),
            nvirA * (ndf_ + 3), 0.0, qAA->get_pointer(), aoccA);

    for (int a = 0; a < aoccA; a++) {
        C_DGEMM('N', 'T', nvirA, nvirA, ndf_ + 3, 1.0, B_p_AR[a * nvirA], ndf_ + 3, T_p_AR->get_pointer(a * nvirA), ndf_ + 3, 1.0,
                qRR->get_pointer(), nvirA);
    }

    free_block(B_p_AR);

    double energy = -4.0 * C_DDOT(aoccA * aoccA, pAA->get_pointer(), 1, qAA->get_pointer(), 1);
    energy -= 4.0 * C_DDOT(nvirA * nvirA, pRR->get_pointer(), 1, qRR->get_pointer(), 1);

    if (debug_) {
        outfile->Printf("    Disp22q_2           = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}

double SAPT2p::disp220q_3(int ampfile, const char *tARARlabel, const char *tARBSlabel, const char trans, int intfile,
                          const char *ARlabel, size_t foccA, size_t noccA, size_t nvirA, size_t foccB, size_t noccB, size_t nvirB) {
    size_t aoccA = noccA - foccA;
    size_t aoccB = noccB - foccB;

    auto xARAR = std::make_shared<Matrix>("xARAR", aoccA * nvirA, aoccA * nvirA);
    double **xARARp = xARAR->pointer();

    if (trans == 'n' || trans == 'N') {
        auto tARBS = std::make_shared<Matrix>("tARBS", aoccA * nvirA, aoccB * nvirB);
        psio_->read_entry(ampfile, tARBSlabel, (char *)tARBS->get_pointer(), sizeof(double) * aoccA * nvirA * aoccB * nvirB);

        C_DGEMM('N', 'T', aoccA * nvirA, aoccA * nvirA, aoccB * nvirB, 1.0, tARBS->get_pointer(), aoccB * nvirB, tARBS->get_pointer(),
                aoccB * nvirB, 0.0, xARAR->get_pointer(), aoccA * nvirA);
    } else if (trans == 't' || trans == 'T') {
        auto tBSAR = std::make_shared<Matrix>("tBSAR", aoccB * nvirB, aoccA * nvirA);
        psio_->read_entry(ampfile, tARBSlabel, (char *)tBSAR->get_pointer(), sizeof(double) * aoccA * nvirA * aoccB * nvirB);

        C_DGEMM('T', 'N', aoccA * nvirA, aoccA * nvirA, aoccB * nvirB, 1.0, tBSAR->get_pointer(), aoccA * nvirA, tBSAR->get_pointer(),
                aoccA * nvirA, 0.0, xARAR->get_pointer(), aoccA * nvirA);
    } else
        throw PsiException("You want me to do what to that matrix?", __FILE__, __LINE__);

    auto tARAR = std::make_shared<Matrix>("tARAR", aoccA * nvirA, aoccA * nvirA);
    psio_->read_entry(ampfile, tARARlabel, (char *)tARAR->get_pointer(), sizeof(double) * aoccA * nvirA * aoccA * nvirA);
    double **tARARp = tARAR->pointer();
    antisym(tARARp, aoccA, nvirA);

    auto yARAR = std::make_shared<Matrix>("yARAR", aoccA * nvirA, aoccA * nvirA);

    C_DGEMM('N', 'T', aoccA * nvirA, aoccA * nvirA, aoccA * nvirA, 1.0, xARAR->get_pointer(), aoccA * nvirA, tARAR->get_pointer(),
            aoccA * nvirA, 0.0, yARAR->get_pointer(), aoccA * nvirA);

    double **B_p_AR = get_DF_ints(intfile, ARlabel, foccA, noccA, 0, nvirA);
    C_DGEMM('N', 'T', aoccA * nvirA, aoccA * nvirA, ndf_ + 3, 1.0, B_p_AR[0], ndf_ + 3, B_p_AR[0], ndf_ + 3, 0.0,
            xARAR->get_pointer(), aoccA * nvirA);
    antisym(xARARp, aoccA, nvirA);

    double energy = 4.0 * C_DDOT((long int)aoccA * nvirA * aoccA * nvirA, xARAR->get_pointer(), 1, yARAR->get_pointer(), 1);

    free_block(B_p_AR);

    if (debug_) {
        outfile->Printf("    Disp22q_3           = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}

double SAPT2p::disp220q_4(int ampfile, const char *tARARlabel, const char *tARBSlabel, const char trans, int intfile,
                          const char *ARlabel, size_t foccA, size_t noccA, size_t nvirA, size_t foccB, size_t noccB, size_t nvirB) {
    size_t aoccA = noccA - foccA;
    int aoccB = noccB - foccB;

    auto rAA = std::make_shared<Matrix>("rAA", aoccA, aoccA);
    auto rRR = std::make_shared<Matrix>("rRR", nvirA, nvirA);

    auto tARAR = std::make_shared<Matrix>("tARAR", aoccA * nvirA, aoccA * nvirA);
    psio_->read_entry(ampfile, tARARlabel, (char *)tARAR->get_pointer(), sizeof(double) * aoccA * nvirA * aoccA * nvirA);

    auto gARAR = std::make_shared<Matrix>("gARAR", aoccA * nvirA, aoccA * nvirA);
    double **B_p_AR = get_DF_ints(intfile, ARlabel, foccA, noccA, 0, nvirA);
    C_DGEMM('N', 'T', aoccA * nvirA, aoccA * nvirA, ndf_ + 3, 1.0, B_p_AR[0], ndf_ + 3, B_p_AR[0], ndf_ + 3, 0.0,
            gARAR->get_pointer(), aoccA * nvirA);
    double **gARARp = gARAR->pointer();
    antisym(gARARp, aoccA, nvirA);

    C_DGEMM('N', 'T', aoccA, aoccA, nvirA * aoccA * nvirA, 1.0, tARAR->get_pointer(), nvirA * aoccA * nvirA, gARAR->get_pointer(),
            nvirA * aoccA * nvirA, 0.0, rAA->get_pointer(), aoccA);

    C_DGEMM('T', 'N', nvirA, nvirA, aoccA * nvirA * aoccA, 1.0, tARAR->get_pointer(), nvirA, gARAR->get_pointer(), nvirA, 0.0, rRR->get_pointer(), nvirA);

    free_block(B_p_AR);

    auto sAA = std::make_shared<Matrix>("sAA", aoccA, aoccA);
    auto sRR = std::make_shared<Matrix>("sRR", nvirA, nvirA);

    if (trans == 'n' || trans == 'N') {
        auto tARBS = std::make_shared<Matrix>("tARBS", aoccA * nvirA, aoccB * nvirB);
        psio_->read_entry(ampfile, tARBSlabel, (char *)tARBS->get_pointer(), sizeof(double) * aoccA * nvirA * aoccB * nvirB);

        C_DGEMM('N', 'T', aoccA, aoccA, nvirA * aoccB * nvirB, 1.0, tARBS->get_pointer(), nvirA * aoccB * nvirB, tARBS->get_pointer(),
                nvirA * aoccB * nvirB, 0.0, sAA->get_pointer(), aoccA);

        for (int a = 0; a < aoccA; a++) {
            C_DGEMM('N', 'T', nvirA, nvirA, aoccB * nvirB, 1.0, tARBS->get_pointer(a * nvirA), aoccB * nvirB, tARBS->get_pointer(a * nvirA),
                    aoccB * nvirB, 1.0, sRR->get_pointer(), nvirA);
        }
    } else if (trans == 't' || trans == 'T') {
        auto tBSAR = std::make_shared<Matrix>("tBSAR", aoccB * nvirB, aoccA * nvirA);
        double **tBSARp = tBSAR->pointer();
        psio_->read_entry(ampfile, tARBSlabel, (char *)tBSAR->get_pointer(), sizeof(double) * aoccA * nvirA * aoccB * nvirB);

        for (int b = 0, bs = 0; b < aoccB; b++) {
            for (int s = 0; s < nvirB; s++, bs++) {
                C_DGEMM('N', 'T', aoccA, aoccA, nvirA, 1.0, tBSARp[bs], nvirA, tBSARp[bs], nvirA, 1.0, sAA->get_pointer(), aoccA);
            }
        }

        C_DGEMM('T', 'N', nvirA, nvirA, aoccA * aoccB * nvirB, 1.0, tBSAR->get_pointer(), nvirA, tBSAR->get_pointer(), nvirA, 0.0, sRR->get_pointer(),
                nvirA);
    } else
        throw PsiException("You want me to do what to that matrix?", __FILE__, __LINE__);

    double energy = -4.0 * C_DDOT(aoccA * aoccA, rAA->get_pointer(), 1, sAA->get_pointer(), 1);
    energy -= 4.0 * C_DDOT(nvirA * nvirA, rRR->get_pointer(), 1, sRR->get_pointer(), 1);

    if (debug_) {
        outfile->Printf("    Disp22q_4           = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}
}  // namespace sapt
}  // namespace psi
