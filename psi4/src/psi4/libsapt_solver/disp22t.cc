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

#include <ctime>

#include "sapt2p.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
namespace sapt {

void SAPT2p::disp22t() {
    if (print_) {
        outfile->Printf("\n");
    }

    double e_disp220t;

    if (nat_orbs_t3_) {
        e_disp220t = disp220t(PSIF_SAPT_AA_DF_INTS, "AA RI Integrals", "AR NO RI Integrals", "RR NO RI Integrals",
                              PSIF_SAPT_BB_DF_INTS, "BS NO RI Integrals", PSIF_SAPT_AMPS, "tARAR NO Amplitudes", foccA_,
                              noccA_, no_nvirA_, foccB_, noccB_, no_nvirB_, no_evalsA_, no_evalsB_);
    } else {
        e_disp220t = disp220t(PSIF_SAPT_AA_DF_INTS, "AA RI Integrals", "AR RI Integrals", "RR RI Integrals",
                              PSIF_SAPT_BB_DF_INTS, "BS RI Integrals", PSIF_SAPT_AMPS, "tARAR Amplitudes", foccA_,
                              noccA_, nvirA_, foccB_, noccB_, nvirB_, evalsA_, evalsB_);
    }

    if (print_) {
        outfile->Printf("\n    Disp220 (T)         = %18.12lf [Eh]\n\n", e_disp220t);
    }

    double e_disp202t;

    if (nat_orbs_t3_) {
        e_disp202t = disp220t(PSIF_SAPT_BB_DF_INTS, "BB RI Integrals", "BS NO RI Integrals", "SS NO RI Integrals",
                              PSIF_SAPT_AA_DF_INTS, "AR NO RI Integrals", PSIF_SAPT_AMPS, "tBSBS NO Amplitudes", foccB_,
                              noccB_, no_nvirB_, foccA_, noccA_, no_nvirA_, no_evalsB_, no_evalsA_);
    } else {
        e_disp202t = disp220t(PSIF_SAPT_BB_DF_INTS, "BB RI Integrals", "BS RI Integrals", "SS RI Integrals",
                              PSIF_SAPT_AA_DF_INTS, "AR RI Integrals", PSIF_SAPT_AMPS, "tBSBS Amplitudes", foccB_,
                              noccB_, nvirB_, foccA_, noccA_, nvirA_, evalsB_, evalsA_);
    }

    if (print_) {
        outfile->Printf("\n    Disp202 (T)         = %18.12lf [Eh]\n\n", e_disp202t);
    }

    e_disp22t_ = e_disp220t + e_disp202t;

    if (print_) {
        outfile->Printf("    Disp22 (T)          = %18.12lf [Eh]\n", e_disp22t_);
    }

    if (nat_orbs_t3_) {
        double scale = e_disp20_ / e_no_disp20_;
        e_disp220t *= scale;
        e_disp202t *= scale;
        e_est_disp22t_ = e_disp220t + e_disp202t;

        if (print_) {
            outfile->Printf("\n    Est. Disp220 (T)    = %18.12lf [Eh]\n", e_disp220t);
            outfile->Printf("    Est. Disp202 (T)    = %18.12lf [Eh]\n\n", e_disp202t);
            outfile->Printf("    Est. Disp22 (T)     = %18.12lf [Eh]\n", e_est_disp22t_);
        }
    }
}

void SAPT2p::disp22tccd() {
    if (print_) {
        outfile->Printf("\n");
    }

    if (nat_orbs_t3_) {
        natural_orbitalify_ccd();
    }

    double e_disp220t;

    if (nat_orbs_t3_) {
        e_disp220t = disp220tccd(PSIF_SAPT_AA_DF_INTS, "AA RI Integrals", PSIF_SAPT_AA_DF_INTS, "AR NO RI Integrals",
                                 "RR NO RI Integrals", PSIF_SAPT_BB_DF_INTS, "BS NO RI Integrals", PSIF_SAPT_CCD,
                                 "T ARAR Natorb Amplitudes", "T BSAR Natorb Amplitudes", no_evalsA_, no_evalsB_, noccA_,
                                 no_nvirA_, foccA_, noccB_, no_nvirB_, foccB_);
    } else {
        e_disp220t =
            disp220tccd(PSIF_SAPT_AA_DF_INTS, "AA RI Integrals", PSIF_SAPT_AA_DF_INTS, "AR RI Integrals",
                        "RR RI Integrals", PSIF_SAPT_BB_DF_INTS, "BS RI Integrals", PSIF_SAPT_CCD, "T ARAR Amplitudes",
                        "T BSAR Amplitudes", evalsA_, evalsB_, noccA_, nvirA_, foccA_, noccB_, nvirB_, foccB_);
    }

    if (print_) {
        outfile->Printf("\n    Disp220 (T)         = %18.12lf [Eh]\n\n", e_disp220t);
    }

    double e_disp202t;

    if (nat_orbs_t3_) {
        e_disp202t = disp220tccd(PSIF_SAPT_BB_DF_INTS, "BB RI Integrals", PSIF_SAPT_BB_DF_INTS, "BS NO RI Integrals",
                                 "SS NO RI Integrals", PSIF_SAPT_AA_DF_INTS, "AR NO RI Integrals", PSIF_SAPT_CCD,
                                 "T BSBS Natorb Amplitudes", "T ARBS Natorb Amplitudes", no_evalsB_, no_evalsA_, noccB_,
                                 no_nvirB_, foccB_, noccA_, no_nvirA_, foccA_);
    } else {
        e_disp202t =
            disp220tccd(PSIF_SAPT_BB_DF_INTS, "BB RI Integrals", PSIF_SAPT_BB_DF_INTS, "BS RI Integrals",
                        "SS RI Integrals", PSIF_SAPT_AA_DF_INTS, "AR RI Integrals", PSIF_SAPT_CCD, "T BSBS Amplitudes",
                        "T ARBS Amplitudes", evalsB_, evalsA_, noccB_, nvirB_, foccB_, noccA_, nvirA_, foccA_);
    }

    if (print_) {
        outfile->Printf("\n    Disp202 (T)         = %18.12lf [Eh]\n\n", e_disp202t);
    }

    e_disp22t_ccd_ = e_disp220t + e_disp202t;

    if (print_) {
        outfile->Printf("    Disp22 (T)          = %18.12lf [Eh]\n", e_disp22t_ccd_);
    }

    if (nat_orbs_t3_) {
        double scale = e_disp20_ / e_no_disp20_;
        e_disp220t *= scale;
        e_disp202t *= scale;
        e_est_disp22t_ccd_ = e_disp220t + e_disp202t;

        if (print_) {
            outfile->Printf("\n    Est. Disp220 (T)    = %18.12lf [Eh]\n", e_disp220t);
            outfile->Printf("    Est. Disp202 (T)    = %18.12lf [Eh]\n\n", e_disp202t);
            outfile->Printf("    Est. Disp22 (T)     = %18.12lf [Eh]\n", e_est_disp22t_ccd_);
        }
    }
}

void SAPT2p::natural_orbitalify_ccd() {
    int occA = noccA_ - foccA_;
    int occB = noccB_ - foccB_;

    auto tARAR = std::make_shared<Matrix>("tARAR", occA * nvirA_, occA * nvirA_);

    psio_->read_entry(PSIF_SAPT_CCD, "T ARAR Amplitudes", (char *)tARAR->get_pointer(),
                      occA * nvirA_ * occA * nvirA_ * (size_t)sizeof(double));

    auto tARAr = std::make_shared<Matrix>("tARAr", occA * nvirA_, occA * no_nvirA_);

    C_DGEMM('N', 'N', occA * nvirA_ * occA, no_nvirA_, nvirA_, 1.0, tARAR->get_pointer(), nvirA_, no_CA_[0], no_nvirA_, 0.0,
            tARAr->get_pointer(), no_nvirA_);

    auto tArAr = std::make_shared<Matrix>("tArAr", occA * no_nvirA_, occA * no_nvirA_);
    double **tARArp = tARAr->pointer();
    double **tArArp = tArAr->pointer();

    for (int a = 0; a < occA; a++) {
        C_DGEMM('T', 'N', no_nvirA_, occA * no_nvirA_, nvirA_, 1.0, no_CA_[0], no_nvirA_, tARArp[a * nvirA_],
                occA * no_nvirA_, 0.0, tArArp[a * no_nvirA_], occA * no_nvirA_);
    }

    psio_->write_entry(PSIF_SAPT_CCD, "T ARAR Natorb Amplitudes", (char *)tArAr->get_pointer(),
                       occA * no_nvirA_ * occA * no_nvirA_ * (size_t)sizeof(double));

    auto tBSBS = std::make_shared<Matrix>("tBSBS", occB * nvirB_, occB * nvirB_);

    psio_->read_entry(PSIF_SAPT_CCD, "T BSBS Amplitudes", (char *)tBSBS->get_pointer(),
                      occB * nvirB_ * occB * nvirB_ * (size_t)sizeof(double));

    auto tBSBs = std::make_shared<Matrix>("tBSBs", occB * nvirB_, occB * no_nvirB_);

    C_DGEMM('N', 'N', occB * nvirB_ * occB, no_nvirB_, nvirB_, 1.0, tBSBS->get_pointer(), nvirB_, no_CB_[0], no_nvirB_, 0.0,
            tBSBs->get_pointer(), no_nvirB_);

    auto tBsBs = std::make_shared<Matrix>("tBsBs", occB * no_nvirB_, occB * no_nvirB_);
    double **tBSBsp = tBSBs->pointer();
    double **tBsBsp = tBsBs->pointer();

    for (int b = 0; b < occB; b++) {
        C_DGEMM('T', 'N', no_nvirB_, occB * no_nvirB_, nvirB_, 1.0, no_CB_[0], no_nvirB_, tBSBsp[b * nvirB_],
                occB * no_nvirB_, 0.0, tBsBsp[b * no_nvirB_], occB * no_nvirB_);
    }

    psio_->write_entry(PSIF_SAPT_CCD, "T BSBS Natorb Amplitudes", (char *)tBsBs->get_pointer(),
                       occB * no_nvirB_ * occB * no_nvirB_ * (size_t)sizeof(double));

    auto tARBS = std::make_shared<Matrix>("tARBS", occA * nvirA_, occB * nvirB_);

    psio_->read_entry(PSIF_SAPT_CCD, "T ARBS Amplitudes", (char *)tARBS->get_pointer(),
                      occA * nvirA_ * occB * nvirB_ * (size_t)sizeof(double));

    auto tARBs = std::make_shared<Matrix>("tARBs", occA * nvirA_, occB * no_nvirB_);

    C_DGEMM('N', 'N', occA * nvirA_ * occB, no_nvirB_, nvirB_, 1.0, tARBS->get_pointer(), nvirB_, no_CB_[0], no_nvirB_, 0.0,
            tARBs->get_pointer(), no_nvirB_);

    auto tArBs = std::make_shared<Matrix>("tArBs", occA * no_nvirA_, occB * no_nvirB_);
    double **tARBsp = tARBs->pointer();
    double **tArBsp = tArBs->pointer();

    for (int a = 0; a < occA; a++) {
        C_DGEMM('T', 'N', no_nvirA_, occB * no_nvirB_, nvirA_, 1.0, no_CA_[0], no_nvirA_, tARBsp[a * nvirA_],
                occB * no_nvirB_, 0.0, tArBsp[a * no_nvirA_], occB * no_nvirB_);
    }

    auto tBsAr = std::make_shared<Matrix>("tBsAr", occB * no_nvirB_, occA * no_nvirA_);
    double **tBsArp = tBsAr->pointer();

    for (int a1 = 0, a1r1 = 0; a1 < occA; a1++) {
        for (int r1 = 0; r1 < no_nvirA_; r1++, a1r1++) {
            for (int b1 = 0, b1s1 = 0; b1 < occB; b1++) {
                for (int s1 = 0; s1 < no_nvirB_; s1++, b1s1++) {
                    tBsArp[b1s1][a1r1] = tArBsp[a1r1][b1s1];
                }
            }
        }
    }

    psio_->write_entry(PSIF_SAPT_CCD, "T ARBS Natorb Amplitudes", (char *)tArBs->get_pointer(),
                       occA * no_nvirA_ * occB * no_nvirB_ * (size_t)sizeof(double));
    psio_->write_entry(PSIF_SAPT_CCD, "T BSAR Natorb Amplitudes", (char *)tBsAr->get_pointer(),
                       occA * no_nvirA_ * occB * no_nvirB_ * (size_t)sizeof(double));
}

double SAPT2p::disp220t(int AAfile, const char *AAlabel, const char *ARlabel, const char *RRlabel, int BBfile,
                        const char *BSlabel, int ampfile, const char *tlabel, size_t foccA, size_t noccA, size_t nvirA,
                        size_t foccB, size_t noccB, size_t nvirB, double *evalsA, double *evalsB) {
    double energy = 0.0;

    size_t aoccA = noccA - foccA;
    size_t aoccB = noccB - foccB;

    auto wARAR = std::make_shared<Matrix>("wARAR", aoccA * nvirA, aoccA * nvirA);
    double **wARARp = wARAR->pointer();

    auto vbsAA = std::make_shared<Matrix>("vbsAA", aoccA, aoccA);
    auto vbsRR = std::make_shared<Matrix>("vbsRR", nvirA, nvirA);
    auto vARAA = std::make_shared<Matrix>("vARAA", aoccA * nvirA, aoccA * aoccA);

    auto tARAR = std::make_shared<Matrix>("tARAR", aoccA * nvirA, aoccA * nvirA);
    psio_->read_entry(ampfile, tlabel, (char *)tARAR->get_pointer(), sizeof(double) * aoccA * nvirA * aoccA * nvirA);
    auto tbsAR = std::make_shared<Matrix>("tbsAR", aoccA, nvirA);
    double **tbsARp = tbsAR->pointer();

    double **B_p_AA = get_DF_ints(AAfile, AAlabel, foccA, noccA, foccA, noccA);
    double **B_p_AR = get_DF_ints(AAfile, ARlabel, foccA, noccA, 0, nvirA);
    double **B_p_RR = get_DF_ints(AAfile, RRlabel, 0, nvirA, 0, nvirA);
    double *B_p_bs = init_array(ndf_ + 3);

    auto C_p_AR = std::make_shared<Matrix>("C_p_AR", aoccA * nvirA, ndf_ + 3);

    C_DGEMM('N', 'T', aoccA * nvirA, aoccA * aoccA, ndf_ + 3, 1.0, &(B_p_AR[0][0]), ndf_ + 3, &(B_p_AA[0][0]), ndf_ + 3,
            0.0, vARAA->get_pointer(), aoccA * aoccA);

    std::time_t start = std::time(nullptr);
    std::time_t stop;

    for (size_t b = 0, bs = 0; b < aoccB; b++) {
        for (int s = 0; s < nvirB; s++, bs++) {
            psio_address next_DF_BS = psio_get_address(
                PSIO_ZERO, sizeof(double) * (b + foccB) * nvirB * (ndf_ + 3) + sizeof(double) * s * (ndf_ + 3));
            psio_->read(BBfile, BSlabel, (char *)&(B_p_bs[0]), sizeof(double) * (ndf_ + 3), next_DF_BS, &next_DF_BS);

            C_DGEMV('n', aoccA * nvirA, ndf_ + 3, 1.0, B_p_AR[0], ndf_ + 3, B_p_bs, 1, 0.0, tbsAR->get_pointer(), 1);

            for (int a = 0, ar = 0; a < aoccA; a++) {
                for (int r = 0; r < nvirA; r++, ar++) {
                    double denom = evalsA[a + foccA] + evalsB[b + foccB] - evalsA[r + noccA] - evalsB[s + noccB];
                    tbsARp[a][r] /= denom;
                }
            }

            C_DGEMV('n', aoccA * aoccA, ndf_ + 3, 1.0, B_p_AA[0], ndf_ + 3, B_p_bs, 1, 0.0, vbsAA->get_pointer(), 1);
            C_DGEMV('n', nvirA * nvirA, ndf_ + 3, 1.0, B_p_RR[0], ndf_ + 3, B_p_bs, 1, 0.0, vbsRR->get_pointer(), 1);

            C_DGEMM('N', 'N', aoccA * nvirA * aoccA, nvirA, nvirA, 1.0, tARAR->get_pointer(), nvirA, vbsRR->get_pointer(), nvirA,
                    0.0, wARAR->get_pointer(), nvirA);
            C_DGEMM('N', 'N', aoccA, nvirA * aoccA * nvirA, aoccA, -1.0, vbsAA->get_pointer(), aoccA, tARAR->get_pointer(),
                    nvirA * aoccA * nvirA, 1.0, wARAR->get_pointer(), nvirA * aoccA * nvirA);
            C_DGEMM('N', 'N', aoccA * nvirA * aoccA, nvirA, aoccA, -1.0, vARAA->get_pointer(), aoccA, tbsAR->get_pointer(), nvirA,
                    1.0, wARAR->get_pointer(), nvirA);
            C_DGEMM('N', 'N', aoccA, nvirA * (ndf_ + 3), nvirA, 1.0, tbsAR->get_pointer(), nvirA, &(B_p_RR[0][0]),
                    nvirA * (ndf_ + 3), 0.0, C_p_AR->get_pointer(), nvirA * (ndf_ + 3));
            C_DGEMM('N', 'T', aoccA * nvirA, aoccA * nvirA, ndf_ + 3, 1.0, &(B_p_AR[0][0]), ndf_ + 3, C_p_AR->get_pointer(),
                    ndf_ + 3, 1.0, wARAR->get_pointer(), aoccA * nvirA);

            for (int a = 0, ar = 0; a < aoccA; a++) {
                for (int r = 0; r < nvirA; r++, ar++) {
                    for (int a1 = 0, a1r1 = 0; a1 < aoccA; a1++) {
                        for (int r1 = 0; r1 < nvirA; r1++, a1r1++) {
                            int a1r = a1 * nvirA + r;
                            int ar1 = a * nvirA + r1;
                            double tval1 = wARARp[ar][a1r1] + wARARp[a1r1][ar];
                            double tval2 = wARARp[a1r][ar1] + wARARp[ar1][a1r];
                            double denom = evalsA[a + foccA] + evalsA[a1 + foccA] + evalsB[b + foccB] -
                                           evalsA[r + noccA] - evalsA[r1 + noccA] - evalsB[s + noccB];
                            energy += ((4.0 * tval1 - 2.0 * tval2) * tval1) / denom;
                        }
                    }
                }
            }
        }
        stop = std::time(nullptr);
        if (print_) {
            outfile->Printf("    (i = %3zu of %3zu) %10ld seconds\n", b + 1, aoccB, stop - start);
        }
    }

    free(B_p_bs);
    free_block(B_p_AA);
    free_block(B_p_AR);
    free_block(B_p_RR);

    return (energy);
}

double SAPT2p::disp220tccd(int AAnum, const char *AA_label, int Rnum, const char *AR_label, const char *RR_label,
                           int BBnum, const char *BS_label, int ampnum, const char *tarar, const char *tbsar,
                           double *evalsA, double *evalsB, size_t noccA, size_t nvirA, size_t foccA, size_t noccB, size_t nvirB,
                           size_t foccB) {
    double energy = 0.0;

    noccA -= foccA;
    noccB -= foccB;

    auto w_ARAR = std::make_shared<Matrix>("w_ARAR", noccA * nvirA, noccA * nvirA);
    double **w_ARARp = w_ARAR->pointer();

    auto v_bsAA = std::make_shared<Matrix>("v_bsAA", noccA, noccA);
    auto v_bsRR = std::make_shared<Matrix>("v_bsRR", nvirA, nvirA);
    auto v_ARAA = std::make_shared<Matrix>("v_ARAA", noccA * nvirA, noccA * noccA);

    double **B_p_AA = get_DF_ints_nongimp(AAnum, AA_label, foccA, noccA + foccA, foccA, noccA + foccA);
    double **B_p_AR = get_DF_ints_nongimp(Rnum, AR_label, foccA, noccA + foccA, 0, nvirA);
    double **B_p_RR = get_DF_ints_nongimp(Rnum, RR_label, 0, nvirA, 0, nvirA);
    double *B_p_bs = init_array(ndf_);

    auto t_bsAR = std::make_shared<Matrix>("t_bsAR", noccA, nvirA);
    double **t_bsARp = t_bsAR->pointer();
    std::shared_ptr<Matrix> t_ARAR;
    double **t_ARARp;

    psio_address next_ARAR;

    if (ampnum == PSIF_SAPT_CCD) {
        t_ARAR = std::make_shared<Matrix>("t_ARAR", noccA * nvirA, noccA * nvirA);
        psio_->read_entry(ampnum, tarar, (char *)t_ARAR->get_pointer(), noccA * nvirA * noccA * nvirA * (size_t)sizeof(double));
        t_ARARp = t_ARAR->pointer();
    } else if (ampnum) {
        t_ARAR = std::make_shared<Matrix>("t_ARAR", noccA * nvirA, noccA * nvirA);
        t_ARARp = t_ARAR->pointer();
        for (int a = 0, ar = 0; a < noccA; a++) {
            for (int r = 0; r < nvirA; r++, ar++) {
                next_ARAR = psio_get_address(
                    PSIO_ZERO, ((foccA * nvirA + ar) * (noccA + foccA) * nvirA + foccA * nvirA) * sizeof(double));
                psio_->read(ampnum, tarar, (char *)t_ARARp[ar], noccA * nvirA * (size_t)sizeof(double), next_ARAR,
                            &next_ARAR);
            }
        }
    } else {
        t_ARAR = std::make_shared<Matrix>("t_ARAR", noccA * nvirA, noccA * nvirA);
        t_ARARp = t_ARAR->pointer();
        C_DGEMM('N', 'T', noccA * nvirA, noccA * nvirA, ndf_, 1.0, &(B_p_AR[0][0]), ndf_, &(B_p_AR[0][0]), ndf_, 0.0,
                t_ARAR->get_pointer(), noccA * nvirA);

        for (int a = 0, ar = 0; a < noccA; a++) {
            for (int r = 0; r < nvirA; r++, ar++) {
                for (int a1 = 0, a1r1 = 0; a1 < noccA; a1++) {
                    for (int r1 = 0; r1 < nvirA; r1++, a1r1++) {
                        double denom = evalsA[a + foccA] + evalsA[a1 + foccA] - evalsA[r + noccA + foccA] -
                                       evalsA[r1 + noccA + foccA];
                        t_ARARp[ar][a1r1] /= denom;
                    }
                }
            }
        }
    }

    auto C_p_AR = std::make_shared<Matrix>("C_p_AR", noccA * nvirA, ndf_);

    C_DGEMM('N', 'T', noccA * nvirA, noccA * noccA, ndf_, 1.0, &(B_p_AR[0][0]), ndf_, &(B_p_AA[0][0]), ndf_, 0.0,
            v_ARAA->get_pointer(), noccA * noccA);

    psio_address next_BSAR;

    std::time_t start = std::time(nullptr);
    std::time_t stop;

    for (size_t b = 0, bs = 0; b < noccB; b++) {
        for (int s = 0; s < nvirB; s++, bs++) {
            psio_address next_DF_BS =
                psio_get_address(PSIO_ZERO, ((foccB + b) * nvirB + s) * (ndf_ + 3) * (size_t)sizeof(double));
            psio_->read(BBnum, BS_label, (char *)&(B_p_bs[0]), sizeof(double) * ndf_, next_DF_BS, &next_DF_BS);

            if (ampnum == PSIF_SAPT_CCD) {
                next_BSAR = psio_get_address(PSIO_ZERO, bs * noccA * nvirA * sizeof(double));
                psio_->read(ampnum, tbsar, (char *)t_bsAR->get_pointer(), sizeof(double) * noccA * nvirA, next_BSAR, &next_BSAR);
            } else if (ampnum) {
                next_BSAR = psio_get_address(
                    PSIO_ZERO, ((foccB * nvirB + bs) * (noccA + foccA) * nvirA + foccA * nvirA) * sizeof(double));
                psio_->read(ampnum, tbsar, (char *)t_bsAR->get_pointer(), sizeof(double) * noccA * nvirA, next_BSAR, &next_BSAR);
            } else {
                C_DGEMV('n', noccA * nvirA, ndf_, 1.0, B_p_AR[0], ndf_, B_p_bs, 1, 0.0, t_bsAR->get_pointer(), 1);

                for (int a = 0; a < noccA; a++) {
                    for (int r = 0; r < nvirA; r++) {
                        double denom = evalsA[a + foccA] + evalsB[b + foccB] - evalsA[r + noccA + foccA] -
                                       evalsB[s + noccB + foccB];
                        t_bsARp[a][r] /= denom;
                    }
                }
            }

            C_DGEMV('n', noccA * noccA, ndf_, 1.0, B_p_AA[0], ndf_, B_p_bs, 1, 0.0, v_bsAA->get_pointer(), 1);
            C_DGEMV('n', nvirA * nvirA, ndf_, 1.0, B_p_RR[0], ndf_, B_p_bs, 1, 0.0, v_bsRR->get_pointer(), 1);

            C_DGEMM('N', 'N', noccA * nvirA * noccA, nvirA, nvirA, 1.0, t_ARAR->get_pointer(), nvirA, v_bsRR->get_pointer(), nvirA,
                    0.0, w_ARAR->get_pointer(), nvirA);
            C_DGEMM('N', 'N', noccA, nvirA * noccA * nvirA, noccA, -1.0, v_bsAA->get_pointer(), noccA, t_ARAR->get_pointer(),
                    nvirA * noccA * nvirA, 1.0, w_ARAR->get_pointer(), nvirA * noccA * nvirA);
            C_DGEMM('N', 'N', noccA * nvirA * noccA, nvirA, noccA, -1.0, v_ARAA->get_pointer(), noccA, t_bsAR->get_pointer(), nvirA,
                    1.0, w_ARAR->get_pointer(), nvirA);
            C_DGEMM('N', 'N', noccA, nvirA * ndf_, nvirA, 1.0, t_bsAR->get_pointer(), nvirA, &(B_p_RR[0][0]), nvirA * ndf_,
                    0.0, C_p_AR->get_pointer(), nvirA * ndf_);
            C_DGEMM('N', 'T', noccA * nvirA, noccA * nvirA, ndf_, 1.0, &(B_p_AR[0][0]), ndf_, C_p_AR->get_pointer(), ndf_,
                    1.0, w_ARAR->get_pointer(), noccA * nvirA);

            for (int a = 0, ar = 0; a < noccA; a++) {
                for (int r = 0; r < nvirA; r++, ar++) {
                    for (int a1 = 0, a1r1 = 0; a1 < noccA; a1++) {
                        for (int r1 = 0; r1 < nvirA; r1++, a1r1++) {
                            int a1r = a1 * nvirA + r;
                            int ar1 = a * nvirA + r1;
                            double tval1 = w_ARARp[ar][a1r1] + w_ARARp[a1r1][ar];
                            double tval2 = w_ARARp[a1r][ar1] + w_ARARp[ar1][a1r];
                            double denom = evalsA[a + foccA] + evalsA[a1 + foccA] + evalsB[b + foccB] -
                                           evalsA[r + noccA + foccA] - evalsA[r1 + noccA + foccA] -
                                           evalsB[s + noccB + foccB];
                            energy += ((4.0 * tval1 - 2.0 * tval2) * tval1) / denom;
                        }
                    }
                }
            }
        }
        stop = std::time(nullptr);
        outfile->Printf("    (i = %3zu of %3zu) %10ld seconds\n", b + 1, noccB, stop - start);
    }

    free(B_p_bs);
    free_block(B_p_AA);
    free_block(B_p_AR);
    free_block(B_p_RR);

    return (energy);
}
}  // namespace sapt
}  // namespace psi
