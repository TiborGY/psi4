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

#include <cmath>
#include <ctime>
#include <memory>

#include "psi4/pragma.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/aiohandler.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"

#include "sapt0.h"
#include "sapt2.h"

namespace psi {
namespace sapt {

void SAPT0::ind20() {
    auto tAR = std::make_shared<Matrix>("tAR", noccA_, nvirA_);
    double **tARp = tAR->pointer();

    for (int a = 0; a < noccA_; a++) {
        for (int r = 0; r < nvirA_; r++) {
            tARp[a][r] = wBAR_[a][r] / (evalsA_[a] - evalsA_[r + noccA_]);
        }
    }

    double indA_B = 2.0 * C_DDOT(noccA_ * nvirA_, tAR->get_pointer(), 1, wBAR_[0], 1);

    if (no_response_) {
        CHFA_ = block_matrix(noccA_, nvirA_);  // TODO: Defer - class member requires API redesign
        C_DCOPY(noccA_ * nvirA_, tAR->get_pointer(), 1, CHFA_[0], 1);
    }

    auto tBS = std::make_shared<Matrix>("tBS", noccB_, nvirB_);
    double **tBSp = tBS->pointer();

    for (int b = 0; b < noccB_; b++) {
        for (int s = 0; s < nvirB_; s++) {
            tBSp[b][s] = wABS_[b][s] / (evalsB_[b] - evalsB_[s + noccB_]);
        }
    }

    double indB_A = 2.0 * C_DDOT(noccB_ * nvirB_, tBS->get_pointer(), 1, wABS_[0], 1);

    if (no_response_) {
        CHFB_ = block_matrix(noccB_, nvirB_);  // TODO: Defer - class member requires API redesign
        C_DCOPY(noccB_ * nvirB_, tBS->get_pointer(), 1, CHFB_[0], 1);
    }

    e_ind20_ = indA_B + indB_A;

    if (print_) {
        outfile->Printf("    Ind20 (A<-B)        = %18.12lf [Eh]\n", indA_B);
        outfile->Printf("    Ind20 (B<-A)        = %18.12lf [Eh]\n", indB_A);
        outfile->Printf("    Ind20               = %18.12lf [Eh]\n", e_ind20_);
    }
}

void SAPT0::ind20r() {
    if (aio_cphf_) {
        ind20rA_B_aio();
        ind20rB_A_aio();
    } else {
        ind20rA_B();
        ind20rB_A();
    }

    double indA_B, indB_A;

    indA_B = 2.0 * C_DDOT(noccA_ * nvirA_, CHFA_[0], 1, wBAR_[0], 1);
    indB_A = 2.0 * C_DDOT(noccB_ * nvirB_, CHFB_[0], 1, wABS_[0], 1);

    e_ind20_ = indA_B + indB_A;

    if (print_) {
        outfile->Printf("    Ind20,r (A<-B)      = %18.12lf [Eh]\n", indA_B);
        outfile->Printf("    Ind20,r (B<-A)      = %18.12lf [Eh]\n", indB_A);
        outfile->Printf("    Ind20,r             = %18.12lf [Eh]\n", e_ind20_);
    }
}

void SAPT0::ind20rA_B() {
    std::time_t start = std::time(nullptr);
    std::time_t stop;
    int iter = 0;
    double E_old, E;
    double conv, dE;
    double *tAR_old = init_array(noccA_ * nvirA_);
    double *tAR_new = init_array(noccA_ * nvirA_);

    double alpha, beta;
    double *R_old = init_array(noccA_ * nvirA_);
    double *R_new = init_array(noccA_ * nvirA_);
    double *z_old = init_array(noccA_ * nvirA_);
    double *z_new = init_array(noccA_ * nvirA_);
    double *Ax = init_array(noccA_ * nvirA_);

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif
    int rank = 0;

    for (int a = 0, ar = 0; a < noccA_; a++) {
        for (int r = 0; r < nvirA_; r++, ar++) {
            tAR_old[ar] = wBAR_[a][r] / (evalsA_[a] - evalsA_[r + noccA_]);
        }
    }

    E_old = 0.0;
    conv = 1.0;

    outfile->Printf("\n    Maxiter = %d\n", maxiter_);
    outfile->Printf("    CPHF R converge = %lE\n", cphf_r_conv_);

    if (print_) outfile->Printf("\n    Iter     Energy [mEh]          dE [mEh]         Residual      Time [s]\n");

    SAPTDFInts C_p_AA = set_C_AA();
    SAPTDFInts C_p_AR = set_C_AR();
    SAPTDFInts C_p_RR = set_C_RR();

    double *X = init_array(ndf_);
    auto xAA = std::make_shared<Matrix>("xAA", nthreads, noccA_ * noccA_);
    double **xAAp = xAA->pointer();
    auto xAR = std::make_shared<Matrix>("xAR", nthreads, noccA_ * nvirA_);
    double **xARp = xAR->pointer();
    auto xRR = std::make_shared<Matrix>("xRR", nthreads, nvirA_ * nvirA_);
    double **xRRp = xRR->pointer();
    auto tAR_dump = std::make_shared<Matrix>("tAR_dump", nthreads, noccA_ * nvirA_);
    double **tAR_dumpp = tAR_dump->pointer();

    do {
        memset(&(Ax[0]), '\0', sizeof(double) * noccA_ * nvirA_);
        memset(tAR_dumpp[0], '\0', sizeof(double) * nthreads * noccA_ * nvirA_);

        Iterator AR_iter = get_iterator(mem_, &C_p_AR);

        for (int i = 0, off = 0; i < AR_iter.num_blocks; i++) {
            read_block(&AR_iter, &C_p_AR);

            C_DGEMV('n', AR_iter.curr_size, noccA_ * nvirA_, 1.0, &(C_p_AR.B_p_[0][0]), noccA_ * nvirA_, tAR_old, 1,
                    0.0, &(X[0]), 1);
            C_DGEMV('t', AR_iter.curr_size, noccA_ * nvirA_, -4.0, &(C_p_AR.B_p_[0][0]), noccA_ * nvirA_, &(X[0]), 1,
                    1.0, Ax, 1);

#pragma omp parallel
            {
#pragma omp for private(rank)
                for (int j = 0; j < AR_iter.curr_size; j++) {
#ifdef _OPENMP
                    rank = omp_get_thread_num();
#endif

                    C_DGEMM('N', 'T', noccA_, noccA_, nvirA_, 1.0, &(C_p_AR.B_p_[j][0]), nvirA_, tAR_old, nvirA_, 0.0,
                            xAAp[rank], noccA_);
                    C_DGEMM('N', 'N', noccA_, nvirA_, noccA_, 1.0, xAAp[rank], noccA_, &(C_p_AR.B_p_[j][0]), nvirA_, 1.0,
                            tAR_dumpp[rank], nvirA_);
                }
            }
            off += AR_iter.curr_size;
        }

        C_p_AR.clear();

        Iterator RR_iter = get_iterator(mem_, &C_p_AA, &C_p_RR);

        for (int i = 0, off = 0; i < RR_iter.num_blocks; i++) {
            read_block(&RR_iter, &C_p_AA, &C_p_RR);

#pragma omp parallel
            {
#pragma omp for private(rank)
                for (int j = 0; j < RR_iter.curr_size; j++) {
#ifdef _OPENMP
                    rank = omp_get_thread_num();
#endif

                    for (int a = 0, ab = 0; a < nvirA_; a++) {
                        for (int b = 0; b <= a; b++) {
                            xRRp[rank][a * nvirA_ + b] = C_p_RR.B_p_[j][ab];
                            xRRp[rank][b * nvirA_ + a] = C_p_RR.B_p_[j][ab++];
                        }
                    }

                    C_DGEMM('N', 'N', noccA_, nvirA_, nvirA_, 1.0, tAR_old, nvirA_, &(xRRp[rank][0]), nvirA_, 0.0,
                            xARp[rank], nvirA_);
                    C_DGEMM('N', 'N', noccA_, nvirA_, noccA_, 1.0, &(C_p_AA.B_p_[j][0]), noccA_, xARp[rank], nvirA_, 1.0,
                            tAR_dumpp[rank], nvirA_);
                }
            }
            off += RR_iter.curr_size;
        }

        C_p_AA.clear();
        C_p_RR.clear();

        for (int n = 0; n < nthreads; n++) C_DAXPY(noccA_ * nvirA_, 1.0, tAR_dumpp[n], 1, Ax, 1);

        for (int a = 0, ar = 0; a < noccA_; a++) {
            for (int r = 0; r < nvirA_; r++, ar++) {
                Ax[ar] += (evalsA_[a] - evalsA_[r + noccA_]) * tAR_old[ar];
            }
        }

        if (!iter) {
            C_DCOPY(noccA_ * nvirA_, wBAR_[0], 1, R_old, 1);
            C_DAXPY(noccA_ * nvirA_, -1.0, Ax, 1, R_old, 1);
            for (int a = 0, ar = 0; a < noccA_; a++) {
                for (int r = 0; r < nvirA_; r++, ar++) {
                    z_old[ar] = R_old[ar] / (evalsA_[a] - evalsA_[r + noccA_]);
                }
            }
            C_DCOPY(noccA_ * nvirA_, tAR_old, 1, tAR_new, 1);
            C_DCOPY(noccA_ * nvirA_, z_old, 1, tAR_old, 1);
        } else {
            alpha = C_DDOT(noccA_ * nvirA_, R_old, 1, z_old, 1);
            alpha /= C_DDOT(noccA_ * nvirA_, tAR_old, 1, Ax, 1);
            C_DAXPY(noccA_ * nvirA_, alpha, tAR_old, 1, tAR_new, 1);
            C_DCOPY(noccA_ * nvirA_, R_old, 1, R_new, 1);
            C_DAXPY(noccA_ * nvirA_, -alpha, Ax, 1, R_new, 1);
            for (int a = 0, ar = 0; a < noccA_; a++) {
                for (int r = 0; r < nvirA_; r++, ar++) {
                    z_new[ar] = R_new[ar] / (evalsA_[a] - evalsA_[r + noccA_]);
                }
            }
            beta = C_DDOT(noccA_ * nvirA_, R_new, 1, z_new, 1);
            beta /= C_DDOT(noccA_ * nvirA_, R_old, 1, z_old, 1);
            C_DSCAL(noccA_ * nvirA_, beta, tAR_old, 1);
            C_DAXPY(noccA_ * nvirA_, 1.0, z_new, 1, tAR_old, 1);
            C_DCOPY(noccA_ * nvirA_, z_new, 1, z_old, 1);
            C_DCOPY(noccA_ * nvirA_, R_new, 1, R_old, 1);
        }

        E = 2.0 * C_DDOT(noccA_ * nvirA_, tAR_new, 1, &(wBAR_[0][0]), 1);

        conv = C_DDOT(noccA_ * nvirA_, R_old, 1, R_old, 1);
        conv = std::sqrt(conv);
        dE = E_old - E;

        iter++;
        stop = std::time(nullptr);
        if (print_) {
            outfile->Printf("    %4d %16.8lf %17.9lf %16.6e    %10ld\n", iter, E * 1000.0, dE * 1000.0, conv,
                            stop - start);
        }
        E_old = E;
    } while (conv > cphf_r_conv_ && iter < maxiter_);

    if (conv <= cphf_r_conv_) {
        if (print_) outfile->Printf("\n    CHF Iterations converged\n\n");
    } else {
        outfile->Printf("\n    CHF Iterations did not converge\n\n");
    }

    CHFA_ = block_matrix(noccA_, nvirA_);  // TODO: Defer - class member requires API redesign
    C_DCOPY(noccA_ * nvirA_, tAR_new, 1, CHFA_[0], 1);

    free(tAR_new);
    free(tAR_old);
    free(X);
    free(R_old);
    free(R_new);
    free(z_old);
    free(z_new);
    free(Ax);
}

void SAPT0::ind20rB_A() {
    std::time_t start = std::time(nullptr);
    std::time_t stop;
    int iter = 0;
    double E_old, E;
    double conv, dE;
    double *tBS_old = init_array(noccB_ * nvirB_);
    double *tBS_new = init_array(noccB_ * nvirB_);

    double alpha, beta;
    double *R_old = init_array(noccB_ * nvirB_);
    double *R_new = init_array(noccB_ * nvirB_);
    double *z_old = init_array(noccB_ * nvirB_);
    double *z_new = init_array(noccB_ * nvirB_);
    double *Ax = init_array(noccB_ * nvirB_);

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif
    int rank = 0;

    for (int b = 0, bs = 0; b < noccB_; b++) {
        for (int s = 0; s < nvirB_; s++, bs++) {
            tBS_old[bs] = wABS_[b][s] / (evalsB_[b] - evalsB_[s + noccB_]);
        }
    }

    E_old = 0.0;
    conv = 1.0;

    outfile->Printf("    Maxiter = %d\n", maxiter_);
    outfile->Printf("    CPHF R converge = %lE\n", cphf_r_conv_);

    if (print_) outfile->Printf("\n    Iter     Energy [mEh]          dE [mEh]         Residual      Time [s]\n");

    SAPTDFInts C_p_BB = set_C_BB();
    SAPTDFInts C_p_BS = set_C_BS();
    SAPTDFInts C_p_SS = set_C_SS();

    double *X = init_array(ndf_);
    auto xBB = std::make_shared<Matrix>("xBB", nthreads, noccB_ * noccB_);
    double **xBBp = xBB->pointer();
    auto xBS = std::make_shared<Matrix>("xBS", nthreads, noccB_ * nvirB_);
    double **xBSp = xBS->pointer();
    auto xSS = std::make_shared<Matrix>("xSS", nthreads, nvirB_ * nvirB_);
    double **xSSp = xSS->pointer();
    auto tBS_dump = std::make_shared<Matrix>("tBS_dump", nthreads, noccB_ * nvirB_);
    double **tBS_dumpp = tBS_dump->pointer();

    do {
        memset(&(Ax[0]), '\0', sizeof(double) * noccB_ * nvirB_);
        memset(tBS_dumpp[0], '\0', sizeof(double) * nthreads * noccB_ * nvirB_);

        Iterator BS_iter = get_iterator(mem_, &C_p_BS);

        for (int i = 0, off = 0; i < BS_iter.num_blocks; i++) {
            read_block(&BS_iter, &C_p_BS);

            C_DGEMV('n', BS_iter.curr_size, noccB_ * nvirB_, 1.0, &(C_p_BS.B_p_[0][0]), noccB_ * nvirB_, tBS_old, 1,
                    0.0, &(X[0]), 1);
            C_DGEMV('t', BS_iter.curr_size, noccB_ * nvirB_, -4.0, &(C_p_BS.B_p_[0][0]), noccB_ * nvirB_, &(X[0]), 1,
                    1.0, Ax, 1);

#pragma omp parallel
            {
#pragma omp for private(rank)
                for (int j = 0; j < BS_iter.curr_size; j++) {
#ifdef _OPENMP
                    rank = omp_get_thread_num();
#endif

                    C_DGEMM('N', 'T', noccB_, noccB_, nvirB_, 1.0, &(C_p_BS.B_p_[j][0]), nvirB_, tBS_old, nvirB_, 0.0,
                            xBBp[rank], noccB_);
                    C_DGEMM('N', 'N', noccB_, nvirB_, noccB_, 1.0, xBBp[rank], noccB_, &(C_p_BS.B_p_[j][0]), nvirB_, 1.0,
                            tBS_dumpp[rank], nvirB_);
                }
            }
            off += BS_iter.curr_size;
        }

        C_p_BS.clear();

        Iterator SS_iter = get_iterator(mem_, &C_p_BB, &C_p_SS);

        for (int i = 0, off = 0; i < SS_iter.num_blocks; i++) {
            read_block(&SS_iter, &C_p_BB, &C_p_SS);

#pragma omp parallel
            {
#pragma omp for private(rank)
                for (int j = 0; j < SS_iter.curr_size; j++) {
#ifdef _OPENMP
                    rank = omp_get_thread_num();
#endif

                    for (int a = 0, ab = 0; a < nvirB_; a++) {
                        for (int b = 0; b <= a; b++) {
                            xSSp[rank][a * nvirB_ + b] = C_p_SS.B_p_[j][ab];
                            xSSp[rank][b * nvirB_ + a] = C_p_SS.B_p_[j][ab++];
                        }
                    }

                    C_DGEMM('N', 'N', noccB_, nvirB_, nvirB_, 1.0, tBS_old, nvirB_, &(xSSp[rank][0]), nvirB_, 0.0,
                            xBSp[rank], nvirB_);
                    C_DGEMM('N', 'N', noccB_, nvirB_, noccB_, 1.0, &(C_p_BB.B_p_[j][0]), noccB_, xBSp[rank], nvirB_, 1.0,
                            tBS_dumpp[rank], nvirB_);
                }
            }
            off += SS_iter.curr_size;
        }

        C_p_BB.clear();
        C_p_SS.clear();

        for (int n = 0; n < nthreads; n++) C_DAXPY(noccB_ * nvirB_, 1.0, tBS_dumpp[n], 1, Ax, 1);

        for (int b = 0, bs = 0; b < noccB_; b++) {
            for (int s = 0; s < nvirB_; s++, bs++) {
                Ax[bs] += (evalsB_[b] - evalsB_[s + noccB_]) * tBS_old[bs];
            }
        }

        if (!iter) {
            C_DCOPY(noccB_ * nvirB_, wABS_[0], 1, R_old, 1);
            C_DAXPY(noccB_ * nvirB_, -1.0, Ax, 1, R_old, 1);
            for (int b = 0, bs = 0; b < noccB_; b++) {
                for (int s = 0; s < nvirB_; s++, bs++) {
                    z_old[bs] = R_old[bs] / (evalsB_[b] - evalsB_[s + noccB_]);
                }
            }
            C_DCOPY(noccB_ * nvirB_, tBS_old, 1, tBS_new, 1);
            C_DCOPY(noccB_ * nvirB_, z_old, 1, tBS_old, 1);
        } else {
            alpha = C_DDOT(noccB_ * nvirB_, R_old, 1, z_old, 1);
            alpha /= C_DDOT(noccB_ * nvirB_, tBS_old, 1, Ax, 1);
            C_DAXPY(noccB_ * nvirB_, alpha, tBS_old, 1, tBS_new, 1);
            C_DCOPY(noccB_ * nvirB_, R_old, 1, R_new, 1);
            C_DAXPY(noccB_ * nvirB_, -alpha, Ax, 1, R_new, 1);
            for (int b = 0, bs = 0; b < noccB_; b++) {
                for (int s = 0; s < nvirB_; s++, bs++) {
                    z_new[bs] = R_new[bs] / (evalsB_[b] - evalsB_[s + noccB_]);
                }
            }
            beta = C_DDOT(noccB_ * nvirB_, R_new, 1, z_new, 1);
            beta /= C_DDOT(noccB_ * nvirB_, R_old, 1, z_old, 1);
            C_DSCAL(noccB_ * nvirB_, beta, tBS_old, 1);
            C_DAXPY(noccB_ * nvirB_, 1.0, z_new, 1, tBS_old, 1);
            C_DCOPY(noccB_ * nvirB_, z_new, 1, z_old, 1);
            C_DCOPY(noccB_ * nvirB_, R_new, 1, R_old, 1);
        }

        E = 2.0 * C_DDOT(noccB_ * nvirB_, tBS_new, 1, &(wABS_[0][0]), 1);

        conv = C_DDOT(noccB_ * nvirB_, R_old, 1, R_old, 1);
        conv = std::sqrt(conv);
        dE = E_old - E;

        iter++;
        stop = std::time(nullptr);
        if (print_) {
            outfile->Printf("    %4d %16.8lf %17.9lf %16.6e    %10ld\n", iter, E * 1000.0, dE * 1000.0, conv,
                            stop - start);
        }
        E_old = E;
    } while (conv > cphf_r_conv_ && iter < maxiter_);

    if (conv <= cphf_r_conv_) {
        if (print_) outfile->Printf("\n    CHF Iterations converged\n\n");
    } else {
        outfile->Printf("\n    CHF Iterations did not converge\n\n");
    }

    CHFB_ = block_matrix(noccB_, nvirB_);  // TODO: Defer - class member requires API redesign
    C_DCOPY(noccB_ * nvirB_, tBS_new, 1, CHFB_[0], 1);

    free(tBS_new);
    free(tBS_old);
    free(X);
    free(R_old);
    free(R_new);
    free(z_old);
    free(z_new);
    free(Ax);
}

void SAPT0::ind20rA_B_aio() {
    std::time_t start = std::time(nullptr);
    std::time_t stop;
    int iter = 0;
    double E_old, E;
    double conv, dE;
    double *tAR_old = init_array(noccA_ * nvirA_);
    double *tAR_new = init_array(noccA_ * nvirA_);

    double alpha, beta;
    double *R_old = init_array(noccA_ * nvirA_);
    double *R_new = init_array(noccA_ * nvirA_);
    double *z_old = init_array(noccA_ * nvirA_);
    double *z_new = init_array(noccA_ * nvirA_);
    double *Ax = init_array(noccA_ * nvirA_);

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif
    int rank = 0;

    for (int a = 0, ar = 0; a < noccA_; a++) {
        for (int r = 0; r < nvirA_; r++, ar++) {
            tAR_old[ar] = wBAR_[a][r] / (evalsA_[a] - evalsA_[r + noccA_]);
        }
    }

    E_old = 0.0;
    conv = 1.0;

    outfile->Printf("\n    Maxiter = %d\n", maxiter_);
    outfile->Printf("    CPHF R converge = %lE\n", cphf_r_conv_);

    if (print_) outfile->Printf("\n    Iter     Energy [mEh]          dE [mEh]         Residual      Time [s]\n");

    SAPTDFInts C_p_AR = set_C_AR();

    double *X = init_array(ndf_);
    auto xAA = std::make_shared<Matrix>("xAA", nthreads, noccA_ * noccA_);
    double **xAAp = xAA->pointer();
    auto xAR = std::make_shared<Matrix>("xAR", nthreads, noccA_ * nvirA_);
    double **xARp = xAR->pointer();
    auto xRR = std::make_shared<Matrix>("xRR", nthreads, nvirA_ * nvirA_);
    double **xRRp = xRR->pointer();
    auto tAR_dump = std::make_shared<Matrix>("tAR_dump", nthreads, noccA_ * nvirA_);
    double **tAR_dumpp = tAR_dump->pointer();

    long int block_length = mem_ / (2L * noccA_ * noccA_ + 2L * nvirA_ * (nvirA_ + 1) / 2);
    if (block_length > ndf_) block_length = ndf_;

    int num_blocks = ndf_ / block_length;
    if (ndf_ % block_length) num_blocks++;

    auto aio = std::make_shared<AIOHandler>(psio_);

    std::shared_ptr<Matrix> C_p_AA[2];
    double **C_p_AAp[2];
    std::shared_ptr<Matrix> C_p_RR[2];
    double **C_p_RRp[2];

    psio_address next_C_p_AA;
    psio_address next_C_p_RR;

    do {
        memset(&(Ax[0]), '\0', sizeof(double) * noccA_ * nvirA_);
        memset(tAR_dumpp[0], '\0', sizeof(double) * nthreads * noccA_ * nvirA_);

        Iterator AR_iter = get_iterator(mem_, &C_p_AR);

        for (int i = 0, off = 0; i < AR_iter.num_blocks; i++) {
            read_block(&AR_iter, &C_p_AR);

            C_DGEMV('n', AR_iter.curr_size, noccA_ * nvirA_, 1.0, &(C_p_AR.B_p_[0][0]), noccA_ * nvirA_, tAR_old, 1,
                    0.0, &(X[0]), 1);
            C_DGEMV('t', AR_iter.curr_size, noccA_ * nvirA_, -4.0, &(C_p_AR.B_p_[0][0]), noccA_ * nvirA_, &(X[0]), 1,
                    1.0, Ax, 1);

#pragma omp parallel
            {
#pragma omp for private(rank)
                for (int j = 0; j < AR_iter.curr_size; j++) {
#ifdef _OPENMP
                    rank = omp_get_thread_num();
#endif

                    C_DGEMM('N', 'T', noccA_, noccA_, nvirA_, 1.0, &(C_p_AR.B_p_[j][0]), nvirA_, tAR_old, nvirA_, 0.0,
                            xAAp[rank], noccA_);
                    C_DGEMM('N', 'N', noccA_, nvirA_, noccA_, 1.0, xAAp[rank], noccA_, &(C_p_AR.B_p_[j][0]), nvirA_, 1.0,
                            tAR_dumpp[rank], nvirA_);
                }
            }
            off += AR_iter.curr_size;
        }

        C_p_AR.clear();

        C_p_AA[0] = std::make_shared<Matrix>("C_p_AA[0]", block_length, noccA_ * noccA_);
        C_p_AAp[0] = C_p_AA[0]->pointer();
        C_p_AA[1] = std::make_shared<Matrix>("C_p_AA[1]", block_length, noccA_ * noccA_);
        C_p_AAp[1] = C_p_AA[1]->pointer();
        C_p_RR[0] = std::make_shared<Matrix>("C_p_RR[0]", block_length, nvirA_ * (nvirA_ + 1) / 2);
        C_p_RRp[0] = C_p_RR[0]->pointer();
        C_p_RR[1] = std::make_shared<Matrix>("C_p_RR[1]", block_length, nvirA_ * (nvirA_ + 1) / 2);
        C_p_RRp[1] = C_p_RR[1]->pointer();

        next_C_p_AA = PSIO_ZERO;
        next_C_p_RR = PSIO_ZERO;

        psio_->read(PSIF_SAPT_AA_DF_INTS, "AA RI Integrals", (char *)C_p_AA[0]->get_pointer(),
                    sizeof(double) * block_length * noccA_ * noccA_, next_C_p_AA, &next_C_p_AA);
        psio_->read(PSIF_SAPT_AA_DF_INTS, "RR RI Integrals", (char *)C_p_RR[0]->get_pointer(),
                    sizeof(double) * block_length * nvirA_ * (nvirA_ + 1) / 2, next_C_p_RR, &next_C_p_RR);

        for (int i = 0; i < num_blocks; i++) {
            if (i < num_blocks - 1) {
                int read_length = block_length;
                if (i == num_blocks - 2 && ndf_ % block_length) read_length = ndf_ % block_length;
                aio->read(PSIF_SAPT_AA_DF_INTS, "AA RI Integrals", (char *)C_p_AA[(i + 1) % 2]->get_pointer(),
                          sizeof(double) * read_length * noccA_ * noccA_, next_C_p_AA, &next_C_p_AA);
                aio->read(PSIF_SAPT_AA_DF_INTS, "RR RI Integrals", (char *)C_p_RR[(i + 1) % 2]->get_pointer(),
                          sizeof(double) * read_length * nvirA_ * (nvirA_ + 1) / 2, next_C_p_RR, &next_C_p_RR);
            }

            int loopsize = block_length;
            if (i == num_blocks - 1 && ndf_ % block_length) loopsize = ndf_ % block_length;

#pragma omp parallel
            {
#pragma omp for private(rank)
                for (int j = 0; j < loopsize; j++) {
#ifdef _OPENMP
                    rank = omp_get_thread_num();
#endif

                    for (int a = 0, ab = 0; a < nvirA_; a++) {
                        for (int b = 0; b <= a; b++) {
                            xRRp[rank][a * nvirA_ + b] = C_p_RRp[i % 2][j][ab];
                            xRRp[rank][b * nvirA_ + a] = C_p_RRp[i % 2][j][ab++];
                        }
                    }

                    C_DGEMM('N', 'N', noccA_, nvirA_, nvirA_, 1.0, tAR_old, nvirA_, &(xRRp[rank][0]), nvirA_, 0.0,
                            xARp[rank], nvirA_);
                    C_DGEMM('N', 'N', noccA_, nvirA_, noccA_, 1.0, &(C_p_AAp[i % 2][j][0]), noccA_, xARp[rank], nvirA_,
                            1.0, tAR_dumpp[rank], nvirA_);
                }
            }

            if (i < num_blocks - 1) aio->synchronize();
        }

        for (int n = 0; n < nthreads; n++) C_DAXPY(noccA_ * nvirA_, 1.0, tAR_dumpp[n], 1, Ax, 1);

        for (int a = 0, ar = 0; a < noccA_; a++) {
            for (int r = 0; r < nvirA_; r++, ar++) {
                Ax[ar] += (evalsA_[a] - evalsA_[r + noccA_]) * tAR_old[ar];
            }
        }

        if (!iter) {
            C_DCOPY(noccA_ * nvirA_, wBAR_[0], 1, R_old, 1);
            C_DAXPY(noccA_ * nvirA_, -1.0, Ax, 1, R_old, 1);
            for (int a = 0, ar = 0; a < noccA_; a++) {
                for (int r = 0; r < nvirA_; r++, ar++) {
                    z_old[ar] = R_old[ar] / (evalsA_[a] - evalsA_[r + noccA_]);
                }
            }
            C_DCOPY(noccA_ * nvirA_, tAR_old, 1, tAR_new, 1);
            C_DCOPY(noccA_ * nvirA_, z_old, 1, tAR_old, 1);
        } else {
            alpha = C_DDOT(noccA_ * nvirA_, R_old, 1, z_old, 1);
            alpha /= C_DDOT(noccA_ * nvirA_, tAR_old, 1, Ax, 1);
            C_DAXPY(noccA_ * nvirA_, alpha, tAR_old, 1, tAR_new, 1);
            C_DCOPY(noccA_ * nvirA_, R_old, 1, R_new, 1);
            C_DAXPY(noccA_ * nvirA_, -alpha, Ax, 1, R_new, 1);
            for (int a = 0, ar = 0; a < noccA_; a++) {
                for (int r = 0; r < nvirA_; r++, ar++) {
                    z_new[ar] = R_new[ar] / (evalsA_[a] - evalsA_[r + noccA_]);
                }
            }
            beta = C_DDOT(noccA_ * nvirA_, R_new, 1, z_new, 1);
            beta /= C_DDOT(noccA_ * nvirA_, R_old, 1, z_old, 1);
            C_DSCAL(noccA_ * nvirA_, beta, tAR_old, 1);
            C_DAXPY(noccA_ * nvirA_, 1.0, z_new, 1, tAR_old, 1);
            C_DCOPY(noccA_ * nvirA_, z_new, 1, z_old, 1);
            C_DCOPY(noccA_ * nvirA_, R_new, 1, R_old, 1);
        }

        E = 2.0 * C_DDOT(noccA_ * nvirA_, tAR_new, 1, &(wBAR_[0][0]), 1);

        conv = C_DDOT(noccA_ * nvirA_, R_old, 1, R_old, 1);
        conv = std::sqrt(conv);
        dE = E_old - E;

        iter++;
        stop = std::time(nullptr);
        if (print_) {
            outfile->Printf("    %4d %16.8lf %17.9lf %16.6e    %10ld\n", iter, E * 1000.0, dE * 1000.0, conv,
                            stop - start);
        }
        E_old = E;
    } while (conv > cphf_r_conv_ && iter < maxiter_);

    if (conv <= cphf_r_conv_) {
        if (print_) outfile->Printf("\n    CHF Iterations converged\n\n");
    } else {
        outfile->Printf("\n    CHF Iterations did not converge\n\n");
    }

    CHFA_ = block_matrix(noccA_, nvirA_);  // TODO: Defer - class member requires API redesign
    C_DCOPY(noccA_ * nvirA_, tAR_new, 1, CHFA_[0], 1);

    free(tAR_new);
    free(tAR_old);
    free(X);
    free(R_old);
    free(R_new);
    free(z_old);
    free(z_new);
    free(Ax);
}

void SAPT0::ind20rB_A_aio() {
    std::time_t start = std::time(nullptr);
    std::time_t stop;
    int iter = 0;
    double E_old, E;
    double conv, dE;
    double *tBS_old = init_array(noccB_ * nvirB_);
    double *tBS_new = init_array(noccB_ * nvirB_);

    double alpha, beta;
    double *R_old = init_array(noccB_ * nvirB_);
    double *R_new = init_array(noccB_ * nvirB_);
    double *z_old = init_array(noccB_ * nvirB_);
    double *z_new = init_array(noccB_ * nvirB_);
    double *Ax = init_array(noccB_ * nvirB_);

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif
    int rank = 0;

    for (int b = 0, bs = 0; b < noccB_; b++) {
        for (int s = 0; s < nvirB_; s++, bs++) {
            tBS_old[bs] = wABS_[b][s] / (evalsB_[b] - evalsB_[s + noccB_]);
        }
    }

    E_old = 0.0;
    conv = 1.0;

    outfile->Printf("    Maxiter = %d\n", maxiter_);
    outfile->Printf("    CPHF R converge = %lE\n", cphf_r_conv_);

    if (print_) outfile->Printf("\n    Iter     Energy [mEh]          dE [mEh]         Residual      Time [s]\n");

    SAPTDFInts C_p_BS = set_C_BS();

    double *X = init_array(ndf_);
    auto xBB = std::make_shared<Matrix>("xBB", nthreads, noccB_ * noccB_);
    double **xBBp = xBB->pointer();
    auto xBS = std::make_shared<Matrix>("xBS", nthreads, noccB_ * nvirB_);
    double **xBSp = xBS->pointer();
    auto xSS = std::make_shared<Matrix>("xSS", nthreads, nvirB_ * nvirB_);
    double **xSSp = xSS->pointer();
    auto tBS_dump = std::make_shared<Matrix>("tBS_dump", nthreads, noccB_ * nvirB_);
    double **tBS_dumpp = tBS_dump->pointer();

    long int block_length = mem_ / (2L * noccB_ * noccB_ + 2L * nvirB_ * (nvirB_ + 1) / 2);
    if (block_length > ndf_) block_length = ndf_;

    int num_blocks = ndf_ / block_length;
    if (ndf_ % block_length) num_blocks++;

    auto aio = std::make_shared<AIOHandler>(psio_);

    std::shared_ptr<Matrix> C_p_BB[2];
    double **C_p_BBp[2];
    std::shared_ptr<Matrix> C_p_SS[2];
    double **C_p_SSp[2];

    psio_address next_C_p_BB;
    psio_address next_C_p_SS;

    do {
        memset(&(Ax[0]), '\0', sizeof(double) * noccB_ * nvirB_);
        memset(tBS_dumpp[0], '\0', sizeof(double) * nthreads * noccB_ * nvirB_);

        Iterator BS_iter = get_iterator(mem_, &C_p_BS);

        for (int i = 0, off = 0; i < BS_iter.num_blocks; i++) {
            read_block(&BS_iter, &C_p_BS);

            C_DGEMV('n', BS_iter.curr_size, noccB_ * nvirB_, 1.0, &(C_p_BS.B_p_[0][0]), noccB_ * nvirB_, tBS_old, 1,
                    0.0, &(X[0]), 1);
            C_DGEMV('t', BS_iter.curr_size, noccB_ * nvirB_, -4.0, &(C_p_BS.B_p_[0][0]), noccB_ * nvirB_, &(X[0]), 1,
                    1.0, Ax, 1);

#pragma omp parallel
            {
#pragma omp for private(rank)
                for (int j = 0; j < BS_iter.curr_size; j++) {
#ifdef _OPENMP
                    rank = omp_get_thread_num();
#endif

                    C_DGEMM('N', 'T', noccB_, noccB_, nvirB_, 1.0, &(C_p_BS.B_p_[j][0]), nvirB_, tBS_old, nvirB_, 0.0,
                            xBBp[rank], noccB_);
                    C_DGEMM('N', 'N', noccB_, nvirB_, noccB_, 1.0, xBBp[rank], noccB_, &(C_p_BS.B_p_[j][0]), nvirB_, 1.0,
                            tBS_dumpp[rank], nvirB_);
                }
            }
            off += BS_iter.curr_size;
        }

        C_p_BS.clear();

        C_p_BB[0] = std::make_shared<Matrix>("C_p_BB[0]", block_length, noccB_ * noccB_);
        C_p_BBp[0] = C_p_BB[0]->pointer();
        C_p_BB[1] = std::make_shared<Matrix>("C_p_BB[1]", block_length, noccB_ * noccB_);
        C_p_BBp[1] = C_p_BB[1]->pointer();
        C_p_SS[0] = std::make_shared<Matrix>("C_p_SS[0]", block_length, nvirB_ * (nvirB_ + 1) / 2);
        C_p_SSp[0] = C_p_SS[0]->pointer();
        C_p_SS[1] = std::make_shared<Matrix>("C_p_SS[1]", block_length, nvirB_ * (nvirB_ + 1) / 2);
        C_p_SSp[1] = C_p_SS[1]->pointer();

        next_C_p_BB = PSIO_ZERO;
        next_C_p_SS = PSIO_ZERO;

        psio_->read(PSIF_SAPT_BB_DF_INTS, "BB RI Integrals", (char *)C_p_BB[0]->get_pointer(),
                    sizeof(double) * block_length * noccB_ * noccB_, next_C_p_BB, &next_C_p_BB);
        psio_->read(PSIF_SAPT_BB_DF_INTS, "SS RI Integrals", (char *)C_p_SS[0]->get_pointer(),
                    sizeof(double) * block_length * nvirB_ * (nvirB_ + 1) / 2, next_C_p_SS, &next_C_p_SS);

        for (int i = 0; i < num_blocks; i++) {
            if (i < num_blocks - 1) {
                int read_length = block_length;
                if (i == num_blocks - 2 && ndf_ % block_length) read_length = ndf_ % block_length;
                aio->read(PSIF_SAPT_BB_DF_INTS, "BB RI Integrals", (char *)C_p_BB[(i + 1) % 2]->get_pointer(),
                          sizeof(double) * read_length * noccB_ * noccB_, next_C_p_BB, &next_C_p_BB);
                aio->read(PSIF_SAPT_BB_DF_INTS, "SS RI Integrals", (char *)C_p_SS[(i + 1) % 2]->get_pointer(),
                          sizeof(double) * read_length * nvirB_ * (nvirB_ + 1) / 2, next_C_p_SS, &next_C_p_SS);
            }

            int loopsize = block_length;
            if (i == num_blocks - 1 && ndf_ % block_length) loopsize = ndf_ % block_length;

#pragma omp parallel
            {
#pragma omp for private(rank)
                for (int j = 0; j < loopsize; j++) {
#ifdef _OPENMP
                    rank = omp_get_thread_num();
#endif

                    for (int a = 0, ab = 0; a < nvirB_; a++) {
                        for (int b = 0; b <= a; b++) {
                            xSSp[rank][a * nvirB_ + b] = C_p_SSp[i % 2][j][ab];
                            xSSp[rank][b * nvirB_ + a] = C_p_SSp[i % 2][j][ab++];
                        }
                    }

                    C_DGEMM('N', 'N', noccB_, nvirB_, nvirB_, 1.0, tBS_old, nvirB_, &(xSSp[rank][0]), nvirB_, 0.0,
                            xBSp[rank], nvirB_);
                    C_DGEMM('N', 'N', noccB_, nvirB_, noccB_, 1.0, &(C_p_BBp[i % 2][j][0]), noccB_, xBSp[rank], nvirB_,
                            1.0, tBS_dumpp[rank], nvirB_);
                }
            }

            if (i < num_blocks - 1) aio->synchronize();
        }

        for (int n = 0; n < nthreads; n++) C_DAXPY(noccB_ * nvirB_, 1.0, tBS_dumpp[n], 1, Ax, 1);

        for (int b = 0, bs = 0; b < noccB_; b++) {
            for (int s = 0; s < nvirB_; s++, bs++) {
                Ax[bs] += (evalsB_[b] - evalsB_[s + noccB_]) * tBS_old[bs];
            }
        }

        if (!iter) {
            C_DCOPY(noccB_ * nvirB_, wABS_[0], 1, R_old, 1);
            C_DAXPY(noccB_ * nvirB_, -1.0, Ax, 1, R_old, 1);
            for (int b = 0, bs = 0; b < noccB_; b++) {
                for (int s = 0; s < nvirB_; s++, bs++) {
                    z_old[bs] = R_old[bs] / (evalsB_[b] - evalsB_[s + noccB_]);
                }
            }
            C_DCOPY(noccB_ * nvirB_, tBS_old, 1, tBS_new, 1);
            C_DCOPY(noccB_ * nvirB_, z_old, 1, tBS_old, 1);
        } else {
            alpha = C_DDOT(noccB_ * nvirB_, R_old, 1, z_old, 1);
            alpha /= C_DDOT(noccB_ * nvirB_, tBS_old, 1, Ax, 1);
            C_DAXPY(noccB_ * nvirB_, alpha, tBS_old, 1, tBS_new, 1);
            C_DCOPY(noccB_ * nvirB_, R_old, 1, R_new, 1);
            C_DAXPY(noccB_ * nvirB_, -alpha, Ax, 1, R_new, 1);
            for (int b = 0, bs = 0; b < noccB_; b++) {
                for (int s = 0; s < nvirB_; s++, bs++) {
                    z_new[bs] = R_new[bs] / (evalsB_[b] - evalsB_[s + noccB_]);
                }
            }
            beta = C_DDOT(noccB_ * nvirB_, R_new, 1, z_new, 1);
            beta /= C_DDOT(noccB_ * nvirB_, R_old, 1, z_old, 1);
            C_DSCAL(noccB_ * nvirB_, beta, tBS_old, 1);
            C_DAXPY(noccB_ * nvirB_, 1.0, z_new, 1, tBS_old, 1);
            C_DCOPY(noccB_ * nvirB_, z_new, 1, z_old, 1);
            C_DCOPY(noccB_ * nvirB_, R_new, 1, R_old, 1);
        }

        E = 2.0 * C_DDOT(noccB_ * nvirB_, tBS_new, 1, &(wABS_[0][0]), 1);

        conv = C_DDOT(noccB_ * nvirB_, R_old, 1, R_old, 1);
        conv = std::sqrt(conv);
        dE = E_old - E;

        iter++;
        stop = std::time(nullptr);
        if (print_) {
            outfile->Printf("    %4d %16.8lf %17.9lf %16.6e    %10ld\n", iter, E * 1000.0, dE * 1000.0, conv,
                            stop - start);
        }
        E_old = E;
    } while (conv > cphf_r_conv_ && iter < maxiter_);

    if (conv <= cphf_r_conv_) {
        if (print_) outfile->Printf("\n    CHF Iterations converged\n\n");
    } else {
        outfile->Printf("\n    CHF Iterations did not converge\n\n");
    }

    CHFB_ = block_matrix(noccB_, nvirB_);  // TODO: Defer - class member requires API redesign
    C_DCOPY(noccB_ * nvirB_, tBS_new, 1, CHFB_[0], 1);

    free(tBS_new);
    free(tBS_old);
    free(X);
    free(R_old);
    free(R_new);
    free(z_old);
    free(z_new);
    free(Ax);
}

void SAPT2::ind20r() {
    CHFA_ = block_matrix(noccA_, nvirA_);  // TODO: Defer - class member requires API redesign

    cphf_solver(CHFA_, wBAR_, evalsA_, PSIF_SAPT_AA_DF_INTS, "AA RI Integrals", "AR RI Integrals", "RR RI Integrals",
                noccA_, nvirA_);

    CHFB_ = block_matrix(noccB_, nvirB_);  // TODO: Defer - class member requires API redesign

    cphf_solver(CHFB_, wABS_, evalsB_, PSIF_SAPT_BB_DF_INTS, "BB RI Integrals", "BS RI Integrals", "SS RI Integrals",
                noccB_, nvirB_);

    double indA_B, indB_A;

    indA_B = 2.0 * C_DDOT(noccA_ * nvirA_, CHFA_[0], 1, wBAR_[0], 1);
    indB_A = 2.0 * C_DDOT(noccB_ * nvirB_, CHFB_[0], 1, wABS_[0], 1);

    e_ind20_ = indA_B + indB_A;

    if (print_) {
        outfile->Printf("    Ind20,r (A<-B)      = %18.12lf [Eh]\n", indA_B);
        outfile->Printf("    Ind20,r (B<-A)      = %18.12lf [Eh]\n", indB_A);
        outfile->Printf("    Ind20,r             = %18.12lf [Eh]\n", e_ind20_);
    }
}

void SAPT2::cphf_solver(double **tAR, double **wBAR, double *evals, int intfile, const char *AAints, const char *ARints,
                        const char *RRints, size_t nocc, size_t nvir) {
    auto B_p_AR = std::make_shared<Matrix>("B_p_AR", nocc * nvir, ndf_ + 3);
    double **B_p_ARp = B_p_AR->pointer();

    psio_->read_entry(intfile, ARints, (char *)B_p_AR->get_pointer(), sizeof(double) * nocc * nvir * (ndf_ + 3));

    auto Amat = std::make_shared<Matrix>("Amat", nocc * nvir, nocc * nvir);
    double **Amatp = Amat->pointer();

    C_DGEMM('N', 'T', nocc * nvir, nocc * nvir, ndf_, -4.0, B_p_ARp[0], ndf_ + 3, B_p_ARp[0], ndf_ + 3, 0.0, Amatp[0],
            nocc * nvir);

    for (int a = 0, ar = 0; a < nocc; a++) {
        for (int r = 0; r < nvir; r++, ar++) {
            C_DGEMM('N', 'T', nocc, nvir, ndf_, 1.0, B_p_ARp[r], nvir * (ndf_ + 3), B_p_ARp[a * nvir], ndf_ + 3, 1.0,
                    Amatp[ar], nvir);
        }
    }

    auto B_p_AA = std::make_shared<Matrix>("B_p_AA", nocc * nocc, ndf_ + 3);
    double **B_p_AAp = B_p_AA->pointer();
    auto B_p_R = std::make_shared<Matrix>("B_p_R", nvir, ndf_ + 3);

    psio_->read_entry(intfile, AAints, (char *)B_p_AA->get_pointer(), sizeof(double) * nocc * nocc * (ndf_ + 3));

    psio_address next_PSIF = PSIO_ZERO;

    for (int r = 0; r < nvir; r++) {
        psio_->read(intfile, RRints, (char *)B_p_R->get_pointer(), sizeof(double) * nvir * (ndf_ + 3), next_PSIF, &next_PSIF);
        for (int a = 0; a < nocc; a++) {
            int ar = a * nvir + r;
            C_DGEMM('N', 'T', nocc, nvir, ndf_, 1.0, B_p_AAp[a * nocc], ndf_ + 3, B_p_R->get_pointer(), ndf_ + 3, 1.0, Amatp[ar],
                    nvir);
        }
    }

    for (int a = 0, ar = 0; a < nocc; a++) {
        for (int r = 0; r < nvir; r++, ar++) {
            Amatp[ar][ar] += (evals[a] - evals[r + nocc]);
        }
    }

    int *ipiv = init_int_array(nocc * nvir);

    C_DCOPY(nocc * nvir, wBAR[0], 1, tAR[0], 1);
    C_DGESV(nocc * nvir, 1, Amatp[0], nocc * nvir, ipiv, tAR[0], nocc * nvir);

    free(ipiv);
}
}  // namespace sapt
}  // namespace psi
