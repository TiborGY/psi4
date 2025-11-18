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

/**
 * @file hamiltonian_usapt.cc
 * @brief Unrestricted CPKS Hamiltonian for USAPT0 calculations
 *
 * This file implements the CPKSUSAPTHamiltonian class, which provides
 * Hamiltonian-vector products for unrestricted SAPT0 coupled response calculations.
 * The implementation consolidates product logic previously duplicated in
 * CPKS_USAPT0::product() (libsapt_solver/usapt0.cc).
 *
 * Key features:
 * - Spin-unrestricted two-monomer problem (4 response vectors: Aa, Ab, Ba, Bb)
 * - Map-based interface for flexible perturbation handling
 * - Product formula: (2J_total - K - K^T + diagonal) where J_total = J_alpha + J_beta
 * - Handles systems with no alpha/beta electrons gracefully
 * - Static preconditioner helper for Jacobi preconditioning
 *
 * The key difference from restricted CPHF is the 2J scaling (vs 4J) which arises
 * from the unrestricted formalism where alpha and beta densities are treated
 * independently but couple through the total Coulomb potential.
 *
 * This class is used internally by CPKS_USAPT0 for backwards compatibility
 * while providing a clean, reusable implementation of the product logic.
 *
 * @see libfock/CPHF_ARCHITECTURE.md for architecture documentation
 * @see hamiltonian_usapt.h for class interface
 */

#include "hamiltonian_usapt.h"
#include "jk.h"

#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {

CPKSUSAPTHamiltonian::CPKSUSAPTHamiltonian(
    std::shared_ptr<JK> jk,
    SharedMatrix Cocca_A, SharedMatrix Cvira_A,
    std::shared_ptr<Vector> eps_occa_A, std::shared_ptr<Vector> eps_vira_A,
    SharedMatrix Coccb_A, SharedMatrix Cvirb_A,
    std::shared_ptr<Vector> eps_occb_A, std::shared_ptr<Vector> eps_virb_A,
    SharedMatrix Cocca_B, SharedMatrix Cvira_B,
    std::shared_ptr<Vector> eps_occa_B, std::shared_ptr<Vector> eps_vira_B,
    SharedMatrix Coccb_B, SharedMatrix Cvirb_B,
    std::shared_ptr<Vector> eps_occb_B, std::shared_ptr<Vector> eps_virb_B)
    : RHamiltonian(jk),
      Cocca_A_(Cocca_A), Cvira_A_(Cvira_A), eps_occa_A_(eps_occa_A), eps_vira_A_(eps_vira_A),
      Coccb_A_(Coccb_A), Cvirb_A_(Cvirb_A), eps_occb_A_(eps_occb_A), eps_virb_A_(eps_virb_A),
      Cocca_B_(Cocca_B), Cvira_B_(Cvira_B), eps_occa_B_(eps_occa_B), eps_vira_B_(eps_vira_B),
      Coccb_B_(Coccb_B), Cvirb_B_(Cvirb_B), eps_occb_B_(eps_occb_B), eps_virb_B_(eps_virb_B) {}

CPKSUSAPTHamiltonian::~CPKSUSAPTHamiltonian() {}

void CPKSUSAPTHamiltonian::print_header() const {
    if (print_) {
        outfile->Printf("  ==> CPKSUSAPTHamiltonian <== \n\n");
        outfile->Printf("  Unrestricted two-monomer CPKS solver for USAPT0\n");
        outfile->Printf("  Product: (2J_total - K - K^T + diagonal)\n\n");
    }
}

std::shared_ptr<Vector> CPKSUSAPTHamiltonian::diagonal() {
    // Create combined diagonal for all monomer/spin combinations
    int noa_A = eps_occa_A_->dim();
    int nva_A = eps_vira_A_->dim();
    int nob_A = eps_occb_A_->dim();
    int nvb_A = eps_virb_A_->dim();
    int noa_B = eps_occa_B_->dim();
    int nva_B = eps_vira_B_->dim();
    int nob_B = eps_occb_B_->dim();
    int nvb_B = eps_virb_B_->dim();

    int total_dim = noa_A * nva_A + nob_A * nvb_A + noa_B * nva_B + nob_B * nvb_B;

    Dimension dim(1);
    dim[0] = total_dim;
    auto diag = std::make_shared<Vector>("USAPT CPKS Diagonal", dim);
    double* dp = diag->pointer();

    int offset = 0;

    // Monomer A alpha
    double* eoa_A = eps_occa_A_->pointer();
    double* eva_A = eps_vira_A_->pointer();
    for (int i = 0; i < noa_A; i++) {
        for (int a = 0; a < nva_A; a++) {
            dp[offset++] = eva_A[a] - eoa_A[i];
        }
    }

    // Monomer A beta
    double* eob_A = eps_occb_A_->pointer();
    double* evb_A = eps_virb_A_->pointer();
    for (int i = 0; i < nob_A; i++) {
        for (int a = 0; a < nvb_A; a++) {
            dp[offset++] = evb_A[a] - eob_A[i];
        }
    }

    // Monomer B alpha
    double* eoa_B = eps_occa_B_->pointer();
    double* eva_B = eps_vira_B_->pointer();
    for (int i = 0; i < noa_B; i++) {
        for (int a = 0; a < nva_B; a++) {
            dp[offset++] = eva_B[a] - eoa_B[i];
        }
    }

    // Monomer B beta
    double* eob_B = eps_occb_B_->pointer();
    double* evb_B = eps_virb_B_->pointer();
    for (int i = 0; i < nob_B; i++) {
        for (int a = 0; a < nvb_B; a++) {
            dp[offset++] = evb_B[a] - eob_B[i];
        }
    }

    return diag;
}

void CPKSUSAPTHamiltonian::product(const std::vector<std::shared_ptr<Vector> >& x,
                                   std::vector<std::shared_ptr<Vector> >& b) {
    throw PSIEXCEPTION("CPKSUSAPTHamiltonian: Use product_map interface for USAPT calculations");
}

std::map<std::string, SharedMatrix> CPKSUSAPTHamiltonian::product_map(
    std::map<std::string, SharedMatrix>& b) {

    std::map<std::string, SharedMatrix> s;

    bool do_A = b.count("Aa") || b.count("Ab");
    bool do_B = b.count("Ba") || b.count("Bb");

    if (!do_A && !do_B) {
        return s;  // Empty result
    }

    // Determine which spin components are active and non-empty
    bool alpha_A = false, beta_A = false, alpha_B = false, beta_B = false;
    if (do_A) {
        if (b.count("Aa")) {
            alpha_A = (b["Aa"]->nrow() > 0) && (b["Aa"]->ncol() > 0);
        }
        if (b.count("Ab")) {
            beta_A = (b["Ab"]->nrow() > 0) && (b["Ab"]->ncol() > 0);
        }
    }
    if (do_B) {
        if (b.count("Ba")) {
            alpha_B = (b["Ba"]->nrow() > 0) && (b["Ba"]->ncol() > 0);
        }
        if (b.count("Bb")) {
            beta_B = (b["Bb"]->nrow() > 0) && (b["Bb"]->ncol() > 0);
        }
    }

    // Set up JK calculation
    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    // Prepare density matrices for JK computation
    // Order: Aa, Ab, Ba, Bb (only if active)
    if (do_A && alpha_A) {
        Cl.push_back(Cocca_A_);
        size_t no = b["Aa"]->nrow();
        size_t nv = b["Aa"]->ncol();
        size_t nso = Cvira_A_->nrow();
        double** Cp = Cvira_A_->pointer();
        double** bp = b["Aa"]->pointer();
        auto T = std::make_shared<Matrix>("T", nso, no);
        double** Tp = T->pointer();
        C_DGEMM('N', 'T', nso, no, nv, 1.0, Cp[0], nv, bp[0], nv, 0.0, Tp[0], no);
        Cr.push_back(T);
    }

    if (do_A && beta_A) {
        Cl.push_back(Coccb_A_);
        size_t no = b["Ab"]->nrow();
        size_t nv = b["Ab"]->ncol();
        size_t nso = Cvirb_A_->nrow();
        double** Cp = Cvirb_A_->pointer();
        double** bp = b["Ab"]->pointer();
        auto T = std::make_shared<Matrix>("T", nso, no);
        double** Tp = T->pointer();
        C_DGEMM('N', 'T', nso, no, nv, 1.0, Cp[0], nv, bp[0], nv, 0.0, Tp[0], no);
        Cr.push_back(T);
    }

    if (do_B && alpha_B) {
        Cl.push_back(Cocca_B_);
        size_t no = b["Ba"]->nrow();
        size_t nv = b["Ba"]->ncol();
        size_t nso = Cvira_B_->nrow();
        double** Cp = Cvira_B_->pointer();
        double** bp = b["Ba"]->pointer();
        auto T = std::make_shared<Matrix>("T", nso, no);
        double** Tp = T->pointer();
        C_DGEMM('N', 'T', nso, no, nv, 1.0, Cp[0], nv, bp[0], nv, 0.0, Tp[0], no);
        Cr.push_back(T);
    }

    if (do_B && beta_B) {
        Cl.push_back(Coccb_B_);
        size_t no = b["Bb"]->nrow();
        size_t nv = b["Bb"]->ncol();
        size_t nso = Cvirb_B_->nrow();
        double** Cp = Cvirb_B_->pointer();
        double** bp = b["Bb"]->pointer();
        auto T = std::make_shared<Matrix>("T", nso, no);
        double** Tp = T->pointer();
        C_DGEMM('N', 'T', nso, no, nv, 1.0, Cp[0], nv, bp[0], nv, 0.0, Tp[0], no);
        Cr.push_back(T);
    }

    // Compute J and K matrices
    jk_->compute();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    // Track indices into J/K arrays
    int ind_jk = 0;

    // Process Monomer A
    if (do_A) {
        SharedMatrix Jva, Jvb, Kva, Kvb, T;

        // Get J and K for each active spin
        if (alpha_A) {
            Jva = J[ind_jk];
            Kva = K[ind_jk];
            ++ind_jk;
            Jva->scale(2.0);
            T = Jva->clone();
        }

        if (beta_A) {
            Jvb = J[ind_jk];
            Kvb = K[ind_jk];
            ++ind_jk;
            Jvb->scale(2.0);
            if (alpha_A) {
                T->add(Jvb);
            } else {
                T = Jvb->clone();
            }
        }

        // Form products for alpha spin
        if (alpha_A) {
            Jva->copy(T);  // J_total = 2(J_alpha + J_beta)
            Jva->subtract(Kva);
            Jva->subtract(Kva->transpose());

            int no = b["Aa"]->nrow();
            int nv = b["Aa"]->ncol();
            s["Aa"] = std::make_shared<Matrix>("SAa", no, nv);

            int nso = Cvira_A_->nrow();
            auto Temp = std::make_shared<Matrix>("T", no, nso);
            double** Cop = Cocca_A_->pointer();
            double** Cvp = Cvira_A_->pointer();
            double** Jp = Jva->pointer();
            double** Tp = Temp->pointer();
            double** Sp = s["Aa"]->pointer();

            C_DGEMM('T', 'N', no, nso, nso, 1.0, Cop[0], no, Jp[0], nso, 0.0, Tp[0], nso);
            C_DGEMM('N', 'N', no, nv, nso, 1.0, Tp[0], nso, Cvp[0], nv, 0.0, Sp[0], nv);

            // Add diagonal term
            double** bp = b["Aa"]->pointer();
            double* eop = eps_occa_A_->pointer();
            double* evp = eps_vira_A_->pointer();
            for (int i = 0; i < no; i++) {
                for (int a = 0; a < nv; a++) {
                    Sp[i][a] += bp[i][a] * (evp[a] - eop[i]);
                }
            }
        }

        // Form products for beta spin
        if (beta_A) {
            Jvb->copy(T);  // J_total = 2(J_alpha + J_beta)
            Jvb->subtract(Kvb);
            Jvb->subtract(Kvb->transpose());

            int no = b["Ab"]->nrow();
            int nv = b["Ab"]->ncol();
            s["Ab"] = std::make_shared<Matrix>("SAb", no, nv);

            int nso = Cvirb_A_->nrow();
            auto Temp = std::make_shared<Matrix>("T", no, nso);
            double** Cop = Coccb_A_->pointer();
            double** Cvp = Cvirb_A_->pointer();
            double** Jp = Jvb->pointer();
            double** Tp = Temp->pointer();
            double** Sp = s["Ab"]->pointer();

            C_DGEMM('T', 'N', no, nso, nso, 1.0, Cop[0], no, Jp[0], nso, 0.0, Tp[0], nso);
            C_DGEMM('N', 'N', no, nv, nso, 1.0, Tp[0], nso, Cvp[0], nv, 0.0, Sp[0], nv);

            // Add diagonal term
            double** bp = b["Ab"]->pointer();
            double* eop = eps_occb_A_->pointer();
            double* evp = eps_virb_A_->pointer();
            for (int i = 0; i < no; i++) {
                for (int a = 0; a < nv; a++) {
                    Sp[i][a] += bp[i][a] * (evp[a] - eop[i]);
                }
            }
        }
    }

    // Process Monomer B (similar logic)
    if (do_B) {
        SharedMatrix Jva, Jvb, Kva, Kvb, T;

        // Get J and K for each active spin
        if (alpha_B) {
            Jva = J[ind_jk];
            Kva = K[ind_jk];
            ++ind_jk;
            Jva->scale(2.0);
            T = Jva->clone();
        }

        if (beta_B) {
            Jvb = J[ind_jk];
            Kvb = K[ind_jk];
            ++ind_jk;
            Jvb->scale(2.0);
            if (alpha_B) {
                T->add(Jvb);
            } else {
                T = Jvb->clone();
            }
        }

        // Form products for alpha spin
        if (alpha_B) {
            Jva->copy(T);  // J_total = 2(J_alpha + J_beta)
            Jva->subtract(Kva);
            Jva->subtract(Kva->transpose());

            int no = b["Ba"]->nrow();
            int nv = b["Ba"]->ncol();
            s["Ba"] = std::make_shared<Matrix>("SBa", no, nv);

            int nso = Cvira_B_->nrow();
            auto Temp = std::make_shared<Matrix>("T", no, nso);
            double** Cop = Cocca_B_->pointer();
            double** Cvp = Cvira_B_->pointer();
            double** Jp = Jva->pointer();
            double** Tp = Temp->pointer();
            double** Sp = s["Ba"]->pointer();

            C_DGEMM('T', 'N', no, nso, nso, 1.0, Cop[0], no, Jp[0], nso, 0.0, Tp[0], nso);
            C_DGEMM('N', 'N', no, nv, nso, 1.0, Tp[0], nso, Cvp[0], nv, 0.0, Sp[0], nv);

            // Add diagonal term
            double** bp = b["Ba"]->pointer();
            double* eop = eps_occa_B_->pointer();
            double* evp = eps_vira_B_->pointer();
            for (int i = 0; i < no; i++) {
                for (int a = 0; a < nv; a++) {
                    Sp[i][a] += bp[i][a] * (evp[a] - eop[i]);
                }
            }
        }

        // Form products for beta spin
        if (beta_B) {
            Jvb->copy(T);  // J_total = 2(J_alpha + J_beta)
            Jvb->subtract(Kvb);
            Jvb->subtract(Kvb->transpose());

            int no = b["Bb"]->nrow();
            int nv = b["Bb"]->ncol();
            s["Bb"] = std::make_shared<Matrix>("SBb", no, nv);

            int nso = Cvirb_B_->nrow();
            auto Temp = std::make_shared<Matrix>("T", no, nso);
            double** Cop = Coccb_B_->pointer();
            double** Cvp = Cvirb_B_->pointer();
            double** Jp = Jvb->pointer();
            double** Tp = Temp->pointer();
            double** Sp = s["Bb"]->pointer();

            C_DGEMM('T', 'N', no, nso, nso, 1.0, Cop[0], no, Jp[0], nso, 0.0, Tp[0], nso);
            C_DGEMM('N', 'N', no, nv, nso, 1.0, Tp[0], nso, Cvp[0], nv, 0.0, Sp[0], nv);

            // Add diagonal term
            double** bp = b["Bb"]->pointer();
            double* eop = eps_occb_B_->pointer();
            double* evp = eps_virb_B_->pointer();
            for (int i = 0; i < no; i++) {
                for (int a = 0; a < nv; a++) {
                    Sp[i][a] += bp[i][a] * (evp[a] - eop[i]);
                }
            }
        }
    }

    return s;
}

void CPKSUSAPTHamiltonian::preconditioner(SharedMatrix r, SharedMatrix z,
                                          std::shared_ptr<Vector> eps_occ,
                                          std::shared_ptr<Vector> eps_vir) {
    int no = eps_occ->dim();
    int nv = eps_vir->dim();

    double** rp = r->pointer();
    double** zp = z->pointer();
    double* op = eps_occ->pointer();
    double* vp = eps_vir->pointer();

    for (int i = 0; i < no; i++) {
        for (int a = 0; a < nv; a++) {
            zp[i][a] = rp[i][a] / (vp[a] - op[i]);
        }
    }
}

}  // namespace psi
