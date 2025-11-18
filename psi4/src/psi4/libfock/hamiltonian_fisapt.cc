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
 * @file hamiltonian_fisapt.cc
 * @brief CPHF Hamiltonian for F-SAPT two-monomer calculations
 *
 * This file implements the CPHFFISAPTHamiltonian class, which provides
 * Hamiltonian-vector products for F-SAPT coupled response calculations.
 * The implementation consolidates product logic previously duplicated in
 * CPHF_FISAPT::product() (fisapt/fisapt.cc).
 *
 * Key features:
 * - Two-monomer problem structure (monomers A and B)
 * - Map-based interface for flexible perturbation handling
 * - Product formula: (4J - K - K^T + diagonal)
 * - Static preconditioner helper for Jacobi preconditioning
 *
 * This class is used internally by CPHF_FISAPT for backwards compatibility
 * while providing a clean, reusable implementation of the product logic.
 *
 * @see libfock/CPHF_ARCHITECTURE.md for architecture documentation
 * @see hamiltonian_fisapt.h for class interface
 */

#include "hamiltonian_fisapt.h"
#include "jk.h"

#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {

CPHFFISAPTHamiltonian::CPHFFISAPTHamiltonian(std::shared_ptr<JK> jk,
                                             SharedMatrix Cocc_A, SharedMatrix Cvir_A,
                                             std::shared_ptr<Vector> eps_occ_A, std::shared_ptr<Vector> eps_vir_A,
                                             SharedMatrix Cocc_B, SharedMatrix Cvir_B,
                                             std::shared_ptr<Vector> eps_occ_B, std::shared_ptr<Vector> eps_vir_B)
    : RHamiltonian(jk),
      Cocc_A_(Cocc_A), Cvir_A_(Cvir_A), eps_occ_A_(eps_occ_A), eps_vir_A_(eps_vir_A),
      Cocc_B_(Cocc_B), Cvir_B_(Cvir_B), eps_occ_B_(eps_occ_B), eps_vir_B_(eps_vir_B) {}

CPHFFISAPTHamiltonian::~CPHFFISAPTHamiltonian() {}

void CPHFFISAPTHamiltonian::print_header() const {
    if (print_) {
        outfile->Printf("  ==> CPHFFISAPTHamiltonian <== \n\n");
        outfile->Printf("  Two-monomer CPHF solver for F-SAPT\n");
        outfile->Printf("  Product: (4J - K - K^T + diagonal)\n\n");
    }
}

std::shared_ptr<Vector> CPHFFISAPTHamiltonian::diagonal() {
    // For FISAPT, create combined diagonal for both monomers
    // This is primarily used for initial guess in CGRSolver
    int no_A = eps_occ_A_->dim();
    int nv_A = eps_vir_A_->dim();
    int no_B = eps_occ_B_->dim();
    int nv_B = eps_vir_B_->dim();

    int total_dim = no_A * nv_A + no_B * nv_B;

    Dimension dim(1);
    dim[0] = total_dim;
    auto diag = std::make_shared<Vector>("FISAPT CPHF Diagonal", dim);
    double* dp = diag->pointer();

    // Monomer A diagonal: eps_vir - eps_occ
    double* eop_A = eps_occ_A_->pointer();
    double* evp_A = eps_vir_A_->pointer();
    int offset = 0;
    for (int i = 0; i < no_A; i++) {
        for (int a = 0; a < nv_A; a++) {
            dp[offset++] = evp_A[a] - eop_A[i];
        }
    }

    // Monomer B diagonal: eps_vir - eps_occ
    double* eop_B = eps_occ_B_->pointer();
    double* evp_B = eps_vir_B_->pointer();
    for (int i = 0; i < no_B; i++) {
        for (int a = 0; a < nv_B; a++) {
            dp[offset++] = evp_B[a] - eop_B[i];
        }
    }

    return diag;
}

void CPHFFISAPTHamiltonian::product(const std::vector<std::shared_ptr<Vector> >& x,
                                    std::vector<std::shared_ptr<Vector> >& b) {
    throw PSIEXCEPTION("CPHFFISAPTHamiltonian: Use product_map interface for FISAPT calculations");
}

std::map<std::string, SharedMatrix> CPHFFISAPTHamiltonian::product_map(
    const std::map<std::string, SharedMatrix>& b) {

    std::map<std::string, SharedMatrix> s;

    bool do_A = b.count("A");
    bool do_B = b.count("B");

    if (!do_A && !do_B) {
        return s;  // Empty result
    }

    // Set up JK calculation
    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    // Prepare density matrices for JK computation
    if (do_A) {
        Cl.push_back(Cocc_A_);
        auto T = linalg::doublet(Cvir_A_, b.at("A"), false, true);
        Cr.push_back(T);
    }

    if (do_B) {
        Cl.push_back(Cocc_B_);
        auto T = linalg::doublet(Cvir_B_, b.at("B"), false, true);
        Cr.push_back(T);
    }

    // Compute J and K matrices
    jk_->compute();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    int indA = 0;
    int indB = (do_A ? 1 : 0);

    // Process monomer A
    if (do_A) {
        SharedMatrix Jv = J[indA];
        SharedMatrix Kv = K[indA];

        // Form 4J - K - K^T
        Jv->scale(4.0);
        Jv->subtract(Kv);
        Jv->subtract(Kv->transpose());

        int no = b.at("A")->nrow();
        int nv = b.at("A")->ncol();

        // s = C_occ^T * (4J - K - K^T) * C_vir
        s["A"] = linalg::triplet(Cocc_A_, Jv, Cvir_A_, true, false, false);

        // Add diagonal term: (eps_vir - eps_occ) * b
        double** Sp = s["A"]->pointer();
        double** bp = b.at("A")->pointer();
        double* op = eps_occ_A_->pointer();
        double* vp = eps_vir_A_->pointer();

        for (int i = 0; i < no; i++) {
            for (int a = 0; a < nv; a++) {
                Sp[i][a] += bp[i][a] * (vp[a] - op[i]);
            }
        }
    }

    // Process monomer B
    if (do_B) {
        SharedMatrix Jv = J[indB];
        SharedMatrix Kv = K[indB];

        // Form 4J - K - K^T
        Jv->scale(4.0);
        Jv->subtract(Kv);
        Jv->subtract(Kv->transpose());

        int no = b.at("B")->nrow();
        int nv = b.at("B")->ncol();

        // s = C_occ^T * (4J - K - K^T) * C_vir
        s["B"] = linalg::triplet(Cocc_B_, Jv, Cvir_B_, true, false, false);

        // Add diagonal term: (eps_vir - eps_occ) * b
        double** Sp = s["B"]->pointer();
        double** bp = b.at("B")->pointer();
        double* op = eps_occ_B_->pointer();
        double* vp = eps_vir_B_->pointer();

        for (int i = 0; i < no; i++) {
            for (int a = 0; a < nv; a++) {
                Sp[i][a] += bp[i][a] * (vp[a] - op[i]);
            }
        }
    }

    return s;
}

void CPHFFISAPTHamiltonian::preconditioner(SharedMatrix r, SharedMatrix z,
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

std::map<std::string, SharedVector> CPHFFISAPTHamiltonian::pack(
    const std::map<std::string, SharedMatrix>& mat_map) {

    std::map<std::string, SharedVector> vec_map;

    for (const auto& pair : mat_map) {
        const std::string& key = pair.first;
        SharedMatrix mat = pair.second;

        int no = mat->nrow();
        int nv = mat->ncol();

        Dimension dim(1);
        dim[0] = no * nv;
        auto vec = std::make_shared<Vector>(key, dim);

        double** mp = mat->pointer();
        double* vp = vec->pointer();

        int offset = 0;
        for (int i = 0; i < no; i++) {
            for (int a = 0; a < nv; a++) {
                vp[offset++] = mp[i][a];
            }
        }

        vec_map[key] = vec;
    }

    return vec_map;
}

std::map<std::string, SharedMatrix> CPHFFISAPTHamiltonian::unpack(
    const std::vector<SharedVector>& vec_list) {

    std::map<std::string, SharedMatrix> mat_map;

    // For FISAPT, we expect vectors named "A" and/or "B"
    for (size_t i = 0; i < vec_list.size(); i++) {
        SharedVector vec = vec_list[i];
        std::string key = vec->name();

        // Determine dimensions based on key
        int no, nv;
        if (key == "A" || key.find("A") != std::string::npos) {
            no = eps_occ_A_->dim();
            nv = eps_vir_A_->dim();
        } else if (key == "B" || key.find("B") != std::string::npos) {
            no = eps_occ_B_->dim();
            nv = eps_vir_B_->dim();
        } else {
            throw PSIEXCEPTION("CPHFFISAPTHamiltonian::unpack: Unknown key " + key);
        }

        auto mat = std::make_shared<Matrix>(key, no, nv);
        double** mp = mat->pointer();
        double* vp = vec->pointer();

        int offset = 0;
        for (int i = 0; i < no; i++) {
            for (int a = 0; a < nv; a++) {
                mp[i][a] = vp[offset++];
            }
        }

        mat_map[key] = mat;
    }

    return mat_map;
}

}  // namespace psi
