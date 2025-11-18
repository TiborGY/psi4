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

#ifndef HAMILTONIAN_USAPT_H
#define HAMILTONIAN_USAPT_H

#include "hamiltonian.h"
#include <map>
#include <string>

namespace psi {

/**
 * Unrestricted CPKS Hamiltonian for USAPT0 calculations.
 *
 * Handles spin-unrestricted two-monomer CPKS problems where both
 * alpha and beta spins respond independently for each monomer.
 * Uses a map-based interface with keys: "Aa", "Ab", "Ba", "Bb"
 * (monomer A/B × alpha/beta spin).
 *
 * The Hamiltonian product is: (2J_total - K_alpha - K_alpha^T + diagonal)
 * where J_total = J_alpha + J_beta for unrestricted formalism.
 */
class CPKSUSAPTHamiltonian : public RHamiltonian {
   protected:
    // Monomer A alpha spin
    SharedMatrix Cocca_A_;
    SharedMatrix Cvira_A_;
    std::shared_ptr<Vector> eps_occa_A_;
    std::shared_ptr<Vector> eps_vira_A_;

    // Monomer A beta spin
    SharedMatrix Coccb_A_;
    SharedMatrix Cvirb_A_;
    std::shared_ptr<Vector> eps_occb_A_;
    std::shared_ptr<Vector> eps_virb_A_;

    // Monomer B alpha spin
    SharedMatrix Cocca_B_;
    SharedMatrix Cvira_B_;
    std::shared_ptr<Vector> eps_occa_B_;
    std::shared_ptr<Vector> eps_vira_B_;

    // Monomer B beta spin
    SharedMatrix Coccb_B_;
    SharedMatrix Cvirb_B_;
    std::shared_ptr<Vector> eps_occb_B_;
    std::shared_ptr<Vector> eps_virb_B_;

   public:
    /**
     * Constructor for USAPT CPKS Hamiltonian
     *
     * @param jk JK object for computing Coulomb and exchange matrices
     * @param Cocca_A Monomer A alpha occupied orbital coefficients
     * @param Cvira_A Monomer A alpha virtual orbital coefficients
     * @param eps_occa_A Monomer A alpha occupied orbital eigenvalues
     * @param eps_vira_A Monomer A alpha virtual orbital eigenvalues
     * @param Coccb_A Monomer A beta occupied orbital coefficients
     * @param Cvirb_A Monomer A beta virtual orbital coefficients
     * @param eps_occb_A Monomer A beta occupied orbital eigenvalues
     * @param eps_virb_A Monomer A beta virtual orbital eigenvalues
     * @param Cocca_B Monomer B alpha occupied orbital coefficients
     * @param Cvira_B Monomer B alpha virtual orbital coefficients
     * @param eps_occa_B Monomer B alpha occupied orbital eigenvalues
     * @param eps_vira_B Monomer B alpha virtual orbital eigenvalues
     * @param Coccb_B Monomer B beta occupied orbital coefficients
     * @param Cvirb_B Monomer B beta virtual orbital coefficients
     * @param eps_occb_B Monomer B beta occupied orbital eigenvalues
     * @param eps_virb_B Monomer B beta virtual orbital eigenvalues
     */
    CPKSUSAPTHamiltonian(std::shared_ptr<JK> jk,
                         SharedMatrix Cocca_A, SharedMatrix Cvira_A,
                         std::shared_ptr<Vector> eps_occa_A, std::shared_ptr<Vector> eps_vira_A,
                         SharedMatrix Coccb_A, SharedMatrix Cvirb_A,
                         std::shared_ptr<Vector> eps_occb_A, std::shared_ptr<Vector> eps_virb_A,
                         SharedMatrix Cocca_B, SharedMatrix Cvira_B,
                         std::shared_ptr<Vector> eps_occa_B, std::shared_ptr<Vector> eps_vira_B,
                         SharedMatrix Coccb_B, SharedMatrix Cvirb_B,
                         std::shared_ptr<Vector> eps_occb_B, std::shared_ptr<Vector> eps_virb_B);
    ~CPKSUSAPTHamiltonian() override;

    void print_header() const override;

    /**
     * Returns diagonal for preconditioning (orbital energy differences)
     * Combined for all spins and monomers
     */
    std::shared_ptr<Vector> diagonal() override;

    /**
     * Standard vector-based product interface (required by base class)
     * Note: For USAPT usage, prefer the map-based product_map interface
     */
    void product(const std::vector<std::shared_ptr<Vector> >& x,
                 std::vector<std::shared_ptr<Vector> >& b) override;

    /**
     * Map-based product interface for USAPT0
     *
     * Computes Hamiltonian-vector products for unrestricted response vectors.
     * Input map may contain keys "Aa", "Ab", "Ba", "Bb" for monomer/spin combinations.
     * Returns map with products for the same keys.
     *
     * Product formula: s = (2J_total - K - K^T)x + (eps_vir - eps_occ)x
     * where J_total = J_alpha + J_beta
     *
     * @param b Input map of perturbation vectors (keys: "Aa", "Ab", "Ba", "Bb")
     * @return Map of product vectors with same keys
     */
    std::map<std::string, SharedMatrix> product_map(std::map<std::string, SharedMatrix>& b);

    /**
     * Compute preconditioner z = M^-1 * r where M is diagonal approximation
     *
     * @param r Residual matrix
     * @param z Preconditioned residual (output)
     * @param eps_occ Occupied orbital eigenvalues
     * @param eps_vir Virtual orbital eigenvalues
     */
    static void preconditioner(SharedMatrix r, SharedMatrix z,
                              std::shared_ptr<Vector> eps_occ,
                              std::shared_ptr<Vector> eps_vir);
};

}  // namespace psi

#endif  // HAMILTONIAN_USAPT_H
