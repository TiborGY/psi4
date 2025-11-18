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

#ifndef HAMILTONIAN_FISAPT_H
#define HAMILTONIAN_FISAPT_H

#include "hamiltonian.h"
#include <map>
#include <string>

namespace psi {

/**
 * CPHF Hamiltonian for F-SAPT calculations.
 *
 * Handles two-monomer CPHF problems where both monomers A and B
 * need to respond to perturbations simultaneously. Uses a map-based
 * interface to handle multiple perturbations flexibly.
 *
 * The Hamiltonian product is: (4J - K - K^T + diagonal)
 * where the J and K matrices are formed from occupied-virtual density matrices.
 */
class CPHFFISAPTHamiltonian : public RHamiltonian {
   protected:
    // Monomer A orbital coefficients and eigenvalues
    SharedMatrix Cocc_A_;
    SharedMatrix Cvir_A_;
    std::shared_ptr<Vector> eps_occ_A_;
    std::shared_ptr<Vector> eps_vir_A_;

    // Monomer B orbital coefficients and eigenvalues
    SharedMatrix Cocc_B_;
    SharedMatrix Cvir_B_;
    std::shared_ptr<Vector> eps_occ_B_;
    std::shared_ptr<Vector> eps_vir_B_;

   public:
    /**
     * Constructor for FISAPT CPHF Hamiltonian
     *
     * @param jk JK object for computing Coulomb and exchange matrices
     * @param Cocc_A Monomer A occupied orbital coefficients
     * @param Cvir_A Monomer A virtual orbital coefficients
     * @param eps_occ_A Monomer A occupied orbital eigenvalues
     * @param eps_vir_A Monomer A virtual orbital eigenvalues
     * @param Cocc_B Monomer B occupied orbital coefficients
     * @param Cvir_B Monomer B virtual orbital coefficients
     * @param eps_occ_B Monomer B occupied orbital eigenvalues
     * @param eps_vir_B Monomer B virtual orbital eigenvalues
     */
    CPHFFISAPTHamiltonian(std::shared_ptr<JK> jk,
                          SharedMatrix Cocc_A, SharedMatrix Cvir_A,
                          std::shared_ptr<Vector> eps_occ_A, std::shared_ptr<Vector> eps_vir_A,
                          SharedMatrix Cocc_B, SharedMatrix Cvir_B,
                          std::shared_ptr<Vector> eps_occ_B, std::shared_ptr<Vector> eps_vir_B);
    ~CPHFFISAPTHamiltonian() override;

    void print_header() const override;

    /**
     * Returns diagonal for preconditioning (orbital energy differences)
     * Note: For FISAPT with two monomers, this returns combined diagonal for both A and B
     */
    std::shared_ptr<Vector> diagonal() override;

    /**
     * Standard vector-based product interface (required by base class)
     * Note: For FISAPT usage, prefer the map-based product_map interface
     */
    void product(const std::vector<std::shared_ptr<Vector> >& x,
                 std::vector<std::shared_ptr<Vector> >& b) override;

    /**
     * Map-based product interface for FISAPT
     *
     * Computes Hamiltonian-vector products for monomer response vectors.
     * Input map may contain keys "A" and/or "B" for the two monomers.
     * Returns map with products for the same keys.
     *
     * Product formula: s = (4J - K - K^T)x + (eps_vir - eps_occ)x
     *
     * @param b Input map of perturbation vectors (keys: "A", "B")
     * @return Map of product vectors with same keys
     */
    std::map<std::string, SharedMatrix> product_map(const std::map<std::string, SharedMatrix>& b);

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

    /**
     * Pack matrix map into vector format (for use with CGRSolver)
     *
     * @param mat_map Map of matrices (keys: "A", "B")
     * @return Vector of packed vectors
     */
    std::map<std::string, SharedVector> pack(const std::map<std::string, SharedMatrix>& mat_map);

    /**
     * Unpack vectors into matrix map format
     *
     * @param vec_list Vector of solution vectors
     * @return Map of matrices (keys: "A", "B")
     */
    std::map<std::string, SharedMatrix> unpack(const std::vector<SharedVector>& vec_list);
};

}  // namespace psi

#endif  // HAMILTONIAN_FISAPT_H
