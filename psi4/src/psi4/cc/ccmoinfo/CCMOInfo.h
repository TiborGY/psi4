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

#ifndef CCMOINFO_H
#define CCMOINFO_H

/*! \file
    \ingroup CC
    \brief Unified MOInfo structure for all CC modules

    This class consolidates the duplicate MOInfo structures previously
    defined separately in ccenergy, cclambda, cceom, ccdensity, ccresponse,
    cchbar, and cctriples modules.

    Design goals:
    - Eliminate ~700 lines of duplicate code
    - Use modern C++ (vectors, smart pointers, RAII)
    - Support common fields needed by all modules
    - Support optional module-specific fields
    - Maintain backward compatibility through accessor methods
*/

#include <string>
#include <vector>
#include <memory>
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/matrix.h"

namespace psi {

class Wavefunction;

namespace ccmoinfo {

/**
 * @class CCMOInfo
 * @brief Unified molecular orbital information for coupled cluster modules
 *
 * This class provides a centralized storage for MO information used across
 * all CC modules. It replaces the duplicate MOInfo structs that were
 * previously maintained separately in each module.
 *
 * Field categories:
 * - Basic dimensions (nirreps, nmo, nso, etc.) - always present
 * - Orbital counts per irrep - always present
 * - Symmetry and offset arrays - always present
 * - Energy values - common across modules
 * - Iteration/convergence tracking - optional
 * - Diagnostics - optional (ccenergy, cceom)
 * - Transformation matrices - optional
 * - Reordering arrays - optional (ccdensity)
 * - Density matrices - optional (ccdensity)
 */
class CCMOInfo {
   public:
    /**
     * Default constructor
     */
    CCMOInfo();

    /**
     * Destructor - handles cleanup automatically via RAII
     */
    ~CCMOInfo();

    /**
     * Initialize from wavefunction and PSIO data
     * @param wfn Reference wavefunction
     * @param reference Reference type (0=RHF, 1=ROHF, 2=UHF)
     */
    void initialize(std::shared_ptr<Wavefunction> wfn, int reference);

    //
    // ========== COMMON FIELDS (Present in all or most modules) ==========
    //

    // Basic dimensions
    int nirreps;                     ///< Number of irreducible representations
    int nmo;                         ///< Number of molecular orbitals
    int nso;                         ///< Number of symmetry orbitals
    int nao;                         ///< Number of atomic orbitals

    Dimension sopi;                  ///< Number of SOs per irrep
    Dimension orbspi;                ///< Number of MOs per irrep
    Dimension clsdpi;                ///< Number of closed-shell orbitals per irrep (excl. frozen)
    Dimension openpi;                ///< Number of open-shell orbitals per irrep
    Dimension uoccpi;                ///< Number of unoccupied orbitals per irrep (excl. frozen virt)
    Dimension frdocc;                ///< Number of frozen doubly-occupied orbitals per irrep
    Dimension fruocc;                ///< Number of frozen virtual orbitals per irrep

    int nvirt;                       ///< Total number of (active) virtual orbitals

    std::vector<std::string> labels; ///< Irrep labels

    // Orbital counts per irrep (for active space)
    Dimension occpi;                 ///< Number of occupied orbitals (incl. open) per irrep
    Dimension aoccpi;                ///< Number of alpha occupied orbitals (incl. open) per irrep
    Dimension boccpi;                ///< Number of beta occupied orbitals (incl. open) per irrep
    Dimension virtpi;                ///< Number of virtual orbitals (incl. open) per irrep
    Dimension avirtpi;               ///< Number of alpha virtual orbitals (incl. open) per irrep
    Dimension bvirtpi;               ///< Number of beta virtual orbitals (incl. open) per irrep

    // Symmetry arrays (use vectors for automatic memory management)
    std::vector<int> occ_sym;        ///< Relative occupied index symmetry
    std::vector<int> aocc_sym;       ///< Relative alpha occupied index symmetry
    std::vector<int> bocc_sym;       ///< Relative beta occupied index symmetry
    std::vector<int> vir_sym;        ///< Relative virtual index symmetry
    std::vector<int> avir_sym;       ///< Relative alpha virtual index symmetry
    std::vector<int> bvir_sym;       ///< Relative beta virtual index symmetry
    std::vector<int> sosym;          ///< SO symmetry (Pitzer ordering)

    // Offset arrays (use vectors for automatic memory management)
    std::vector<int> occ_off;        ///< Occupied orbital offsets within each irrep
    std::vector<int> aocc_off;       ///< Alpha occupied orbital offsets within each irrep
    std::vector<int> bocc_off;       ///< Beta occupied orbital offsets within each irrep
    std::vector<int> vir_off;        ///< Virtual orbital offsets within each irrep
    std::vector<int> avir_off;       ///< Alpha virtual orbital offsets within each irrep
    std::vector<int> bvir_off;       ///< Beta virtual orbital offsets within each irrep

    // Common energy values
    double enuc;                     ///< Nuclear repulsion energy
    double escf;                     ///< SCF energy (from wavefunction)
    double eref;                     ///< Reference energy (from PSIO file)
    double ecc;                      ///< Coupled cluster energy

    //
    // ========== OPTIONAL FIELDS (Module-specific) ==========
    //

    // Iteration and convergence tracking
    int iter;                        ///< Current iteration number
    double conv;                     ///< Current convergence level
    int sym;                         ///< Symmetry of converged CCSD state

    // MP2 energies (ccenergy)
    double emp2;                     ///< MP2 energy
    double emp2_ss;                  ///< Same-spin MP2 correlation energy
    double emp2_os;                  ///< Opposite-spin MP2 correlation energy
    double emp2_s;                   ///< Singles MP2 correlation energy

    // CC energy components (ccenergy)
    double ecc_ss;                   ///< Same-spin CC correlation energy
    double ecc_os;                   ///< Opposite-spin CC correlation energy
    double ecc_s;                    ///< Singles CC correlation energy

    // Diagnostics (ccenergy, cceom)
    double t1diag;                   ///< Standard T1 diagnostic
    double d1diag;                   ///< Janssen-Nielsen D1 diagnostic
    double new_d1diag;               ///< Lee's modified D1 diagnostic
    double d2diag;                   ///< Nielsen-Janssen D2 diagnostic

    // Lambda pseudoenergy (cclambda)
    double lcc;                      ///< Lambda pseudoenergy

    // Triples energy (ccdensity)
    double et;                       ///< (T) energy from cctriples

    // Open-shell flag (cceom)
    int iopen;                       ///< 0=closed shell, >0=open shell

    // Active orbital counts (ccdensity, ccresponse)
    int nactive;                     ///< Number of active orbitals
    int nfzc;                        ///< Total number of frozen core orbitals
    int nfzv;                        ///< Total number of frozen virtual orbitals
    int nclsd;                       ///< Total number of closed shells (excl. frozen)
    int nopen;                       ///< Total number of open shells
    int nuocc;                       ///< Total number of unoccupied shells (excl. frozen)

    // Triangle sizes (ccresponse)
    int noei;                        ///< Number of elements in SOxSO lower triangle
    int ntri;                        ///< Number of elements in MOxMO lower triangle
    int noei_ao;                     ///< Number of elements in AOxAO lower triangle

    // Active dimensions (ccresponse)
    Dimension actpi;                 ///< Number of active orbitals per irrep
    Dimension act_occpi;             ///< Dimension form of occpi

    // Lowercase irrep labels (cceom)
    std::vector<std::string> irr_labs_lowercase;

    //
    // ========== ORBITAL TRANSFORMATION MATRICES ==========
    //
    // These are stored as raw pointers for compatibility with existing code
    // They are allocated and deallocated by helper methods
    //

    // ccenergy naming convention
    double ***Cv;                    ///< Virtual orbital transformation matrix
    double ***Cav;                   ///< Alpha virtual orbital transformation matrix
    double ***Cbv;                   ///< Beta virtual orbital transformation matrix
    double ***Co;                    ///< Occupied orbital transformation matrix
    double ***Cao;                   ///< Alpha occupied orbital transformation matrix
    double ***Cbo;                   ///< Beta occupied orbital transformation matrix

    // cclambda/cceom/ccresponse naming convention
    double ***C;                     ///< Virtual orbital transformation matrix
    double ***Ca;                    ///< Alpha virtual orbital transformation matrix
    double ***Cb;                    ///< Beta virtual orbital transformation matrix

    // ccresponse Matrix version (for active MOs)
    std::shared_ptr<Matrix> Ca_matrix;  ///< Active MO coefficients (ccresponse)

    //
    // ========== POLARIZABILITY/RESPONSE-SPECIFIC (ccresponse) ==========
    //

    int *mu_irreps;                  ///< Irreps of x,y,z dipole components
    int *l_irreps;                   ///< Irreps of x,y,z angular momentum components

    //
    // ========== REORDERING ARRAYS (ccdensity) ==========
    //

    // QT <-> CC reordering for occupied orbitals
    std::vector<int> cc_occ;         ///< QT->CC active occupied reordering
    std::vector<int> cc_aocc;        ///< QT->CC alpha active occupied reordering
    std::vector<int> cc_bocc;        ///< QT->CC beta active occupied reordering
    std::vector<int> qt_occ;         ///< CC->QT active occupied reordering
    std::vector<int> qt_aocc;        ///< CC->QT alpha active occupied reordering
    std::vector<int> qt_bocc;        ///< CC->QT beta active occupied reordering

    // QT <-> CC reordering for virtual orbitals
    std::vector<int> cc_vir;         ///< QT->CC active virtual reordering
    std::vector<int> cc_avir;        ///< QT->CC alpha active virtual reordering
    std::vector<int> cc_bvir;        ///< QT->CC beta active virtual reordering
    std::vector<int> qt_vir;         ///< CC->QT active virtual reordering
    std::vector<int> qt_avir;        ///< CC->QT alpha active virtual reordering
    std::vector<int> qt_bvir;        ///< CC->QT beta active virtual reordering

    // Pitzer <-> QT reordering
    std::vector<int> pitzer2qt;      ///< Pitzer to QT reordering
    std::vector<int> qt2pitzer;      ///< QT to Pitzer reordering

    //
    // ========== DENSITY MATRICES (ccdensity) ==========
    //

    Matrix opdm;                     ///< One-particle density matrix (full space)
    Matrix opdm_a;                   ///< Alpha one-particle density matrix (full space)
    Matrix opdm_b;                   ///< Beta one-particle density matrix (full space)

    Matrix ltd_mat;                  ///< Left transition density <0|O|n>
    Matrix ltd_a_mat;                ///< Alpha left transition density
    Matrix ltd_b_mat;                ///< Beta left transition density

    Matrix rtd_mat;                  ///< Right transition density <n|O|0>
    Matrix rtd_a_mat;                ///< Alpha right transition density
    Matrix rtd_b_mat;                ///< Beta right transition density

    //
    // ========== SCF ORBITALS AND PROPERTY MATRICES ==========
    //

    SharedMatrix Ca_matrix;          ///< SCF orbitals (standard ordering)
    std::vector<SharedMatrix> L;     ///< Angular momentum matrices
    std::vector<SharedMatrix> nabla; ///< Nabla matrices
    std::vector<SharedMatrix> dip;   ///< Dipole matrices

    //
    // ========== PROPERTY-RELATED FIELDS (ccresponse) ==========
    //

    std::vector<int> mu_irreps;      ///< Irreps of x,y,z dipole components
    std::vector<int> l_irreps;       ///< Irreps of x,y,z angular momentum components
    int natom;                       ///< Number of atoms
    std::vector<double> zvals;       ///< Atomic Z values

    //
    // ========== BACKWARD COMPATIBILITY ACCESSORS ==========
    //
    // These provide raw pointer access for code that expects old-style pointers
    //

    /** Get raw pointer to occ_sym array */
    int* get_occ_sym() { return occ_sym.empty() ? nullptr : occ_sym.data(); }
    /** Get raw pointer to aocc_sym array */
    int* get_aocc_sym() { return aocc_sym.empty() ? nullptr : aocc_sym.data(); }
    /** Get raw pointer to bocc_sym array */
    int* get_bocc_sym() { return bocc_sym.empty() ? nullptr : bocc_sym.data(); }
    /** Get raw pointer to vir_sym array */
    int* get_vir_sym() { return vir_sym.empty() ? nullptr : vir_sym.data(); }
    /** Get raw pointer to avir_sym array */
    int* get_avir_sym() { return avir_sym.empty() ? nullptr : avir_sym.data(); }
    /** Get raw pointer to bvir_sym array */
    int* get_bvir_sym() { return bvir_sym.empty() ? nullptr : bvir_sym.data(); }
    /** Get raw pointer to sosym array */
    int* get_sosym() { return sosym.empty() ? nullptr : sosym.data(); }

    /** Get raw pointer to occ_off array */
    int* get_occ_off() { return occ_off.empty() ? nullptr : occ_off.data(); }
    /** Get raw pointer to aocc_off array */
    int* get_aocc_off() { return aocc_off.empty() ? nullptr : aocc_off.data(); }
    /** Get raw pointer to bocc_off array */
    int* get_bocc_off() { return bocc_off.empty() ? nullptr : bocc_off.data(); }
    /** Get raw pointer to vir_off array */
    int* get_vir_off() { return vir_off.empty() ? nullptr : vir_off.data(); }
    /** Get raw pointer to avir_off array */
    int* get_avir_off() { return avir_off.empty() ? nullptr : avir_off.data(); }
    /** Get raw pointer to bvir_off array */
    int* get_bvir_off() { return bvir_off.empty() ? nullptr : bvir_off.data(); }

    /** Get raw pointer to mu_irreps array */
    int* get_mu_irreps() { return mu_irreps.empty() ? nullptr : mu_irreps.data(); }
    /** Get raw pointer to l_irreps array */
    int* get_l_irreps() { return l_irreps.empty() ? nullptr : l_irreps.data(); }
    /** Get raw pointer to zvals array */
    double* get_zvals() { return zvals.empty() ? nullptr : zvals.data(); }

   private:
    /**
     * Read common data from wavefunction
     */
    void read_common_data(std::shared_ptr<Wavefunction> wfn, int reference);

    /**
     * Read orbital information from PSIO
     */
    void read_orbital_info(int reference);

    /**
     * Allocate orbital transformation matrices
     */
    void allocate_transformation_matrices(std::shared_ptr<Wavefunction> wfn, int reference);

    /**
     * Free orbital transformation matrices
     */
    void free_transformation_matrices();

    // Store reference type for cleanup
    int reference_type_;
};

}  // namespace ccmoinfo
}  // namespace psi

#endif  // CCMOINFO_H
