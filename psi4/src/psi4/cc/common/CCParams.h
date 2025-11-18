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

#ifndef CC_COMMON_CCPARAMS_H
#define CC_COMMON_CCPARAMS_H

/*! \file
    \ingroup CC
    \brief Common parameters shared across CC modules
*/

#include <string>

namespace psi {
namespace cc {
namespace common {

/*!
 * \brief Base structure containing parameters common to multiple CC modules
 *
 * This structure contains parameters that are used by multiple CC modules
 * (ccenergy, cclambda, cceom, ccdensity, cchbar, ccresponse, cctriples).
 * Individual modules extend this structure with their module-specific parameters.
 *
 * Universal parameters (present in all 7 modules):
 * - ref: Reference wavefunction type
 * - wfn: Wavefunction/method type
 *
 * Common parameters (present in 6+ modules):
 * - memory: Available memory
 * - cachelev: DPD cache level
 * - dertype: Derivative type
 *
 * Frequently used parameters (present in 3-5 modules):
 * - restart, local, abcd, print, maxiter, convergence, diis, num_amps
 */
struct CCParams {
    // =========================================================================
    // Universal Parameters (ALL 7 modules)
    // =========================================================================

    /*! Reference wavefunction type:
     *  0 = RHF (restricted Hartree-Fock)
     *  1 = ROHF (restricted open-shell Hartree-Fock)
     *  2 = UHF (unrestricted Hartree-Fock)
     */
    int ref;

    /*! Wavefunction/method type (e.g., "CCSD", "CCSD_T", "EOM_CCSD", etc.) */
    std::string wfn;

    // =========================================================================
    // Common Parameters (6/7 modules - missing in cctriples only)
    // =========================================================================

    /*! Available memory in bytes */
    long int memory;

    /*! DPD cache level for disk I/O optimization */
    int cachelev;

    // =========================================================================
    // Common Parameters (6/7 modules - missing in cceom only)
    // =========================================================================

    /*! Derivative type:
     *  0 = NONE (energy only)
     *  1 = FIRST (first derivative/gradient)
     *  3 = RESPONSE (linear response)
     */
    int dertype;

    // =========================================================================
    // Frequently Used Parameters (4/7 modules)
    // =========================================================================

    /*! Restart from previous amplitudes on disk */
    int restart;

    /*! Use local correlation methods */
    int local;

    /*! Algorithm for (AB|CD) integrals (e.g., "NEW", "OLD") */
    std::string abcd;

    /*! Output print level */
    int print;

    // =========================================================================
    // Moderately Used Parameters (3/7 modules)
    // =========================================================================

    /*! Maximum number of iterations */
    int maxiter;

    /*! Convergence criterion for amplitudes */
    double convergence;

    /*! Use DIIS (Direct Inversion in the Iterative Subspace) */
    int diis;

    /*! Number of largest amplitudes to print */
    int num_amps;

    // =========================================================================
    // Constructor with default values
    // =========================================================================

    CCParams()
        : ref(0),
          wfn("NONE"),
          memory(0),
          cachelev(2),
          dertype(0),
          restart(0),
          local(0),
          abcd("NEW"),
          print(0),
          maxiter(50),
          convergence(1e-7),
          diis(1),
          num_amps(10) {}
};

}  // namespace common
}  // namespace cc
}  // namespace psi

#endif  // CC_COMMON_CCPARAMS_H
