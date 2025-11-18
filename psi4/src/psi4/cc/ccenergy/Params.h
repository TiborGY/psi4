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

/*! \file
    \ingroup CCENERGY
    \brief Parameters for ccenergy module
*/

#ifndef _psi_src_bin_ccenergy_params_h
#define _psi_src_bin_ccenergy_params_h

#include <string>
#include "psi4/cc/common/CCParams.h"

namespace psi {
namespace ccenergy {

/*! \brief Parameters for ccenergy module
 *
 * Extends common CC parameters with ccenergy-specific parameters.
 * This is the base CC module, used by cclambda and others.
 */
struct Params : public psi::cc::common::CCParams {
    // Module-specific parameters
    double e_convergence;       // Energy convergence criterion
    std::string aobasis;        // AO basis method (DISK, NONE, etc.)
    int cachetype;              // DPD cache type (LOW=1, LRU=0)
    int print_mp2_amps;         // Print MP2 amplitudes
    int brueckner;              // Use Brueckner orbitals
    double bconv;               // Brueckner orbital convergence
    int analyze;                // Analyze T2 amplitudes
    int print_pair_energies;    // Print pair energies
    int semicanonical;          // Use semicanonical orbitals
    int t2_coupled;             // Use T2-coupled equations
    std::string prop;           // User-selected property
    int just_energy;            // Just compute energy from T amplitudes on disk and quit
    int just_residuals;         // Just compute residuals from T amplitudes on disk and quit
    int t3_Ws_incore;           // Keep T3 W intermediates in core
    int nthreads;               // Number of threads for parallel execution
    int scs;                    // Use SCS-MP2 scaling
    int scsn;                   // Use SCSN-MP2 scaling
    int scscc;                  // Use SCS-CCSD scaling
    double scsmp2_scale_os;     // SCS-MP2 opposite-spin scaling factor
    double scsmp2_scale_ss;     // SCS-MP2 same-spin scaling factor
    double scscc_scale_os;      // SCS-CCSD opposite-spin scaling factor
    double scscc_scale_ss;      // SCS-CCSD same-spin scaling factor
    int newtrips;               // Use new triples algorithm
    int df;                     // Density-fitted CC

    Params() : CCParams(), e_convergence(1e-10), aobasis("NONE"), cachetype(1),
               print_mp2_amps(0), brueckner(0), bconv(1e-7), analyze(0),
               print_pair_energies(0), semicanonical(0), t2_coupled(1),
               prop(""), just_energy(0), just_residuals(0), t3_Ws_incore(0),
               nthreads(1), scs(0), scsn(0), scscc(0), scsmp2_scale_os(1.20),
               scsmp2_scale_ss(0.33), scscc_scale_os(1.27), scscc_scale_ss(1.13),
               newtrips(0), df(0) {}
};

}  // namespace ccenergy
}  // namespace psi

#endif  // _psi_src_bin_ccenergy_params_h
