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
    \ingroup CCEOM
    \brief Parameters for cceom module
*/

#ifndef _psi_src_bin_cceom_params_h
#define _psi_src_bin_cceom_params_h

#include <string>
#include "psi4/cc/common/CCParams.h"

namespace psi {
namespace cceom {

/*! \brief Parameters for cceom module
 *
 * Extends common CC parameters with cceom-specific parameters.
 */
struct Params : public psi::cc::common::CCParams {
    // Module-specific parameters
    int cachetype;       // DPD cache type
    int eom_ref;         // EOM reference type (can differ from ref for RHF references)
    int semicanonical;   // Use semicanonical orbitals
    int full_matrix;     // Include reference rows/cols in diagonalization
    int t3_Ws_incore;    // Keep T3 W intermediates in core
    int nthreads;        // Number of threads for parallel execution
    int newtrips;        // Use new triples algorithm
    int overlap;         // Check for overlaps between current wfn set and older set on disk

    Params() : CCParams(), cachetype(1), eom_ref(0), semicanonical(0), full_matrix(0),
               t3_Ws_incore(0), nthreads(1), newtrips(0), overlap(0) {}
};

struct Eom_params {
    int max_iter;
    int vectors_per_root;
    int *states_per_irrep;
    int *cs_per_irrep;
    int number_of_states;
    double eval_tol;
    double eval_tol_SS;
    double residual_tol;
    int prop_root;
    int prop_sym;
    int save_all;
    int print_singles;
    double complex_tol;
    double schmidt_add_residual_tol;
    int max_iter_SS;
    int vectors_per_root_SS;
    int excitation_range;
    double residual_tol_SS;
    std::string guess;
    int rhf_triplets;
    int mult;
    bool follow_root;
    bool collapse_with_last;
    bool collapse_with_last_cc3 ;
    int skip_diagSS;
    int vectors_cc3;
    int restart_eom_cc3;
    int amps_to_print;

    /* compute overlap of normalized R with L (must run cclambda first) */
    int dot_with_L;
    double L0;
    int L_irr;
};

}  // namespace cceom
}  // namespace psi

#endif  // _psi_src_bin_cceom_params_h
