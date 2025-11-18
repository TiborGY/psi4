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

#ifndef CCDENSITY_PARAMS_H
#define CCDENSITY_PARAMS_H

/*! \file
    \ingroup CCDENSITY
    \brief Parameters for ccdensity module
*/

#include <string>
#include "psi4/cc/common/CCParams.h"

namespace psi {
namespace ccdensity {

/*! \brief Parameters for ccdensity module
 *
 * Extends common CC parameters with ccdensity-specific parameters.
 */
struct Params : public psi::cc::common::CCParams {
    // Module-specific parameters
    double tolerance;       // Convergence tolerance
    int aobasis;            // Use AO basis (boolean)
    bool onepdm;            // Produce ONLY the onepdm for properties
    int relax_opdm;         // Compute relaxed one-particle density matrix
    int use_zeta;           // Use zeta equations
    int calc_xi;            // Calculate xi intermediates
    int connect_xi;         // Connect xi intermediates
    int ground;             // Ground state density
    int transition;         // Transition density
    // TODO: Make meaning of variable invariant.
    //       Initially, this variable means number of states, but later
    //       it changes to number of excited states. This inconsistency
    //       should be fixed. Resolving the problem is not safe until the
    //       ccdensity module is documented.
    int nstates;            // Number of states (meaning varies in code!)
    int prop_sym;           // Property symmetry
    int prop_root;          // Property root
    int prop_all;           // Compute properties for all states
    std::string gauge;      // Gauge choice for properties
    bool write_nos;         // Write natural orbitals
    int debug_;             // Debug level

    // These are used by Xi and twopdm code
    int G_irr;              // G irrep
    int R_irr;              // R irrep
    int L_irr;              // L irrep
    double R0;              // R0 value
    double L0;              // L0 value
    double cceom_energy;    // EOM-CC energy
    double overlap1;        // <L1|R1>
    double overlap2;        // <L2|R2>
    double RD_overlap;      // Rmnef <mn||ef>
    double RZ_overlap;      // <R|zeta>

    Params() : CCParams(), tolerance(1e-7), aobasis(0), onepdm(false), relax_opdm(0),
               use_zeta(0), calc_xi(0), connect_xi(0), ground(0), transition(0),
               nstates(0), prop_sym(0), prop_root(0), prop_all(0), gauge("LENGTH"),
               write_nos(false), debug_(0), G_irr(0), R_irr(0), L_irr(0), R0(0.0),
               L0(0.0), cceom_energy(0.0), overlap1(0.0), overlap2(0.0),
               RD_overlap(0.0), RZ_overlap(0.0) {}
};

struct RHO_Params {
    int L_irr;
    int R_irr;
    int G_irr;
    int L_root;
    int R_root;
    int L_ground;
    int R_ground;
    double R0;
    double L0;
#ifdef EOM_DEBUG
    char R1A_lbl[32];
#endif
    double cceom_energy;
    double overlap1;   /* <L1|R1> */
    double overlap2;   /* <L2|R2> */
    double RD_overlap; /* Rmnef <mn||ef> */
    char DIJ_lbl[10];
    char Dij_lbl[10];
    char DAB_lbl[10];
    char Dab_lbl[10];
    char DIA_lbl[10];
    char Dia_lbl[10];
    char DAI_lbl[10];
    char Dai_lbl[10];
    char opdm_lbl[32];
    char opdm_a_lbl[32];
    char opdm_b_lbl[32];
};

// TODO: Refactor TD_Params and XTD_Params into Root_Params and Transition_Params.
// Map pairs of state identifiers to Transition_Params, and be able to grab the EOM
// from Root_Params. All ex_* files that duplicate another file can be then made obsolete.

// Describes both an excited state and the transition to it, from the ground state.
struct TD_Params {
    int irrep;           // Irrep index of the transition between states.
    int root;            // Index of the target root within irrep, excluding the ground state.
    double R0;
    double cceom_energy;
    char L1A_lbl[32];
    char L1B_lbl[32];
    char L2AA_lbl[32];
    char L2BB_lbl[32];
    char L2AB_lbl[32];
    char R1A_lbl[32];
    char R1B_lbl[32];
    char R2AA_lbl[32];
    char R2BB_lbl[32];
    char R2AB_lbl[32];
    double OS;
    double RS_length;
    double RS_velocity;
    double einstein_a;
    double einstein_b;
};

// Describes a transition between excited states.
struct XTD_Params {
    int irrep1;
    int irrep2;
    int root1;
    int root2;
    double cceom_energy;
    double OS;
    double RS_length;
    double RS_velocity;
    double einstein_a;
    double einstein_b;
};

}  // namespace ccdensity
}  // namespace psi

#endif
