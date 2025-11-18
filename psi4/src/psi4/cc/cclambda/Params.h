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

#ifndef CCLAMBDA_PARAMS_H
#define CCLAMBDA_PARAMS_H

/*! \file
    \ingroup CCLAMBDA
    \brief Parameters for cclambda module
*/

#include <string>
#include "psi4/cc/common/CCParams.h"

namespace psi {
namespace cclambda {

/*! \brief Parameters for cclambda module
 *
 * Extends common CC parameters with cclambda-specific parameters.
 * Inherits from CCParams (via ccenergy indirectly).
 */
struct Params : public psi::cc::common::CCParams {
    // Module-specific parameters
    int aobasis;    // Use AO basis (boolean, int type)
    int nstates;    // Total number of L vectors to compute
    int zeta;       // Boolean for solving zeta equations (implies excited state)
    int sekino;     // Sekino-Bartlett size-extensive models
    int all;        // Find Ls for all excited states plus ground state (obsolete)
    int ground;     // Find L for only ground state (obsolete)

    Params() : CCParams(), aobasis(0), nstates(0), zeta(0), sekino(0), all(0), ground(0) {}
};

struct L_Params {
    int irrep;           /* same as corresponding R */
    double R0;           /* same as corresponding R */
    double cceom_energy; /* same as corresponding R */
    int root;            /* index of root within irrep */
    bool ground;         /* boolean, is this a ground state L ? */
    char L1A_lbl[32];
    char L1B_lbl[32];
    char L2AA_lbl[32];
    char L2BB_lbl[32];
    char L2AB_lbl[32];
    char L2RHF_lbl[32];
};

}  // namespace cclambda
}  // namespace psi

#endif
