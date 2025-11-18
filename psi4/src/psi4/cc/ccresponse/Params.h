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

#ifndef CCRESPONSE_PARAMS_H
#define CCRESPONSE_PARAMS_H

/*! \file
    \ingroup ccresponse
    \brief Parameters for ccresponse module
*/
#include <string>
#include "psi4/cc/common/CCParams.h"

namespace psi {
namespace ccresponse {

/*! \brief Parameters for ccresponse module
 *
 * Extends common CC parameters with ccresponse-specific parameters.
 */
struct Params : public psi::cc::common::CCParams {
    // Module-specific parameters
    double *omega;         // Energy of applied field (a.u) for dynamic polarizabilities
    int nomega;            // Number of field energies desired
    std::string prop;      // User-selected property
    int analyze;           // Analysis flag
    std::string gauge;     // Choice of gauge for optical rotation
    int sekino;            // Sekino-Bartlett size-extensive model-III
    int linear;            // Bartlett size-extensive (?) linear model

    Params() : CCParams(), omega(nullptr), nomega(0), prop(""), analyze(0),
               gauge("LENGTH"), sekino(0), linear(0) {}
};

}  // namespace ccresponse
}  // namespace psi
#endif

