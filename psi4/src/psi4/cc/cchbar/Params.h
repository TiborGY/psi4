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
    \ingroup CCHBAR
    \brief Parameters for cchbar module
*/
#include <string>
#include "psi4/cc/common/CCParams.h"

#ifndef CCHBAR_PARAMS_H
#define CCHBAR_PARAMS_H

namespace psi {
namespace cchbar {

/*! \brief Parameters for cchbar module
 *
 * Extends common CC parameters with cchbar-specific parameters.
 */
struct Params : public psi::cc::common::CCParams {
    // Module-specific parameters
    int Tamplitude;      // Compute T-amplitude equation matrix elements
    int wabei_lowdisk;   // Use minimal-disk algorithm for Wabei (very slow)

    Params() : CCParams(), Tamplitude(0), wabei_lowdisk(0) {}
};

}  // namespace cchbar
}  // namespace psi

#endif
