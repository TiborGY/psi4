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
    \brief Enter brief description of file here
*/

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/psi4-dec.h"
#include "MOInfo.h"
#include "Params.h"
#include "psi4/cc/common/CCParamsParser.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cchbar {

void get_params(Options &options) {
    // Parse common CC parameters using shared parser
    psi::cc::common::parse_memory(params.memory);
    psi::cc::common::parse_cachelev(params.cachelev, options);
    psi::cc::common::parse_print(params.print, options);
    psi::cc::common::parse_wfn(params.wfn, options);
    psi::cc::common::parse_dertype(params.dertype, options);

    // Module-specific parameters
    params.Tamplitude = options.get_bool("T_AMPS");
    params.wabei_lowdisk = options.get_bool("WABEI_LOWDISK");
}

}  // namespace cchbar
}  // namespace psi
