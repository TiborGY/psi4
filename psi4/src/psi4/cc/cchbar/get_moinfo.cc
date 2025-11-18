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
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/cc/ccmoinfo/CCMOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cchbar {

/*
** get_moinfo():  Routine to obtain basic orbital information from
** CC_INFO.
**
** T. Daniel Crawford, October 1996.
** Modified by TDC, March 1999.
** Migrated to use unified CCMOInfo, November 2025.
*/

void get_moinfo(std::shared_ptr<Wavefunction> wfn, Options &options) {
    int reference;

    // Read reference type from PSIO
    psio_read_entry(PSIF_CC_INFO, "Reference Wavefunction", (char *)&(params.ref), sizeof(int));

    // Allow ROHF EOM calculation after RHF energy
    std::string read_eom_ref = options.get_str("EOM_REFERENCE");
    if (read_eom_ref == "ROHF") params.ref = 1;

    // Store reference type for initialization
    reference = params.ref;

    // Initialize MOInfo using unified CCMOInfo class
    moinfo.initialize(wfn, reference);
}

/* Frees memory allocated in get_moinfo() and dumps some info. */
void cleanup() {
    // CCMOInfo handles cleanup automatically via destructor
    // No manual memory freeing needed anymore!
}

}  // namespace cchbar
}  // namespace psi
