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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/psi4-dec.h"
#include "psi4/cc/ccmoinfo/CCMOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {

/* get_moinfo(): Routine to obtain basic orbital information from
** CC_INFO.
**
** T. Daniel Crawford, October 1996
** Modified for ccresponse by TDC May, 2003
** Migrated to use unified CCMOInfo, November 2025
*/

void get_moinfo(std::shared_ptr<Wavefunction> wfn) {
    int reference;

    // Read reference type from PSIO
    psio_read_entry(PSIF_CC_INFO, "Reference Wavefunction", (char *)&(params.ref), sizeof(int));
    reference = params.ref;

    // Initialize MOInfo using unified CCMOInfo class
    moinfo.initialize(wfn, reference);

    // ccresponse-specific: arrange active SCF MO's
    moinfo.Ca_matrix = wfn->Ca_subset("SO", "ACTIVE");
}

/* Frees memory allocated in get_moinfo() and dumps out the energy. */
void cleanup() {
    // mu_irreps and l_irreps are freed if they were allocated (in get_params.cc)
    if (moinfo.mu_irreps != nullptr) {
        free(moinfo.mu_irreps);
        moinfo.mu_irreps = nullptr;
    }
    if (moinfo.l_irreps != nullptr) {
        free(moinfo.l_irreps);
        moinfo.l_irreps = nullptr;
    }

    // CCMOInfo handles cleanup automatically via destructor
}

}  // namespace ccresponse
}  // namespace psi
