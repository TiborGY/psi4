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
    \brief Enter brief description of file here
*/

#include "Params.h"
#include "psi4/cc/ccwave.h"

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace psi {
namespace ccenergy {

/*
** get_moinfo():  Routine to obtain basic orbital information from
** CC_INFO.
**
** T. Daniel Crawford, October 1996
** Modified by TDC, March 1999
** Migrated to use unified CCMOInfo, November 2025
*/

void CCEnergyWavefunction::get_moinfo() {
    int reference;

    // Read reference type from PSIO
    psio_read_entry(PSIF_CC_INFO, "Reference Wavefunction", (char *)&(params_.ref), sizeof(int));
    reference = params_.ref;

    // Initialize MOInfo using unified CCMOInfo class
    // Pass shared_from_this() to provide the wavefunction
    moinfo_.initialize(shared_from_this(), reference);

    // Store active orbital dimensions (ccenergy-specific)
    if (reference == 0 || reference == 1) {
        act_occpi_ = moinfo_.occpi;
        act_virpi_ = moinfo_.virtpi;
    }

    // Print summary
    outfile->Printf("\n    Nuclear Rep. energy (wfn)     = %20.15f\n", moinfo_.enuc);
    outfile->Printf("    SCF energy          (wfn)     = %20.15f\n", moinfo_.escf);
    outfile->Printf("    Reference energy    (file100) = %20.15f\n", moinfo_.eref);
}

/* Frees memory allocated in get_moinfo() and dumps out the energy. */
void CCEnergyWavefunction::cleanup() {
    if (params_.wfn == "CC2" || params_.wfn == "EOM_CC2")
        psio_write_entry(PSIF_CC_INFO, "CC2 Energy", (char *)&(moinfo_.ecc), sizeof(double));
    else if (params_.wfn == "CC3" || params_.wfn == "EOM_CC3")
        psio_write_entry(PSIF_CC_INFO, "CC3 Energy", (char *)&(moinfo_.ecc), sizeof(double));
    else
        psio_write_entry(PSIF_CC_INFO, "CCSD Energy", (char *)&(moinfo_.ecc), sizeof(double));

    // CCMOInfo handles cleanup automatically via destructor
    // No manual memory freeing needed anymore!
}

}  // namespace ccenergy
}  // namespace psi
