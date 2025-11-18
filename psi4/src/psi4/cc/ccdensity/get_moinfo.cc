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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {

/*
** get_moinfo():  Routine to obtain basic orbital information from
** CC_INFO.
**
** T. Daniel Crawford, October 1996.
** Modified by TDC, March 1999.
** Migrated to use unified CCMOInfo, November 2025
*/

void get_moinfo(std::shared_ptr<Wavefunction> wfn) {
    int i, j;

    // Initialize MOInfo using unified CCMOInfo class
    moinfo.initialize(wfn, params.ref);

    // Compute state symmetry (ccdensity-specific)
    moinfo.sym = 0;
    for (i = 0; i < moinfo.nirreps; ++i)
        for (j = 0; j < moinfo.openpi[i]; ++j)
            moinfo.sym = moinfo.sym ^ i;

    // Print energy summary
    outfile->Printf("\n\tNuclear Rep. energy (wfn)     = %20.15f\n", moinfo.enuc);
    outfile->Printf("\tSCF energy          (wfn)     = %20.15f\n", moinfo.escf);
    outfile->Printf("\tReference energy    (file100) = %20.15f\n", moinfo.eref);

    if (params.wfn == "CC2" || params.wfn == "EOM_CC2") {
        psio_read_entry(PSIF_CC_INFO, "CC2 Energy", (char *)&(moinfo.ecc), sizeof(double));
        outfile->Printf("\tCC2 energy          (CC_INFO) = %20.15f\n", moinfo.ecc);
        outfile->Printf("\tTotal CC2 energy    (CC_INFO) = %20.15f\n", moinfo.eref + moinfo.ecc);
    } else if (params.wfn == "CCSD" || params.wfn == "EOM_CCSD") {
        psio_read_entry(PSIF_CC_INFO, "CCSD Energy", (char *)&(moinfo.ecc), sizeof(double));
        outfile->Printf("\tCCSD energy         (CC_INFO) = %20.15f\n", moinfo.ecc);
        outfile->Printf("\tTotal CCSD energy   (CC_INFO) = %20.15f\n", moinfo.eref + moinfo.ecc);
    } else if (params.wfn == "CCSD_T") {
        psio_read_entry(PSIF_CC_INFO, "CCSD Energy", (char *)&(moinfo.ecc), sizeof(double));
        psio_read_entry(PSIF_CC_INFO, "(T) Energy", (char *)&(moinfo.et), sizeof(double));
        outfile->Printf("\tCCSD energy         (CC_INFO) = %20.15f\n", moinfo.ecc);
        outfile->Printf("\t(T) energy          (CC_INFO) = %20.15f\n", moinfo.et);
        outfile->Printf("\tTotal CCSD(T) energy(CC_INFO) = %20.15f\n", moinfo.eref + moinfo.ecc + moinfo.et);
    } else if (params.wfn == "CC3" || params.wfn == "EOM_CC3") {
        psio_read_entry(PSIF_CC_INFO, "CC3 Energy", (char *)&(moinfo.ecc), sizeof(double));
        outfile->Printf("\tCC3 energy          (CC_INFO) = %20.15f\n", moinfo.ecc);
        outfile->Printf("\tTotal CC3 energy    (CC_INFO) = %20.15f\n", moinfo.eref + moinfo.ecc);
    }
}

}  // namespace ccdensity
}  // namespace psi
