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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/cc/ccmoinfo/CCMOInfo.h"
#include "cclambda.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cclambda {

/*
** get_moinfo(wfn):  Routine to obtain basic orbital information from
** the passed wavefunction object.
**
** T. Daniel Crawford, October 1996.
** Modified by TDC, March 1999.
** Modified for UHF references by TDC, June 2002.
** Migrated to use unified CCMOInfo, November 2025.
*/

void CCLambdaWavefunction::get_moinfo(std::shared_ptr<Wavefunction> wfn) {
    int i, j;
    int reference;

    // Read reference type from PSIO
    psio_read_entry(PSIF_CC_INFO, "Reference Wavefunction", (char *)&(params.ref), sizeof(int));
    reference = params.ref;

    // Initialize MOInfo using unified CCMOInfo class
    moinfo.initialize(wfn, reference);

    // Compute state symmetry (cclambda-specific)
    int sym = 0;
    for (i = 0; i < moinfo.nirreps; ++i)
        for (j = 0; j < moinfo.openpi[i]; ++j)
            sym = sym ^ i;
    moinfo.sym = sym;

    // Print summary
    outfile->Printf("\n\tNuclear Rep. energy (wfn)     = %20.15f\n", moinfo.enuc);
    outfile->Printf("\tReference           (wfn)     = %20d\n", params.ref);
    outfile->Printf("\tSCF energy          (wfn)     = %20.15f\n", moinfo.escf);
    outfile->Printf("\tReference energy    (CC_INFO) = %20.15f\n", moinfo.eref);
}

/* Frees memory allocated in get_moinfo() and dumps some info. */
void CCLambdaWavefunction::cleanup() {
    // Write lambda pseudoenergy to PSIO
    psio_write_entry(PSIF_CC_INFO, "Lambda Pseudoenergy", (char *)&(moinfo.lcc), sizeof(double));

    // CCMOInfo handles cleanup automatically via destructor
    // No manual memory freeing needed anymore!
}

}  // namespace cclambda
}  // namespace psi
