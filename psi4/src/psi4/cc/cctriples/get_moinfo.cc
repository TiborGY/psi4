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
    \ingroup CCTRIPLES
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>

#include "psi4/libpsi4util/process.h"
#include "psi4/libpsi4util/PsiOutStream.h"
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
namespace cctriples {

/*
 ** get_moinfo():  Routine to obtain basic orbital information from
 ** CC_INFO.
 **
 ** T. Daniel Crawford, October 1996.
 ** Modified by TDC, March 1999.
 ** Migrated to use unified CCMOInfo, November 2025.
 */

void get_moinfo(std::shared_ptr<Wavefunction> wfn, Options &options) {
    std::string junk;
    int reference;

    // Parse WFN parameter
    params.wfn = options.get_str("WFN");
    if (params.wfn != "CCSD" && params.wfn != "CCSD_T" && params.wfn != "CCSD_AT" && params.wfn != "BCCD" &&
        params.wfn != "BCCD_T") {
        throw PsiException("Invalid value of input keyword WFN", __FILE__, __LINE__);
    }

    // Get number of threads
    params.nthreads = Process::environment.get_n_threads();
    if (options["CC_NUM_THREADS"].has_changed()) {
        params.nthreads = options.get_int("CC_NUM_THREADS");
    }

    // Parse reference type
    params.semicanonical = 0;
    junk = options.get_str("REFERENCE");
    /* if no reference is given, assume rhf */
    if (junk == "RHF")
        params.ref = 0;
    else if (junk == "ROHF" && (params.wfn == "CCSD_T" || params.wfn == "BCCD_T")) {
        params.ref = 2;
        params.semicanonical = 1;
    } else if (junk == "ROHF")
        params.ref = 1;
    else if (junk == "UHF")
        params.ref = 2;
    else {
        throw PsiException("Invalid value of input keyword REFERENCE", __FILE__, __LINE__);
    }

    // Parse derivative type
    junk = options.get_str("DERTYPE");
    if (junk == "NONE")
        params.dertype = 0;
    else if (junk == "FIRST")
        params.dertype = 1;
    else {
        throw PsiException("Value of keyword DERTYPE is not applicable to CCSD(T)", __FILE__, __LINE__);
    }

    // Store reference type for initialization
    reference = params.ref;

    // Initialize MOInfo using unified CCMOInfo class
    moinfo.initialize(wfn, reference);

    // Read CCSD energy from PSIO (CCMOInfo only reads reference energy)
    psio_read_entry(PSIF_CC_INFO, "CCSD Energy", (char *)&(moinfo.ecc), sizeof(double));

    // Print summary
    outfile->Printf("\n\n");
    outfile->Printf("    Wave function   =    %6s\n", params.wfn.c_str());
    if (params.semicanonical) {
        outfile->Printf("    Reference wfn   =    ROHF changed to UHF for Semicanonical Orbitals\n");
    } else {
        outfile->Printf("    Reference wfn   =    %5s\n",
                        (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
    }

    outfile->Printf("\n    Nuclear Rep. energy (wfn)                = %20.15f\n", moinfo.enuc);
    outfile->Printf("    SCF energy          (wfn)                = %20.15f\n", moinfo.escf);
    outfile->Printf("    Reference energy    (file100)            = %20.15f\n", moinfo.eref);
    outfile->Printf("    CCSD energy         (file100)            = %20.15f\n", moinfo.ecc);
    outfile->Printf("    Total CCSD energy   (file100)            = %20.15f\n", moinfo.eref + moinfo.ecc);
}

/* Frees memory allocated in get_moinfo() and dumps some info. */
void cleanup() {
    // CCMOInfo handles cleanup automatically via destructor
    // No manual memory freeing needed anymore!
}

}  // namespace cctriples
}  // namespace psi
