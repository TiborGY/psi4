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
    \ingroup CC
    \brief Implementation of parser functions for common CC parameters
*/

#include "CCParamsParser.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/psi4-dec.h"

#include <string>

namespace psi {
namespace cc {
namespace common {

void parse_common_params(CCParams& params, Options& options) {
    // Universal parameters (all 7 modules)
    parse_wfn(params.wfn, options);
    parse_ref(params.ref, options);

    // Common parameters (6/7 modules)
    parse_memory(params.memory);
    parse_cachelev(params.cachelev, options);

    // Derivative type (6/7 modules - not in cceom)
    if (options["DERTYPE"].has_changed()) {
        parse_dertype(params.dertype, options);
    }

    // Frequently used parameters (parse if keyword is present)
    if (options["RESTART"].has_changed()) {
        parse_restart(params.restart, options);
    }

    if (options["LOCAL"].has_changed()) {
        parse_local(params.local, options);
    }

    if (options["ABCD"].has_changed()) {
        parse_abcd(params.abcd, options);
    }

    if (options["PRINT"].has_changed()) {
        parse_print(params.print, options);
    }

    if (options["MAXITER"].has_changed()) {
        parse_maxiter(params.maxiter, options);
    }

    if (options["R_CONVERGENCE"].has_changed()) {
        parse_convergence(params.convergence, options);
    }

    if (options["DIIS"].has_changed()) {
        parse_diis(params.diis, options);
    }

    if (options["NUM_AMPS_PRINT"].has_changed()) {
        parse_num_amps(params.num_amps, options);
    }
}

void parse_ref(int& ref, Options& options) {
    std::string junk = options.get_str("REFERENCE");

    if (junk == "RHF") {
        ref = 0;
    } else if (junk == "ROHF") {
        ref = 1;
    } else if (junk == "UHF") {
        ref = 2;
    } else {
        throw PsiException("Invalid value of input keyword REFERENCE", __FILE__, __LINE__);
    }
}

void parse_dertype(int& dertype, Options& options) {
    std::string junk = options.get_str("DERTYPE");

    if (junk == "NONE") {
        dertype = 0;
    } else if (junk == "FIRST") {
        dertype = 1;
    } else if (junk == "RESPONSE") {
        dertype = 3;  // linear response
    } else {
        throw PsiException("Invalid value of input keyword DERTYPE", __FILE__, __LINE__);
    }
}

void parse_wfn(std::string& wfn, Options& options) {
    wfn = options.get_str("WFN");
}

void parse_memory(long int& memory) {
    memory = Process::environment.get_memory();
}

void parse_cachelev(int& cachelev, Options& options) {
    cachelev = options.get_int("CACHELEVEL");
}

void parse_restart(int& restart, Options& options) {
    restart = options.get_bool("RESTART");
}

void parse_local(int& local, Options& options) {
    local = options.get_bool("LOCAL");
}

void parse_abcd(std::string& abcd, Options& options) {
    abcd = options.get_str("ABCD");
}

void parse_print(int& print, Options& options) {
    print = options.get_int("PRINT");
}

void parse_maxiter(int& maxiter, Options& options) {
    maxiter = options.get_int("MAXITER");
}

void parse_convergence(double& convergence, Options& options) {
    convergence = options.get_double("R_CONVERGENCE");
}

void parse_diis(int& diis, Options& options) {
    diis = options.get_bool("DIIS");
}

void parse_num_amps(int& num_amps, Options& options) {
    num_amps = options.get_int("NUM_AMPS_PRINT");
}

}  // namespace common
}  // namespace cc
}  // namespace psi
