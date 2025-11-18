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

#ifndef CC_COMMON_CCPARAMSPARSER_H
#define CC_COMMON_CCPARAMSPARSER_H

/*! \file
    \ingroup CC
    \brief Parser functions for common CC parameters
*/

#include "CCParams.h"

namespace psi {
class Options;

namespace cc {
namespace common {

/*!
 * \brief Parse all common parameters from Options
 *
 * This is the main entry point for parsing common CC parameters.
 * It calls all individual parsing functions to populate the CCParams structure.
 *
 * \param params The CCParams structure to populate
 * \param options The Options object containing user input
 */
void parse_common_params(CCParams& params, Options& options);

/*!
 * \brief Parse the REFERENCE keyword and set ref parameter
 *
 * Parses "RHF", "ROHF", or "UHF" and sets ref to 0, 1, or 2 respectively.
 * Throws PsiException if invalid reference type is provided.
 *
 * \param ref Reference to the ref parameter to set
 * \param options The Options object
 */
void parse_ref(int& ref, Options& options);

/*!
 * \brief Parse the DERTYPE keyword and set dertype parameter
 *
 * Parses "NONE", "FIRST", or "RESPONSE" and sets dertype to 0, 1, or 3.
 * Throws PsiException if invalid derivative type is provided.
 *
 * \param dertype Reference to the dertype parameter to set
 * \param options The Options object
 */
void parse_dertype(int& dertype, Options& options);

/*!
 * \brief Parse the WFN keyword and set wfn parameter
 *
 * Gets the wavefunction/method type string (e.g., "CCSD", "CCSD_T").
 * Does not validate the wfn type (validation is module-specific).
 *
 * \param wfn Reference to the wfn parameter to set
 * \param options The Options object
 */
void parse_wfn(std::string& wfn, Options& options);

/*!
 * \brief Parse memory from environment and set memory parameter
 *
 * Gets available memory in bytes from Process::environment.
 *
 * \param memory Reference to the memory parameter to set
 */
void parse_memory(long int& memory);

/*!
 * \brief Parse the CACHELEVEL keyword and set cachelev parameter
 *
 * \param cachelev Reference to the cachelev parameter to set
 * \param options The Options object
 */
void parse_cachelev(int& cachelev, Options& options);

/*!
 * \brief Parse the RESTART keyword and set restart parameter
 *
 * \param restart Reference to the restart parameter to set
 * \param options The Options object
 */
void parse_restart(int& restart, Options& options);

/*!
 * \brief Parse the LOCAL keyword and set local parameter
 *
 * \param local Reference to the local parameter to set
 * \param options The Options object
 */
void parse_local(int& local, Options& options);

/*!
 * \brief Parse the ABCD keyword and set abcd parameter
 *
 * \param abcd Reference to the abcd parameter to set
 * \param options The Options object
 */
void parse_abcd(std::string& abcd, Options& options);

/*!
 * \brief Parse the PRINT keyword and set print parameter
 *
 * \param print Reference to the print parameter to set
 * \param options The Options object
 */
void parse_print(int& print, Options& options);

/*!
 * \brief Parse the MAXITER keyword and set maxiter parameter
 *
 * \param maxiter Reference to the maxiter parameter to set
 * \param options The Options object
 */
void parse_maxiter(int& maxiter, Options& options);

/*!
 * \brief Parse the R_CONVERGENCE keyword and set convergence parameter
 *
 * \param convergence Reference to the convergence parameter to set
 * \param options The Options object
 */
void parse_convergence(double& convergence, Options& options);

/*!
 * \brief Parse the DIIS keyword and set diis parameter
 *
 * \param diis Reference to the diis parameter to set
 * \param options The Options object
 */
void parse_diis(int& diis, Options& options);

/*!
 * \brief Parse the NUM_AMPS_PRINT keyword and set num_amps parameter
 *
 * \param num_amps Reference to the num_amps parameter to set
 * \param options The Options object
 */
void parse_num_amps(int& num_amps, Options& options);

}  // namespace common
}  // namespace cc
}  // namespace psi

#endif  // CC_COMMON_CCPARAMSPARSER_H
