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

#ifndef CCMOINFO_DPD_UTILS_H
#define CCMOINFO_DPD_UTILS_H

/*! \file
    \ingroup CC
    \brief DPD initialization utilities for coupled cluster modules

    This file provides convenience functions for initializing DPD (Distributed
    Packed Data) structures from CCMOInfo, eliminating ~70 lines of duplicate
    boilerplate across the 7 CC modules.
*/

namespace psi {
namespace ccmoinfo {

class CCMOInfo;

/**
 * @brief Initialize DPD for RHF/ROHF references
 *
 * Sets up DPD with occupied and virtual orbital spaces for closed-shell
 * or restricted open-shell references.
 *
 * @param moinfo Unified molecular orbital information
 * @param dpd_num DPD instance number (typically 0)
 * @param memory Memory available for DPD caching (in bytes)
 * @param cachetype Cache type (0 = default)
 * @param cachefiles Array of cache file unit numbers
 * @param cachelist 2D array specifying which integrals to cache
 * @param priority Optional priority list for cache management
 */
void dpd_init_rhf(CCMOInfo& moinfo, int dpd_num, long int memory,
                  int cachetype, int* cachefiles, int** cachelist,
                  int* priority = nullptr);

/**
 * @brief Initialize DPD for UHF references
 *
 * Sets up DPD with separate alpha/beta occupied and virtual orbital spaces
 * for unrestricted references.
 *
 * @param moinfo Unified molecular orbital information
 * @param dpd_num DPD instance number (typically 0)
 * @param memory Memory available for DPD caching (in bytes)
 * @param cachetype Cache type (0 = default)
 * @param cachefiles Array of cache file unit numbers
 * @param cachelist 2D array specifying which integrals to cache
 * @param priority Optional priority list for cache management
 */
void dpd_init_uhf(CCMOInfo& moinfo, int dpd_num, long int memory,
                  int cachetype, int* cachefiles, int** cachelist,
                  int* priority = nullptr);

/**
 * @brief Initialize AO-basis DPD for RHF/ROHF references
 *
 * Sets up a secondary DPD instance for AO-basis algorithms with
 * occupied MO and symmetry orbital (SO) spaces.
 *
 * @param moinfo Unified molecular orbital information
 * @param dpd_num DPD instance number (typically 1 for AO-basis)
 * @param memory Memory available for DPD caching (in bytes)
 * @param cachefiles Array of cache file unit numbers
 * @param cachelist 2D array specifying which integrals to cache
 */
void dpd_init_ao_rhf(CCMOInfo& moinfo, int dpd_num, long int memory,
                     int* cachefiles, int** cachelist);

/**
 * @brief Initialize AO-basis DPD for UHF references
 *
 * Sets up a secondary DPD instance for AO-basis algorithms with
 * separate alpha/beta occupied MO and symmetry orbital (SO) spaces.
 *
 * @param moinfo Unified molecular orbital information
 * @param dpd_num DPD instance number (typically 1 for AO-basis)
 * @param memory Memory available for DPD caching (in bytes)
 * @param cachefiles Array of cache file unit numbers
 * @param cachelist 2D array specifying which integrals to cache
 */
void dpd_init_ao_uhf(CCMOInfo& moinfo, int dpd_num, long int memory,
                     int* cachefiles, int** cachelist);

/**
 * @brief Convenience wrapper for DPD initialization
 *
 * Automatically selects RHF or UHF initialization based on reference type.
 *
 * @param moinfo Unified molecular orbital information
 * @param ref Reference type (0=RHF, 1=ROHF, 2=UHF)
 * @param dpd_num DPD instance number (typically 0)
 * @param memory Memory available for DPD caching (in bytes)
 * @param cachetype Cache type (0 = default)
 * @param cachefiles Array of cache file unit numbers
 * @param cachelist 2D array specifying which integrals to cache
 * @param priority Optional priority list for cache management
 */
void dpd_init_from_moinfo(CCMOInfo& moinfo, int ref, int dpd_num,
                          long int memory, int cachetype, int* cachefiles,
                          int** cachelist, int* priority = nullptr);

}  // namespace ccmoinfo
}  // namespace psi

#endif  // CCMOINFO_DPD_UTILS_H
