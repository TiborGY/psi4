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
    \brief DPD initialization utilities for coupled cluster modules
*/

#include "dpd_utils.h"
#include "CCMOInfo.h"
#include "psi4/libdpd/dpd.h"
#include <vector>

namespace psi {
namespace ccmoinfo {

void dpd_init_rhf(CCMOInfo& moinfo, int dpd_num, long int memory,
                  int cachetype, int* cachefiles, int** cachelist,
                  int* priority) {
    std::vector<int *> spaces;
    spaces.push_back(moinfo.occpi);
    spaces.push_back(moinfo.occ_sym.data());
    spaces.push_back(moinfo.virtpi);
    spaces.push_back(moinfo.vir_sym.data());

    dpd_init(dpd_num, moinfo.nirreps, memory, cachetype, cachefiles,
             cachelist, priority, 2, spaces);
}

void dpd_init_uhf(CCMOInfo& moinfo, int dpd_num, long int memory,
                  int cachetype, int* cachefiles, int** cachelist,
                  int* priority) {
    std::vector<int *> spaces;
    spaces.push_back(moinfo.aoccpi);
    spaces.push_back(moinfo.aocc_sym.data());
    spaces.push_back(moinfo.avirtpi);
    spaces.push_back(moinfo.avir_sym.data());
    spaces.push_back(moinfo.boccpi);
    spaces.push_back(moinfo.bocc_sym.data());
    spaces.push_back(moinfo.bvirtpi);
    spaces.push_back(moinfo.bvir_sym.data());

    dpd_init(dpd_num, moinfo.nirreps, memory, cachetype, cachefiles,
             cachelist, priority, 4, spaces);
}

void dpd_init_ao_rhf(CCMOInfo& moinfo, int dpd_num, long int memory,
                     int* cachefiles, int** cachelist) {
    std::vector<int *> aospaces;
    aospaces.push_back(moinfo.occpi);
    aospaces.push_back(moinfo.occ_sym.data());
    aospaces.push_back(moinfo.sopi);
    aospaces.push_back(moinfo.sosym.data());

    dpd_init(dpd_num, moinfo.nirreps, memory, 0, cachefiles,
             cachelist, nullptr, 2, aospaces);
}

void dpd_init_ao_uhf(CCMOInfo& moinfo, int dpd_num, long int memory,
                     int* cachefiles, int** cachelist) {
    std::vector<int *> aospaces;
    aospaces.push_back(moinfo.aoccpi);
    aospaces.push_back(moinfo.aocc_sym.data());
    aospaces.push_back(moinfo.sopi);
    aospaces.push_back(moinfo.sosym.data());
    aospaces.push_back(moinfo.boccpi);
    aospaces.push_back(moinfo.bocc_sym.data());
    aospaces.push_back(moinfo.sopi);
    aospaces.push_back(moinfo.sosym.data());

    dpd_init(dpd_num, moinfo.nirreps, memory, 0, cachefiles,
             cachelist, nullptr, 4, aospaces);
}

void dpd_init_from_moinfo(CCMOInfo& moinfo, int ref, int dpd_num,
                          long int memory, int cachetype, int* cachefiles,
                          int** cachelist, int* priority) {
    if (ref == 2) {
        // UHF
        dpd_init_uhf(moinfo, dpd_num, memory, cachetype, cachefiles,
                     cachelist, priority);
    } else {
        // RHF or ROHF
        dpd_init_rhf(moinfo, dpd_num, memory, cachetype, cachefiles,
                     cachelist, priority);
    }
}

}  // namespace ccmoinfo
}  // namespace psi
