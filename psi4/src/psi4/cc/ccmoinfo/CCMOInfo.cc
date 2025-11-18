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
    \brief Implementation of unified CCMOInfo class
*/

#include "CCMOInfo.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/psifiles.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <cstring>

namespace psi {
namespace ccmoinfo {

CCMOInfo::CCMOInfo()
    : nirreps(0),
      nmo(0),
      nso(0),
      nao(0),
      nvirt(0),
      enuc(0.0),
      escf(0.0),
      eref(0.0),
      ecc(0.0),
      iter(0),
      conv(0.0),
      sym(0),
      emp2(0.0),
      emp2_ss(0.0),
      emp2_os(0.0),
      emp2_s(0.0),
      ecc_ss(0.0),
      ecc_os(0.0),
      ecc_s(0.0),
      t1diag(0.0),
      d1diag(0.0),
      new_d1diag(0.0),
      d2diag(0.0),
      lcc(0.0),
      et(0.0),
      iopen(0),
      nactive(0),
      nfzc(0),
      nfzv(0),
      nclsd(0),
      nopen(0),
      nuocc(0),
      noei(0),
      ntri(0),
      noei_ao(0),
      natom(0),
      Cv(nullptr),
      Cav(nullptr),
      Cbv(nullptr),
      Co(nullptr),
      Cao(nullptr),
      Cbo(nullptr),
      C(nullptr),
      Ca(nullptr),
      Cb(nullptr),
      reference_type_(0) {}

CCMOInfo::~CCMOInfo() {
    free_transformation_matrices();
}

void CCMOInfo::initialize(std::shared_ptr<Wavefunction> wfn, int reference) {
    reference_type_ = reference;

    // Read common data from wavefunction
    read_common_data(wfn, reference);

    // Read orbital information from PSIO
    read_orbital_info(reference);

    // Allocate transformation matrices if needed
    allocate_transformation_matrices(wfn, reference);
}

void CCMOInfo::read_common_data(std::shared_ptr<Wavefunction> wfn, int reference) {
    // Basic dimensions from wavefunction
    nirreps = wfn->nirrep();
    nmo = wfn->nmo();
    nso = wfn->nso();
    nao = wfn->basisset()->nao();

    // Dimension arrays
    sopi = wfn->nsopi();
    orbspi = wfn->nmopi();
    frdocc = wfn->frzcpi();
    fruocc = wfn->frzvpi();
    openpi = wfn->soccpi();
    clsdpi = wfn->doccpi() - frdocc;
    uoccpi = orbspi - clsdpi - openpi - fruocc - frdocc;

    // Labels
    labels = wfn->molecule()->irrep_labels();

    // Energies from wavefunction
    enuc = wfn->molecule()->nuclear_repulsion_energy(wfn->get_dipole_field_strength());
    escf = wfn->reference_wavefunction() ? wfn->reference_wavefunction()->energy() : wfn->energy();

    // Initialize convergence tracking
    conv = 0.0;

    // Read reference energy from PSIO
    psio_read_entry(PSIF_CC_INFO, "Reference Energy", (char *)&eref, sizeof(double));

    // Compute active orbital counts based on reference type
    if (reference == 2) { // UHF
        aoccpi = clsdpi + wfn->soccpi();
        boccpi = clsdpi;
        avirtpi = uoccpi;
        bvirtpi = uoccpi + wfn->soccpi();
    } else { // RHF or ROHF
        occpi = clsdpi + wfn->soccpi();
        virtpi = uoccpi + wfn->soccpi();
        act_occpi = occpi;  // For ccresponse
    }

    // Compute total virtual orbitals for RHF
    if (reference == 0) {
        nvirt = virtpi.sum();
    }

    // Build sosym array (for AO-basis routines)
    sosym.resize(nso);
    for (int h = 0, q = 0; h < nirreps; h++) {
        for (int p = 0; p < sopi[h]; p++) {
            sosym[q++] = h;
        }
    }
}

void CCMOInfo::read_orbital_info(int reference) {
    int nactive;

    // Read number of active orbitals
    psio_read_entry(PSIF_CC_INFO, "No. of Active Orbitals", (char *)&nactive, sizeof(int));
    this->nactive = nactive;

    if (reference == 2) { // UHF
        // Resize symmetry arrays
        aocc_sym.resize(nactive);
        bocc_sym.resize(nactive);
        avir_sym.resize(nactive);
        bvir_sym.resize(nactive);

        // Read symmetry arrays
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Occ Orb Symmetry", (char *)aocc_sym.data(),
                       sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Occ Orb Symmetry", (char *)bocc_sym.data(),
                       sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Virt Orb Symmetry", (char *)avir_sym.data(),
                       sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Virt Orb Symmetry", (char *)bvir_sym.data(),
                       sizeof(int) * nactive);

        // Resize offset arrays
        aocc_off.resize(nirreps);
        bocc_off.resize(nirreps);
        avir_off.resize(nirreps);
        bvir_off.resize(nirreps);

        // Read offset arrays
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Occ Orb Offsets", (char *)aocc_off.data(),
                       sizeof(int) * nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Occ Orb Offsets", (char *)bocc_off.data(),
                       sizeof(int) * nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Virt Orb Offsets", (char *)avir_off.data(),
                       sizeof(int) * nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Virt Orb Offsets", (char *)bvir_off.data(),
                       sizeof(int) * nirreps);

    } else { // RHF or ROHF
        // Resize symmetry arrays
        occ_sym.resize(nactive);
        vir_sym.resize(nactive);

        // Read symmetry arrays
        psio_read_entry(PSIF_CC_INFO, "Active Occ Orb Symmetry", (char *)occ_sym.data(),
                       sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orb Symmetry", (char *)vir_sym.data(),
                       sizeof(int) * nactive);

        // Resize offset arrays
        occ_off.resize(nirreps);
        vir_off.resize(nirreps);

        // Read offset arrays
        psio_read_entry(PSIF_CC_INFO, "Active Occ Orb Offsets", (char *)occ_off.data(),
                       sizeof(int) * nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orb Offsets", (char *)vir_off.data(),
                       sizeof(int) * nirreps);
    }
}

void CCMOInfo::allocate_transformation_matrices(std::shared_ptr<Wavefunction> wfn, int reference) {
    psio_address next;

    if (reference == 0 || reference == 1) { // RHF or ROHF
        // Allocate and read occupied orbital transformation matrix
        Co = (double ***)malloc(nirreps * sizeof(double **));
        next = PSIO_ZERO;
        for (int h = 0; h < nirreps; h++) {
            Co[h] = nullptr;
            if (sopi[h] && occpi[h]) {
                Co[h] = block_matrix(sopi[h], occpi[h]);
                psio_read(PSIF_CC_INFO, "RHF/ROHF Active Occupied Orbitals", (char *)Co[h][0],
                         sizeof(double) * sopi[h] * occpi[h], next, &next);
            }
        }

        // Allocate and read virtual orbital transformation matrix
        Cv = (double ***)malloc(nirreps * sizeof(double **));
        next = PSIO_ZERO;
        for (int h = 0; h < nirreps; h++) {
            Cv[h] = nullptr;
            if (sopi[h] && virtpi[h]) {
                Cv[h] = block_matrix(sopi[h], virtpi[h]);
                psio_read(PSIF_CC_INFO, "RHF/ROHF Active Virtual Orbitals", (char *)Cv[h][0],
                         sizeof(double) * sopi[h] * virtpi[h], next, &next);
            }
        }

    } else if (reference == 2) { // UHF
        // Allocate and read alpha virtual orbital transformation matrix
        Cav = (double ***)malloc(nirreps * sizeof(double **));
        next = PSIO_ZERO;
        for (int h = 0; h < nirreps; h++) {
            Cav[h] = nullptr;
            if (sopi[h] && avirtpi[h]) {
                Cav[h] = block_matrix(sopi[h], avirtpi[h]);
                psio_read(PSIF_CC_INFO, "UHF Active Alpha Virtual Orbs", (char *)Cav[h][0],
                         sizeof(double) * sopi[h] * avirtpi[h], next, &next);
            }
        }

        // Allocate and read beta virtual orbital transformation matrix
        Cbv = (double ***)malloc(nirreps * sizeof(double **));
        next = PSIO_ZERO;
        for (int h = 0; h < nirreps; h++) {
            Cbv[h] = nullptr;
            if (sopi[h] && bvirtpi[h]) {
                Cbv[h] = block_matrix(sopi[h], bvirtpi[h]);
                psio_read(PSIF_CC_INFO, "UHF Active Beta Virtual Orbs", (char *)Cbv[h][0],
                         sizeof(double) * sopi[h] * bvirtpi[h], next, &next);
            }
        }

        // Note: UHF occupied transformation matrices (Cao, Cbo) would be allocated here
        // if needed by specific modules. For now, they're left as nullptr.
    }
}

void CCMOInfo::free_transformation_matrices() {
    // Free RHF/ROHF matrices
    if (reference_type_ == 0 || reference_type_ == 1) {
        if (Co != nullptr) {
            for (int h = 0; h < nirreps; h++) {
                if (Co[h] != nullptr) {
                    free_block(Co[h]);
                }
            }
            free(Co);
            Co = nullptr;
        }

        if (Cv != nullptr) {
            for (int h = 0; h < nirreps; h++) {
                if (Cv[h] != nullptr) {
                    free_block(Cv[h]);
                }
            }
            free(Cv);
            Cv = nullptr;
        }
    }

    // Free UHF matrices
    if (reference_type_ == 2) {
        if (Cav != nullptr) {
            for (int h = 0; h < nirreps; h++) {
                if (Cav[h] != nullptr) {
                    free_block(Cav[h]);
                }
            }
            free(Cav);
            Cav = nullptr;
        }

        if (Cbv != nullptr) {
            for (int h = 0; h < nirreps; h++) {
                if (Cbv[h] != nullptr) {
                    free_block(Cbv[h]);
                }
            }
            free(Cbv);
            Cbv = nullptr;
        }

        if (Cao != nullptr) {
            for (int h = 0; h < nirreps; h++) {
                if (Cao[h] != nullptr) {
                    free_block(Cao[h]);
                }
            }
            free(Cao);
            Cao = nullptr;
        }

        if (Cbo != nullptr) {
            for (int h = 0; h < nirreps; h++) {
                if (Cbo[h] != nullptr) {
                    free_block(Cbo[h]);
                }
            }
            free(Cbo);
            Cbo = nullptr;
        }
    }

    // Free alternative naming convention matrices (C, Ca, Cb)
    // These typically alias Cv, Cav, Cbv, so we don't free them separately
    // unless they were allocated independently
    C = nullptr;
    Ca = nullptr;
    Cb = nullptr;
}

}  // namespace ccmoinfo
}  // namespace psi
