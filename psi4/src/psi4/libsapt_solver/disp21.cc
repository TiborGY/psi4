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

#include "sapt2p.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
namespace sapt {

void SAPT2p::disp21() {
    double e_disp210 = disp21_1(PSIF_SAPT_AMPS, "gARAR x tARBS", "tARBS Amplitudes", aoccA_, nvirA_, aoccB_, nvirB_);
    e_disp210 += disp21_2(PSIF_SAPT_AMPS, "T AR Intermediates", "Theta AR Intermediates", aoccA_, nvirA_);

    if (debug_) {
        outfile->Printf("    Disp210             = %18.12lf [Eh]\n", e_disp210);
    }

    double e_disp201 = disp21_1(PSIF_SAPT_AMPS, "gBSBS x tARBS", "tARBS Amplitudes", aoccA_, nvirA_, aoccB_, nvirB_);
    e_disp201 += disp21_2(PSIF_SAPT_AMPS, "T BS Intermediates", "Theta BS Intermediates", aoccB_, nvirB_);

    if (debug_) {
        outfile->Printf("    Disp201             = %18.12lf [Eh]\n\n", e_disp201);
    }

    e_disp21_ = e_disp210 + e_disp201;

    if (print_) {
        outfile->Printf("    Disp21              = %18.12lf [Eh]\n", e_disp21_);
    }
}

double SAPT2p::disp21_1(int ampfile, const char *glabel, const char *tlabel, size_t aoccA, size_t nvirA, size_t aoccB,
                        size_t nvirB) {
    auto tARBS = std::make_shared<Matrix>("tARBS", aoccA * nvirA, aoccB * nvirB);
    psio_->read_entry(ampfile, tlabel, (char *)tARBS->get_pointer(), sizeof(double) * aoccA * nvirA * aoccB * nvirB);

    auto gARBS = std::make_shared<Matrix>("gARBS", aoccA * nvirA, aoccB * nvirB);
    psio_->read_entry(ampfile, glabel, (char *)gARBS->get_pointer(), sizeof(double) * aoccA * nvirA * aoccB * nvirB);

    double energy = 4.0 * C_DDOT((long int)aoccA * nvirA * aoccB * nvirB, tARBS->get_pointer(), 1, gARBS->get_pointer(), 1);

    if (debug_) {
        outfile->Printf("\n    Disp21_1            = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}

double SAPT2p::disp21_2(int ampfile, const char *tlabel, const char *thetalabel, size_t aoccA, size_t nvirA) {
    auto T_p_AR = std::make_shared<Matrix>("T_p_AR", aoccA * nvirA, ndf_ + 3);
    psio_->read_entry(ampfile, tlabel, (char *)T_p_AR->get_pointer(), sizeof(double) * aoccA * nvirA * (ndf_ + 3));

    auto theta_p_AR = std::make_shared<Matrix>("theta_p_AR", aoccA * nvirA, ndf_ + 3);
    psio_->read_entry(ampfile, thetalabel, (char *)theta_p_AR->get_pointer(), sizeof(double) * aoccA * nvirA * (ndf_ + 3));

    double energy = 8.0 * C_DDOT((long int)aoccA * nvirA * (ndf_ + 3), T_p_AR->get_pointer(), 1, theta_p_AR->get_pointer(), 1);

    if (debug_) {
        outfile->Printf("    Disp21_2            = %18.12lf [Eh]\n", energy);
    }

    return (energy);
}
}  // namespace sapt
}  // namespace psi
