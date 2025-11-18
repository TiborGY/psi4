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

#ifndef _dfocc_arrays_h_
#define _dfocc_arrays_h_

// This file now serves as a compatibility wrapper for the unified libarray
#include "psi4/libarray/arrays.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

// Import array classes from libarray into dfoccwave namespace for backward compatibility
using libarray::Array1d;
using libarray::Array2d;
using libarray::Array3d;
using libarray::Array1i;
using libarray::Array2i;
using libarray::Array3i;

}  // namespace dfoccwave
}  // namespace psi
#endif  // _dfocc_arrays_h_
