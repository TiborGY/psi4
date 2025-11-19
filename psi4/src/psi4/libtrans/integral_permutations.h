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

#ifndef _PSI_SRC_LIB_LIBTRANS_INTEGRAL_PERMUTATIONS_H_
#define _PSI_SRC_LIB_LIBTRANS_INTEGRAL_PERMUTATIONS_H_

#include "psi4/libdpd/dpd.h"
#include <string>

namespace psi {
namespace libtrans {

/**
 * @brief Utility functions for common integral permutation operations
 *
 * This namespace provides semantic wrappers around DPD buf4_sort operations
 * to make integral transformations more readable and maintainable.
 *
 * **Design Philosophy:**
 * - Functions infer DPD indices and labels from input buffer when possible
 * - Optional parameters allow customization when needed
 * - Standard naming conventions are enforced for consistency
 *
 * **Usage:**
 * @code
 * using libtrans;
 * ovov_to_oovv(&I, PSIF_LIBTRANS_DPD);
 * chemist_to_physicist_and_close(&T, PSIF_OCC_DPD, ID("[O,V]"), ID("[O,V]"), label);
 * @endcode
 */

//===========================================
// Common Integral Space Transformations
//===========================================

/**
 * @brief Transform (OV|OV) → <OO|VV>
 *
 * Common transformation for building 2-electron integrals in
 * occupied-occupied/virtual-virtual space.
 *
 * **Smart Defaults:**
 * - Automatically determines output indices from input buffer
 * - Generates standard label "MO Ints <OO|VV>" (or <oo|vv>, <Oo|Vv> for spin cases)
 * - Handles all spin cases (RHF, UHF-AA, UHF-BB, UHF-AB) automatically
 *
 * @param InBuf Input (OV|OV) buffer (reads DPD parameters to infer output)
 * @param outfilenum Output file number
 * @param label Optional custom label (auto-generated if empty)
 *
 * @code
 * using libtrans;
 * ovov_to_oovv(&I, PSIF_LIBTRANS_DPD);
 * ovov_to_oovv(&I, PSIF_LIBTRANS_DPD, "Custom <OO|VV>");
 * @endcode
 */
void ovov_to_oovv(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label = ""
);

/**
 * @brief Transform (OV|OV) → <OO|VV> and close input buffer
 *
 * Convenience function that performs the transformation and automatically
 * closes the input buffer, following the common pattern of immediate closure
 * after transformation.
 *
 * @param InBuf Input (OV|OV) buffer (will be closed after transformation)
 * @param outfilenum Output file number
 * @param label Optional custom label (auto-generated if empty)
 */
void ovov_to_oovv_and_close(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label = ""
);

/**
 * @brief Transform (VV|OO) → <OV|OV>
 *
 * Creates <OV|OV> integrals from (VV|OO) chemist's notation.
 * Auto-detects spin case and generates appropriate indices/label.
 *
 * @param InBuf Input (VV|OO) buffer
 * @param outfilenum Output file
 * @param label Optional custom label (default: "MO Ints <OV|OV>")
 */
void vvoo_to_ovov(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label = ""
);

/**
 * @brief Transform (VV|OO) → <OV|OV> and close input buffer
 *
 * @param InBuf Input (VV|OO) buffer (will be closed after transformation)
 * @param outfilenum Output file
 * @param label Optional custom label
 */
void vvoo_to_ovov_and_close(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label = ""
);

/**
 * @brief Create transpose: <OO|VV> → <VV|OO>
 *
 * Swaps bra and ket pairs. Auto-detects indices from input.
 *
 * @param InBuf Input <OO|VV> buffer
 * @param outfilenum Output file
 * @param label Optional custom label (default: "MO Ints <VV|OO>")
 */
void transpose_oovv(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label = ""
);

/**
 * @brief Generic transpose: <AB|CD> → <CD|AB>
 *
 * Swaps bra and ket pairs for any integral type.
 *
 * @param InBuf Input buffer
 * @param outfilenum Output file
 * @param label Optional custom label (auto-generated from input if empty)
 */
void transpose_bra_ket(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label = ""
);

//===========================================
// Additional Common Transformations
//===========================================

/**
 * @brief Transform (VO|OO) → <OV|OO>
 *
 * Common for OOOV-type integrals. Auto-detects spin case.
 *
 * @param InBuf Input (VO|OO) buffer
 * @param outfilenum Output file
 * @param label Optional custom label
 */
void vooo_to_ovoo(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label = ""
);

/**
 * @brief Transform (OV|VV) → <OV|VV>
 *
 * Common for OVVV-type integrals. Auto-detects spin case.
 *
 * @param InBuf Input (OV|VV) buffer
 * @param outfilenum Output file
 * @param label Optional custom label
 */
void ovvv_to_ovvv(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label = ""
);

/**
 * @brief Transform (OV|VV) → <OV|VV> and close input buffer
 *
 * @param InBuf Input (OV|VV) buffer (will be closed after transformation)
 * @param outfilenum Output file
 * @param label Optional custom label
 */
void ovvv_to_ovvv_and_close(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label = ""
);

//===========================================
// Notation Conversions (General)
//===========================================

/**
 * @brief Convert chemist's notation (pq|rs) to physicist's notation <pr|qs>
 *
 * Applies the prqs permutation: (pq|rs) → (pr|qs)
 * This is the most common transformation in quantum chemistry codes.
 *
 * Unlike the specific functions above, this requires explicit indices
 * because the space types cannot be inferred from generic (pq|rs).
 *
 * @param InBuf Input buffer in chemist's notation
 * @param outfilenum PSI file number for output
 * @param pq_indices DPD indices for bra (pr)
 * @param rs_indices DPD indices for ket (qs)
 * @param label String label for output buffer
 */
void chemist_to_physicist(
    dpdbuf4 *InBuf,
    int outfilenum,
    int pq_indices,
    int rs_indices,
    const std::string &label
);

/**
 * @brief Convert chemist's notation to physicist's notation and close input buffer
 *
 * Convenience function that performs the transformation and automatically
 * closes the input buffer.
 *
 * @param InBuf Input buffer in chemist's notation (will be closed after transformation)
 * @param outfilenum PSI file number for output
 * @param pq_indices DPD indices for bra (pr)
 * @param rs_indices DPD indices for ket (qs)
 * @param label String label for output buffer
 */
void chemist_to_physicist_and_close(
    dpdbuf4 *InBuf,
    int outfilenum,
    int pq_indices,
    int rs_indices,
    const std::string &label
);

//===========================================
// Batch Operations
//===========================================

/**
 * @brief Perform standard OVOV transformations in batch
 *
 * Common pattern: Generate all permutations needed for OVOV integrals:
 * - (OV|OV) → <OO|VV>
 * - <OO|VV> → <VV|OO> (if requested)
 *
 * All indices and labels are auto-generated.
 *
 * @param InBuf Input (OV|OV) buffer
 * @param outfilenum Output file
 * @param generate_transpose Also create <VV|OO>
 */
void standard_ovov_transformations(
    dpdbuf4 *InBuf,
    int outfilenum,
    bool generate_transpose = true
);

//===========================================
// Low-Level Interface (for advanced use)
//===========================================

/**
 * @brief Apply arbitrary permutation with full control
 *
 * Direct wrapper around DPD buf4_sort for cases where automatic
 * inference isn't suitable. Most code should use the semantic
 * functions above instead.
 *
 * @param InBuf Input buffer
 * @param outfilenum Output file
 * @param permutation DPD indices enum (prqs, rspq, etc.)
 * @param pq_indices New bra indices
 * @param rs_indices New ket indices
 * @param label Output label
 */
void apply_permutation(
    dpdbuf4 *InBuf,
    int outfilenum,
    indices permutation,
    int pq_indices,
    int rs_indices,
    const std::string &label
);

//===========================================
// Utilities
//===========================================

/**
 * @brief Get human-readable name for permutation
 *
 * @param permutation DPD indices enum
 * @return String describing the permutation (e.g., "Chemist to Physicist")
 */
std::string get_permutation_name(indices permutation);

/**
 * @brief Validate permutation is supported
 *
 * @param permutation DPD indices enum
 * @return true if permutation is implemented and tested
 */
bool is_permutation_supported(indices permutation);

//===========================================
// Internal Helpers
//===========================================

/**
 * @brief Derive output OO indices from input OV indices
 *
 * Analyzes input buffer's DPD parameters to determine appropriate
 * output indices for OO space, handling spin cases correctly.
 * @return DPD index string like "[O,O]", "[o,o]", or "[O,o]"
 */
std::string derive_oo_indices(dpdbuf4 *InBuf);

/**
 * @brief Derive output VV indices from input OV indices
 * @return DPD index string like "[V,V]", "[v,v]", or "[V,v]"
 */
std::string derive_vv_indices(dpdbuf4 *InBuf);

/**
 * @brief Generate standard label for integral type
 *
 * @param bra_space Space designation for bra (e.g., "OO", "oo", "Oo")
 * @param ket_space Space designation for ket (e.g., "VV", "vv", "Vv")
 * @return Standard label like "MO Ints <OO|VV>"
 */
std::string generate_standard_label(
    const std::string &bra_space,
    const std::string &ket_space
);

/**
 * @brief Detect spin case from DPD indices
 *
 * @return String: "AA", "BB", or "AB" for unrestricted; "R" for restricted
 */
std::string detect_spin_case(dpdbuf4 *InBuf);

} // namespace libtrans
} // namespace psi

#endif
