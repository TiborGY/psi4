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

#include "integral_permutations.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psi4-dec.h"
#include <stdexcept>
#include <sstream>
#include <cctype>

namespace psi {
namespace libtrans {

//===========================================
// Internal Helper Implementations
//===========================================

std::string detect_spin_case(dpdbuf4 *InBuf) {
    if (InBuf == nullptr) {
        throw std::runtime_error("libtrans::detect_spin_case: null input buffer");
    }

    // Get the DPD parameter indices for bra and ket
    int pqnum = InBuf->params->pqnum;
    int rsnum = InBuf->params->rsnum;

    // Access DPD's pair information to determine spin case
    // The pair labels in DPD encode spin information through case:
    // Upper case (O, V) = alpha
    // Lower case (o, v) = beta

    // Get the row and column labels for the pair
    // This is a simplified approach - we examine the label directly
    std::string label(InBuf->file.label);

    // Count upper and lower case O/V characters
    bool has_upper = false;
    bool has_lower = false;

    for (char c : label) {
        if (c == 'O' || c == 'V') has_upper = true;
        if (c == 'o' || c == 'v') has_lower = true;
    }

    // Determine spin case based on character case
    if (has_upper && has_lower) {
        return "AB";  // Mixed alpha-beta
    } else if (has_lower) {
        return "BB";  // Beta-beta
    } else if (has_upper) {
        return "AA";  // Alpha-alpha
    } else {
        return "R";   // Restricted (no case distinction)
    }
}

std::string derive_oo_indices(dpdbuf4 *InBuf) {
    if (InBuf == nullptr) {
        throw std::runtime_error("libtrans::derive_oo_indices: null input buffer");
    }

    // Detect spin case from input buffer
    std::string spin_case = detect_spin_case(InBuf);

    // Return appropriate DPD index string based on spin case
    if (spin_case == "AA") {
        return "[O,O]";
    } else if (spin_case == "BB") {
        return "[o,o]";
    } else if (spin_case == "AB") {
        return "[O,o]";
    } else {  // Restricted
        return "[O,O]";
    }
}

std::string derive_vv_indices(dpdbuf4 *InBuf) {
    if (InBuf == nullptr) {
        throw std::runtime_error("libtrans::derive_vv_indices: null input buffer");
    }

    // Detect spin case from input buffer
    std::string spin_case = detect_spin_case(InBuf);

    // Return appropriate DPD index string based on spin case
    if (spin_case == "AA") {
        return "[V,V]";
    } else if (spin_case == "BB") {
        return "[v,v]";
    } else if (spin_case == "AB") {
        return "[V,v]";
    } else {  // Restricted
        return "[V,V]";
    }
}

std::string generate_standard_label(
    const std::string &bra_space,
    const std::string &ket_space)
{
    std::ostringstream oss;
    oss << "MO Ints <" << bra_space << "|" << ket_space << ">";
    return oss.str();
}

//===========================================
// Public API Implementations
//===========================================

void ovov_to_oovv(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label)
{
    // Validate input
    if (InBuf == nullptr) {
        throw std::runtime_error("libtrans::ovov_to_oovv: null input buffer");
    }

    // Auto-detect output indices
    std::string oo_indices = derive_oo_indices(InBuf);
    std::string vv_indices = derive_vv_indices(InBuf);

    // Generate label if not provided
    std::string final_label = label;
    if (label.empty()) {
        std::string spin_case = detect_spin_case(InBuf);
        std::string bra, ket;

        if (spin_case == "AA") {
            bra = "OO"; ket = "VV";
        } else if (spin_case == "BB") {
            bra = "oo"; ket = "vv";
        } else if (spin_case == "AB") {
            bra = "Oo"; ket = "Vv";
        } else {  // Restricted
            bra = "OO"; ket = "VV";
        }

        final_label = generate_standard_label(bra, ket);
    }

    // Perform transformation using DPD library
    // (OV|OV) → <OO|VV> using prqs permutation
    global_dpd_->buf4_sort(InBuf, outfilenum, prqs, oo_indices, vv_indices, final_label);
}

void ovov_to_oovv_and_close(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label)
{
    ovov_to_oovv(InBuf, outfilenum, label);
    global_dpd_->buf4_close(InBuf);
}

void vvoo_to_ovov(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label)
{
    // Validate input
    if (InBuf == nullptr) {
        throw std::runtime_error("libtrans::vvoo_to_ovov: null input buffer");
    }

    // For (VV|OO) → <OV|OV>, we need OV indices
    // This is the reverse of ovov_to_oovv
    std::string spin_case = detect_spin_case(InBuf);

    std::string ov_indices;
    std::string final_label = label;

    if (spin_case == "AA") {
        ov_indices = "[O,V]";
        if (label.empty()) final_label = generate_standard_label("OV", "OV");
    } else if (spin_case == "BB") {
        ov_indices = "[o,v]";
        if (label.empty()) final_label = generate_standard_label("ov", "ov");
    } else if (spin_case == "AB") {
        ov_indices = "[O,v]";
        if (label.empty()) final_label = generate_standard_label("Ov", "Ov");
    } else {  // Restricted
        ov_indices = "[O,V]";
        if (label.empty()) final_label = generate_standard_label("OV", "OV");
    }

    // Apply prqs permutation: (VV|OO) → (VO|OV) = <OV|OV>
    global_dpd_->buf4_sort(InBuf, outfilenum, prqs, ov_indices, ov_indices, final_label);
}

void vvoo_to_ovov_and_close(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label)
{
    vvoo_to_ovov(InBuf, outfilenum, label);
    global_dpd_->buf4_close(InBuf);
}

void transpose_oovv(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label)
{
    // Validate input
    if (InBuf == nullptr) {
        throw std::runtime_error("libtrans::transpose_oovv: null input buffer");
    }

    // Swap OO and VV to get <VV|OO>
    // Input should be <OO|VV> format (physicist's notation)

    // Read indices from input (they're already in the correct format)
    int oo_indices = InBuf->params->pqnum;  // Currently OO
    int vv_indices = InBuf->params->rsnum;  // Currently VV

    // Generate label if not provided
    std::string final_label = label;
    if (label.empty()) {
        // Detect spaces from input and reverse them
        std::string spin_case = detect_spin_case(InBuf);
        std::string bra, ket;

        if (spin_case == "AA") {
            bra = "VV"; ket = "OO";
        } else if (spin_case == "BB") {
            bra = "vv"; ket = "oo";
        } else if (spin_case == "AB") {
            bra = "Vv"; ket = "Oo";
        } else {  // Restricted
            bra = "VV"; ket = "OO";
        }

        final_label = generate_standard_label(bra, ket);
    }

    // Swap bra and ket using rspq permutation
    global_dpd_->buf4_sort(InBuf, outfilenum, rspq, vv_indices, oo_indices, final_label);
}

void transpose_bra_ket(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label)
{
    // Generic transpose for any integral type
    if (InBuf == nullptr) {
        throw std::runtime_error("libtrans::transpose_bra_ket: null input buffer");
    }

    int pq_indices = InBuf->params->pqnum;
    int rs_indices = InBuf->params->rsnum;

    // Auto-generate label if not provided by appending "transposed"
    std::string final_label = label;
    if (label.empty()) {
        std::ostringstream oss;
        oss << InBuf->file.label << " (transposed)";
        final_label = oss.str();
    }

    global_dpd_->buf4_sort(InBuf, outfilenum, rspq, rs_indices, pq_indices, final_label);
}

void vooo_to_ovoo(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label)
{
    // Validate input
    if (InBuf == nullptr) {
        throw std::runtime_error("libtrans::vooo_to_ovoo: null input buffer");
    }

    // For (VO|OO) → <OV|OO>, we swap within bra using qprs
    std::string spin_case = detect_spin_case(InBuf);

    std::string ov_indices, oo_indices;
    std::string final_label = label;

    if (spin_case == "AA") {
        ov_indices = "[O,V]";
        oo_indices = "[O,O]";
        if (label.empty()) final_label = generate_standard_label("OV", "OO");
    } else if (spin_case == "BB") {
        ov_indices = "[o,v]";
        oo_indices = "[o,o]";
        if (label.empty()) final_label = generate_standard_label("ov", "oo");
    } else if (spin_case == "AB") {
        ov_indices = "[O,v]";
        oo_indices = "[O,o]";
        if (label.empty()) final_label = generate_standard_label("Ov", "Oo");
    } else {  // Restricted
        ov_indices = "[O,V]";
        oo_indices = "[O,O]";
        if (label.empty()) final_label = generate_standard_label("OV", "OO");
    }

    // Apply qprs permutation: (VO|OO) → (OV|OO)
    global_dpd_->buf4_sort(InBuf, outfilenum, qprs, ov_indices, oo_indices, final_label);
}

void ovvv_to_ovvv(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label)
{
    // Validate input
    if (InBuf == nullptr) {
        throw std::runtime_error("libtrans::ovvv_to_ovvv: null input buffer");
    }

    // For (OV|VV) → <OV|VV>, this is typically a chemist to physicist conversion
    std::string spin_case = detect_spin_case(InBuf);

    std::string ov_indices, vv_indices;
    std::string final_label = label;

    if (spin_case == "AA") {
        ov_indices = "[O,V]";
        vv_indices = "[V,V]";
        if (label.empty()) final_label = generate_standard_label("OV", "VV");
    } else if (spin_case == "BB") {
        ov_indices = "[o,v]";
        vv_indices = "[v,v]";
        if (label.empty()) final_label = generate_standard_label("ov", "vv");
    } else if (spin_case == "AB") {
        ov_indices = "[O,v]";
        vv_indices = "[V,v]";
        if (label.empty()) final_label = generate_standard_label("Ov", "Vv");
    } else {  // Restricted
        ov_indices = "[O,V]";
        vv_indices = "[V,V]";
        if (label.empty()) final_label = generate_standard_label("OV", "VV");
    }

    // Apply prqs permutation: (OV|VV) → (OV|VV)
    // Note: This might be identity in some cases, but generally handles notation conversion
    global_dpd_->buf4_sort(InBuf, outfilenum, prqs, ov_indices, vv_indices, final_label);
}

void ovvv_to_ovvv_and_close(
    dpdbuf4 *InBuf,
    int outfilenum,
    const std::string &label)
{
    ovvv_to_ovvv(InBuf, outfilenum, label);
    global_dpd_->buf4_close(InBuf);
}

void chemist_to_physicist(
    dpdbuf4 *InBuf,
    int outfilenum,
    int pq_indices,
    int rs_indices,
    const std::string &label)
{
    // Low-level interface - requires explicit indices
    if (InBuf == nullptr) {
        throw std::runtime_error("libtrans::chemist_to_physicist: null input buffer");
    }

    // Perform transformation using DPD library
    // (pq|rs) → (pr|qs) using prqs permutation
    global_dpd_->buf4_sort(InBuf, outfilenum, prqs, pq_indices, rs_indices, label);
}

void chemist_to_physicist_and_close(
    dpdbuf4 *InBuf,
    int outfilenum,
    int pq_indices,
    int rs_indices,
    const std::string &label)
{
    chemist_to_physicist(InBuf, outfilenum, pq_indices, rs_indices, label);
    global_dpd_->buf4_close(InBuf);
}

void standard_ovov_transformations(
    dpdbuf4 *InBuf,
    int outfilenum,
    bool generate_transpose)
{
    // Validate input
    if (InBuf == nullptr) {
        throw std::runtime_error("libtrans::standard_ovov_transformations: null input buffer");
    }

    // Transform (OV|OV) → <OO|VV>
    ovov_to_oovv(InBuf, outfilenum);  // Uses auto-generated label

    if (generate_transpose) {
        // Open the just-created <OO|VV> buffer
        dpdbuf4 OOVVBuf;
        std::string oo_idx = derive_oo_indices(InBuf);
        std::string vv_idx = derive_vv_indices(InBuf);

        // Determine standard label for the intermediate buffer
        std::string spin_case = detect_spin_case(InBuf);
        std::string bra, ket;

        if (spin_case == "AA") {
            bra = "OO"; ket = "VV";
        } else if (spin_case == "BB") {
            bra = "oo"; ket = "vv";
        } else if (spin_case == "AB") {
            bra = "Oo"; ket = "Vv";
        } else {  // Restricted
            bra = "OO"; ket = "VV";
        }

        std::string oovv_label = generate_standard_label(bra, ket);

        global_dpd_->buf4_init(&OOVVBuf, outfilenum, InBuf->file.my_irrep,
                              oo_idx, vv_idx, 0, oovv_label);

        // Create <VV|OO> transpose
        transpose_oovv(&OOVVBuf, outfilenum);  // Auto-generates label

        global_dpd_->buf4_close(&OOVVBuf);
    }
}

void apply_permutation(
    dpdbuf4 *InBuf,
    int outfilenum,
    indices permutation,
    int pq_indices,
    int rs_indices,
    const std::string &label)
{
    // Direct wrapper for advanced use
    if (InBuf == nullptr) {
        throw std::runtime_error("libtrans::apply_permutation: null input buffer");
    }

    global_dpd_->buf4_sort(InBuf, outfilenum, permutation, pq_indices, rs_indices, label);
}

//===========================================
// Utility Functions
//===========================================

std::string get_permutation_name(indices permutation) {
    switch(permutation) {
        case prqs: return "Chemist to Physicist notation (prqs)";
        case rspq: return "Transpose bra and ket (rspq)";
        case psrq: return "Mixed permutation (ps|rq)";
        case qprs: return "Swap within bra (qprs)";
        case pqsr: return "Swap within ket (pqsr)";
        case rqps: return "Complex permutation (rqps)";
        default: return "Unknown permutation";
    }
}

bool is_permutation_supported(indices permutation) {
    // List of tested and supported permutations
    switch(permutation) {
        case prqs:
        case rspq:
        case psrq:
        case qprs:
        case pqsr:
        case rqps:
            return true;
        default:
            return false;
    }
}

} // namespace libtrans
} // namespace psi
