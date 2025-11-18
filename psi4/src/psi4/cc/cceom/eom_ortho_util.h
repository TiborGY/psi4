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

/*!
** \file eom_ortho_util.h
** \brief DPD-aware orthogonalization utilities for EOM-CC calculations
** \ingroup CCEOM
**
** This file provides orthogonalization utilities specifically for EOM-CC
** calculations that work with DPD (Distributed Paired Data) objects.
** These functions follow the same algorithmic patterns as the general
** OrthoUtil library (psi4/libqt/ortho_util.h) but are adapted for the
** specialized data structures used in coupled cluster theory.
**
** Relationship to OrthoUtil:
** - dot_product_C() is analogous to OrthoUtil::dot_product()
** - orthogonalize_C() is analogous to OrthoUtil::orthogonalize_vector()
** - normalize_C() is analogous to OrthoUtil::normalize_vector()
** - schmidt_add() is analogous to OrthoUtil::schmidt_add()
**
** These functions implement Modified Gram-Schmidt orthogonalization for
** EOM vectors stored as combinations of dpdfile2 (singles) and dpdbuf4
** (doubles) objects.
*/

#pragma once

#include "psi4/libdpd/dpd.h"

namespace psi {
namespace cceom {

/**
 * @brief Compute dot product of two EOM vectors (DPD version)
 *
 * Computes the dot product of two EOM vectors, where each vector consists of
 * singles (file2) and doubles (buf4) components for both alpha and beta spins.
 *
 * This is the DPD equivalent of OrthoUtil::dot_product().
 *
 * For ROHF/UHF: result = <RIA|CME> + <Ria|Cme> + <RIJAB|CMNEF> + <Rijab|Cmnef> + <RIjAb|CMnEf>
 *
 * @param RIA Alpha singles amplitudes (first vector)
 * @param Ria Beta singles amplitudes (first vector)
 * @param RIJAB Alpha-alpha doubles amplitudes (first vector)
 * @param Rijab Beta-beta doubles amplitudes (first vector)
 * @param RIjAb Alpha-beta doubles amplitudes (first vector)
 * @param CME Alpha singles amplitudes (second vector)
 * @param Cme Beta singles amplitudes (second vector)
 * @param CMNEF Alpha-alpha doubles amplitudes (second vector)
 * @param Cmnef Beta-beta doubles amplitudes (second vector)
 * @param CMnEf Alpha-beta doubles amplitudes (second vector)
 * @return Dot product <R|C>
 */
double dot_product_C(dpdfile2 *RIA, dpdfile2 *Ria, dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb,
                     dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);

/**
 * @brief Orthogonalize EOM vector against a single basis vector
 *
 * Performs R = R - <R|C> * C, removing the component of R along C.
 * This is one step in the Modified Gram-Schmidt process.
 *
 * This is the DPD equivalent of the orthogonalization step in
 * OrthoUtil::orthogonalize_vector().
 *
 * @param RIA Alpha singles amplitudes (vector to orthogonalize, modified in place)
 * @param Ria Beta singles amplitudes (vector to orthogonalize, modified in place)
 * @param RIJAB Alpha-alpha doubles amplitudes (vector to orthogonalize, modified in place)
 * @param Rijab Beta-beta doubles amplitudes (vector to orthogonalize, modified in place)
 * @param RIjAb Alpha-beta doubles amplitudes (vector to orthogonalize, modified in place)
 * @param CME Alpha singles amplitudes (basis vector)
 * @param Cme Beta singles amplitudes (basis vector)
 * @param CMNEF Alpha-alpha doubles amplitudes (basis vector)
 * @param Cmnef Beta-beta doubles amplitudes (basis vector)
 * @param CMnEf Alpha-beta doubles amplitudes (basis vector)
 * @param overlap The computed overlap <R|C> (output)
 */
void orthogonalize_C(dpdfile2 *RIA, dpdfile2 *Ria, dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb,
                     dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf,
                     double *overlap);

/**
 * @brief RHF version of dot product for spin-adapted vectors
 *
 * Computes dot product for RHF EOM vectors with proper spin adaptation.
 * For RHF: result = 2<RIA|CME> + <2RIjAb - RIjbA|CMnEf>
 *
 * @param RIA Singles amplitudes (first vector)
 * @param RIjAb Doubles amplitudes (first vector)
 * @param CME Singles amplitudes (second vector)
 * @param CMnEf Doubles amplitudes (second vector)
 * @param irrep Irrep for temporary storage
 * @return Dot product
 */
double dot_product_C_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, dpdfile2 *CME, dpdbuf4 *CMnEf, int irrep);

/**
 * @brief RHF version of orthogonalization
 *
 * Orthogonalizes RHF EOM vector against a basis vector with proper spin adaptation.
 *
 * @param RIA Singles amplitudes (vector to orthogonalize, modified in place)
 * @param RIjAb Doubles amplitudes (vector to orthogonalize, modified in place)
 * @param CME Singles amplitudes (basis vector)
 * @param CMnEf Doubles amplitudes (basis vector)
 * @param overlap The computed overlap (output)
 * @param irrep Irrep for temporary storage
 */
void orthogonalize_C_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, dpdfile2 *CME, dpdbuf4 *CMnEf, double *overlap, int irrep);

}  // namespace cceom
}  // namespace psi
