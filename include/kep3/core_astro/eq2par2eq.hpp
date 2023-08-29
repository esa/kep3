// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// From cartesian to osculating Keplerian
/**
 * Transforms cartesian coordinates (r,v) to Keplerian elements (a,e,i,W,w,E).
 * Note that we use the eccentric anomaly (or Gudermannian if e > 1)
 */

#ifndef kep3_EQ2PAR_H
#define kep3_EQ2PAR_H

#include <array>

#include <kep3/detail/visibility.hpp>

namespace kep3 {

kep3_DLL_PUBLIC std::array<double, 6> eq2par(const std::array<double, 6> &eq,
                                             bool retrogade = false);

kep3_DLL_PUBLIC std::array<double, 6> par2eq(const std::array<double, 6> &par,
                                             bool retrogade = false);

} // namespace kep3
#endif // kep3_EQ2PAR_H