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

#ifndef kep3_FLYBY_H
#define kep3_FLYBY_H

#include <array>
#include <utility>

#include <kep3/detail/visibility.hpp>
#include <kep3/planet.hpp>

namespace kep3
{

kep3_DLL_PUBLIC std::pair<double, double> fb_con(const std::array<double, 3> &v_rel_in,
                                                 const std::array<double, 3> &v_rel_out, double, double);

kep3_DLL_PUBLIC double fb_dv(const std::array<double, 3> &v_rel_in, const std::array<double, 3> &v_rel_out, double,
                             double);

kep3_DLL_PUBLIC double fb_vout(const std::array<double, 3> &v_in, const std::array<double, 3> &v_pla, double rp,
                               double beta, double mu);

} // namespace kep3
#endif // kep3_FLYBY_H