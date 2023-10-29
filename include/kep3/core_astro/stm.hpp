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

#ifndef kep3_STM_H
#define kep3_STM_H

#include <array>
#include <utility>

#include <kep3/detail/visibility.hpp>

namespace kep3
{
// From:
// Reynolds, Reid G. "Direct Solution of the Keplerian State Transition Matrix." Journal of Guidance, Control, and
// Dynamics 45, no. 6 (2022): 1162-1165.
kep3_DLL_PUBLIC std::array<double, 36> stm_reynolds(const std::array<std::array<double, 3>, 2> &pos_vel0,
                                           const std::array<std::array<double, 3>, 2> &pos_vel, double tof,
                                           double mu = 1.);

kep3_DLL_PUBLIC std::pair<std::array<std::array<double, 3>, 2>, std::array<double, 36>>
propagate_stm_reynolds(const std::array<std::array<double, 3>, 2> &pos_vel0, double tof, double mu = 1.);

kep3_DLL_PUBLIC std::pair<std::array<std::array<double, 3>, 2>, std::array<double, 36>>
propagate_stm2(const std::array<std::array<double, 3>, 2> &pos_vel0, double tof, double mu = 1.);

} // namespace kep3
#endif // kep3_IC2EQ2IC_H