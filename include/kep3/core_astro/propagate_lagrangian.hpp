// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_PROPAGATE_LAGRANGIAN_H
#define kep3_PROPAGATE_LAGRANGIAN_H

#include <array>
#include <optional>

#include <kep3/detail/visibility.hpp>
#include <kep3/core_astro/stm.hpp>


namespace kep3
{

/// Lagrangian propagation
/**
 * This function propagates an initial Cartesian state for a time t assuming a
 * central body and a keplerian motion. Lagrange coefficients are used as basic
 * numerical technique. All units systems can be used, as long
 * as the input parameters are all expressed in the same system.
 */

kep3_DLL_PUBLIC std::optional<std::array<double, 36>>
propagate_lagrangian(std::array<std::array<double, 3>, 2> &pos_vel, double tof, double mu, bool stm = false);

kep3_DLL_PUBLIC void propagate_lagrangian_u(std::array<std::array<double, 3>, 2> &pos_vel, double dt, double mu);

kep3_DLL_PUBLIC void propagate_keplerian(std::array<std::array<double, 3>, 2> &pos_vel, double dt, double mu);
} // namespace kep3

#endif // KEP_TOOLBOX_PROPAGATE_LAGRANGIAN_H