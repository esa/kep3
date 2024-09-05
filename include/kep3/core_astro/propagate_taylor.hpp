// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_PROPAGATE_TAYLOR_H
#define kep3_PROPAGATE_TAYLOR_H

#include <array>
#include <optional>
#include <utility>
#include <vector>

#include <kep3/detail/visibility.hpp>
#include <kep3/ta/lt_kepler.hpp>

namespace kep3
{

/// Taylor propagation for the Stark problem
/**
 * This function propagates an initial Cartesian state for a time t assuming a
 * central body and a keplerian motion perturbed by a constant inertial thrust. Heyoka Taylor adaptive integration is
 * used as a basic numerical technique. All units systems can be used, as long as the input parameters are all expressed
 * in the same system. The information on the state transition matrix can be optionally asked, in which case also the
 * sensitivities to the throttles are returned.
 */
kep3_DLL_PUBLIC std::tuple<std::array<std::array<double, 3>, 2>, double,
                          std::optional<std::pair<std::array<double, 49>, std::array<double, 21>>>>
propagate_taylor(const std::array<std::array<double, 3>, 2> &pos_vel, double m, std::array<double, 3> thrust,
                 double tof, double mu, double veff, double tol, bool stm = false);

// kep3_DLL_PUBLIC std::vector<std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>>
// propagate_taylor_v(const std::array<std::array<double, 3>, 2> &pos_vel, std::vector<double> tof, double mu,
// bool stm = false);
} // namespace kep3

#endif // kep3_PROPAGATE_TAYLOR_H