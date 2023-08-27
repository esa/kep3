// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdexcept>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <boost/math/constants/constants.hpp>

#include <kep3/core_astro/kepler_equations.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>

#include "catch.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using kep3::propagate_lagrangian;
using kep3::propagate_lagrangian_u;

constexpr double pi{boost::math::constants::pi<double>()};

TEST_CASE("propagate_lagrangian") {
  { // Elliptical circular motion xy
    std::array<std::array<double, 3>, 2> pos_vel = {1., 0, 0., 0., 1., 0.};
    propagate_lagrangian(pos_vel, pi / 2., 1.);
    auto &[pos, vel] = pos_vel;

    REQUIRE_THAT(pos[0], WithinAbs(0., 1e-14));
    REQUIRE_THAT(pos[1], WithinAbs(1., 1e-14));
    REQUIRE_THAT(pos[2], WithinAbs(0., 1e-14));
    REQUIRE_THAT(vel[0], WithinAbs(-1., 1e-14));
    REQUIRE_THAT(vel[1], WithinAbs(0., 1e-14));
    REQUIRE_THAT(vel[2], WithinAbs(0., 1e-14));
  }
  { // Elliptical circular motion xy
    std::array<std::array<double, 3>, 2> pos_vel = {1., 0, 0., 0., 1., 0.};
    propagate_lagrangian(pos_vel, - pi / 2., 1.);
    auto &[pos, vel] = pos_vel;

    REQUIRE_THAT(pos[0], WithinAbs(0., 1e-14));
    REQUIRE_THAT(pos[1], WithinAbs(-1., 1e-14));
    REQUIRE_THAT(pos[2], WithinAbs(0., 1e-14));
    REQUIRE_THAT(vel[0], WithinAbs(1., 1e-14));
    REQUIRE_THAT(vel[1], WithinAbs(0., 1e-14));
    REQUIRE_THAT(vel[2], WithinAbs(0., 1e-14));
  }
  { // Elliptical circular motion xz
    std::array<std::array<double, 3>, 2> pos_vel = {1., 0, 0., 0., 0., 1.};
    propagate_lagrangian(pos_vel, pi / 2., 1.);
    auto &[pos, vel] = pos_vel;

    REQUIRE_THAT(pos[0], WithinAbs(0., 1e-14));
    REQUIRE_THAT(pos[1], WithinAbs(0., 1e-14));
    REQUIRE_THAT(pos[2], WithinAbs(1., 1e-14));
    REQUIRE_THAT(vel[0], WithinAbs(-1., 1e-14));
    REQUIRE_THAT(vel[1], WithinAbs(0., 1e-14));
    REQUIRE_THAT(vel[2], WithinAbs(0., 1e-14));
  }
  { // Elliptical circular motion yz
    std::array<std::array<double, 3>, 2> pos_vel = {0., 1, 0., 0., 0., 1.};
    propagate_lagrangian(pos_vel, pi / 2., 1.);
    auto &[pos, vel] = pos_vel;

    REQUIRE_THAT(pos[0], WithinAbs(0., 1e-14));
    REQUIRE_THAT(pos[1], WithinAbs(0., 1e-14));
    REQUIRE_THAT(pos[2], WithinAbs(1., 1e-14));
    REQUIRE_THAT(vel[0], WithinAbs(0., 1e-14));
    REQUIRE_THAT(vel[1], WithinAbs(-1., 1e-14));
    REQUIRE_THAT(vel[2], WithinAbs(0., 1e-14));
  }
}