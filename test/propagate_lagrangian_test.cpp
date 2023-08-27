// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <functional>
#include <stdexcept>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <boost/math/constants/constants.hpp>

#include <kep3/core_astro/ic2par.hpp>
#include <kep3/core_astro/kepler_equations.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>

#include "catch.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using kep3::propagate_lagrangian;
using kep3::propagate_lagrangian_u;

constexpr double pi{boost::math::constants::pi<double>()};

void test_propagate_lagrangian(
    const std::function<void(std::array<std::array<double, 3>, 2> &, double,
                             double)> &propagate,
    unsigned int N = 10000) {
  { // Elliptical circular motion xy
    std::array<std::array<double, 3>, 2> pos_vel = {1., 0, 0., 0., 1., 0.};
    propagate(pos_vel, pi / 2., 1.);
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
    propagate(pos_vel, -pi / 2., 1.);
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
    propagate(pos_vel, pi / 2., 1.);
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
    propagate(pos_vel, pi / 2., 1.);
    auto &[pos, vel] = pos_vel;

    REQUIRE_THAT(pos[0], WithinAbs(0., 1e-14));
    REQUIRE_THAT(pos[1], WithinAbs(0., 1e-14));
    REQUIRE_THAT(pos[2], WithinAbs(1., 1e-14));
    REQUIRE_THAT(vel[0], WithinAbs(0., 1e-14));
    REQUIRE_THAT(vel[1], WithinAbs(-1., 1e-14));
    REQUIRE_THAT(vel[2], WithinAbs(0., 1e-14));
  }
  // We test orbital parameters are unchanged for random propagations
  // Engines
  // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
  std::mt19937 rng_engine(12201203u);
  // Distribution for the initial Cartesian state
  std::uniform_real_distribution<double> r_d(-2., 2.);
  std::uniform_real_distribution<double> v_d(-2., 2.);
  std::uniform_real_distribution<double> time(0.1, 20.);
  // Testing on N random calls
  for (auto i = 0u; i < N; ++i) {
    std::array<double, 3> pos = {r_d(rng_engine), r_d(rng_engine),
                                 r_d(rng_engine)};
    std::array<double, 3> vel = {v_d(rng_engine), v_d(rng_engine),
                                 v_d(rng_engine)};
    std::array<std::array<double, 3>, 2> pos_vel = {pos, vel};
    auto par_before = kep3::ic2par(pos_vel, 1.0);
    // We filter out cases of little significance (too close to singularity)
    if (std::abs(par_before[0]) > 0.5 && std::abs(par_before[0]) < 10. &&
        std::abs(1 - par_before[1]) > 1e-1) {
      propagate(pos_vel, time(rng_engine), 1.);
      auto par_after = kep3::ic2par(pos_vel, 1.0);
      fmt::print("\n\n{}", par_before);
      fmt::print("\n{}", par_after);

      REQUIRE_THAT(par_before[0], WithinAbs(par_after[0], 1e-11));
      REQUIRE_THAT(par_before[1], WithinAbs(par_after[1], 1e-11));
      REQUIRE_THAT(par_before[2], WithinAbs(par_after[2], 1e-11));
      REQUIRE_THAT(par_before[3], WithinAbs(par_after[3], 1e-11));
      REQUIRE_THAT(par_before[4], WithinAbs(par_after[4], 1e-11));
    }
  }
}

TEST_CASE("propagate_lagrangian") {
  //test_propagate_lagrangian(&propagate_lagrangian, 10000u);
  // test_propagate_lagrangian(&propagate_lagrangian_u, 10000u);
}