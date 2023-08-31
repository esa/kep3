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

#include <kep3/core_astro/ic2eq2ic.hpp>

#include "catch.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using kep3::eq2ic;
using kep3::ic2eq;

constexpr double pi{boost::math::constants::pi<double>()};

TEST_CASE("ic2eq") {
  // Zero inclination and eccentricity
  {
    auto par = ic2eq({1.0, 0.0, 0.0, 0.0, 1.0, 0.0}, 1.0);
    REQUIRE(par[0] == 1.); // p is 1
    REQUIRE(par[1] == 0.); // f is zero
    REQUIRE(par[2] == 0.); // g is zero
    REQUIRE(par[0] == 1.); // h is 1
    REQUIRE(par[1] == 0.); // k is zero
    REQUIRE(par[2] == 0.); // L is zero
  }
  // Orbit at 90 degrees inclination
  {
    auto par = ic2eq({1.0, 0.0, 0.0, 0.0, 0.0, 1.1}, 1.0);
    REQUIRE_THAT(par[0],
                 WithinRel(1.2658227848101269 * (1. - 0.21 * 0.21), 1e-14));
    REQUIRE_THAT(par[1], WithinRel(0.21, 1e-14));
    REQUIRE(par[2] == 0.); // f is zero
    REQUIRE(par[3] == 1.); // h is 1
    REQUIRE(par[4] == 0.); // k is zero
    REQUIRE(par[5] == 0.); // L is zero
  }
  // Orbit at 90 degrees inclination
  {
    auto par = ic2eq({1.0, 0.0, 0.0, 0.0, 0.0, -1.1}, 1.0);
    REQUIRE_THAT(par[0],
                 WithinRel(1.2658227848101269 * (1. - 0.21 * 0.21), 1e-14));
    REQUIRE_THAT(par[1], WithinRel(0.21, 1e-14));
    REQUIRE(par[2] == 0.);  // f is zero
    REQUIRE(par[3] == -1.); // h is 1
    REQUIRE(par[4] == 0.);  // k is zero
    REQUIRE(par[5] == 0.);  // L is zero
  }
}

TEST_CASE("ic2eq2ic") {
  // Engines
  // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
  std::mt19937 rng_engine(122012203u);
  // Distributions for the elements
  std::uniform_real_distribution<double> sma_d(1.1, 100.);
  std::uniform_real_distribution<double> ecc_d(0, 0.99);
  std::uniform_real_distribution<double> incl_d(0., 3.);
  std::uniform_real_distribution<double> Omega_d(0, 2 * pi);
  std::uniform_real_distribution<double> omega_d(0., pi);
  std::uniform_real_distribution<double> ni_d(0, 2 * pi);
  {
    {
      // Testing on N random calls on ellipses
      unsigned N = 10000;
      for (auto i = 0u; i < N; ++i) {
        // We sample randomly on the Keplerian space
        double sma = sma_d(rng_engine);
        double ecc = ecc_d(rng_engine);
        double incl = incl_d(rng_engine);
        double Omega = Omega_d(rng_engine);
        double omega = omega_d(rng_engine);
        double ni = ni_d(rng_engine);
        // Compute the modified equinoctial
        double p = sma * (1. - ecc * ecc);
        double h = ecc * std::cos(Omega + omega);
        double k = ecc * std::sin(Omega + omega);
        double f = std::tan(incl / 2.) * std::cos(omega);
        double g = std::tan(incl / 2.) * std::sin(omega);
        double L = Omega + omega + f;
        auto pos_vel = eq2ic({p, h, k, f, g, L}, 1.);
        auto par = ic2eq(pos_vel, 1.0);
        REQUIRE_THAT(par[0], WithinAbs(p, 1e-10));
        REQUIRE_THAT(par[1], WithinAbs(h, 1e-10));
        REQUIRE_THAT(par[2], WithinAbs(k, 1e-10));
        REQUIRE_THAT(par[3], WithinRel(f, 1e-10));
        //REQUIRE_THAT(par[4], WithinAbs(g, 1e-10));
        //REQUIRE_THAT(par[5], WithinAbs(L, 1e-10));
      }
    }
  }
}