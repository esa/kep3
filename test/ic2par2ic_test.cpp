// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <boost/math/constants/constants.hpp>

#include <kep3/core_astro/ic2par.hpp>
#include <kep3/core_astro/par2ic.hpp>
#include <stdexcept>

#include "catch.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using kep3::ic2par;
using kep3::par2ic;

constexpr double pi{boost::math::constants::pi<double>()};

TEST_CASE("ic2par") {
  // Singular orbit zero inclination and eccentricity
  {
    auto par = ic2par({1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0);
    REQUIRE(par[0] == 1.); // sma is 1
    REQUIRE(par[1] == 0.); // ecc is zero
    REQUIRE(par[2] == 0.); // incl is zero
    REQUIRE(!std::isfinite(par[3]));
    REQUIRE(!std::isfinite(par[4]));
    REQUIRE(!std::isfinite(par[5]));
  }
  // Orbit at 90 degrees inclination
  {
    auto par = ic2par({1.0, 0.0, 0.0}, {0.0, 0.0, 1.1}, 1.0);
    REQUIRE_THAT(par[0], WithinRel(1.2658227848101269, 1e-14));
    REQUIRE_THAT(par[1], WithinRel(0.21, 1e-14));
    REQUIRE_THAT(par[2], WithinRel(pi / 2, 1e-14)); // incl at 90 degrees
    REQUIRE(par[3] == 0.);                          // Omega is zero
    REQUIRE(par[4] == 0.);                          // omega is zero
    REQUIRE(par[5] == 0.);                          // true anomaly is zero
  }
  // Orbit at 90 degrees inclination
  {
    auto par = ic2par({1.0, 0.0, 0.0}, {0.0, 0.0, -1.1}, 1.0);
    REQUIRE_THAT(par[0], WithinRel(1.2658227848101269, 1e-14));
    REQUIRE_THAT(par[1], WithinRel(0.21, 1e-14));
    REQUIRE_THAT(par[2], WithinRel(pi / 2, 1e-14)); // incl at 90 degrees
    REQUIRE_THAT(par[3], WithinRel(pi, 1e-14));     // Omega at 180 degrees
    REQUIRE_THAT(par[4], WithinRel(pi, 1e-14));     // omeg at 180 degrees
    REQUIRE(par[5] == 0.);                          // true anomaly is zero
  }
  // Retrogade orbit
  {
    auto par = ic2par({1.0, 0.0, 0.0}, {0.0, -1.0, 0.1}, 1.0);
    REQUIRE_THAT(par[0], WithinRel(1.01010101010101, 1e-14));
    REQUIRE_THAT(par[1], WithinRel(0.01, 1e-14));
    REQUIRE_THAT(par[2], WithinRel(174.289406862500 / 180.0 * pi,
                                   1e-14)); // incl
    REQUIRE(par[3] == 0.);                  // Omega is zero
    REQUIRE(par[4] == 0.);                  // omega is zero
    REQUIRE(par[5] == 0.);                  // true anomaly is zero
  }
  // A random orbit
  {
    auto par = ic2par({-1.1823467354129, 0.0247369349235, -0.014848484784},
                      {0.00232349642367, 1.1225625625625, -0.34678634567}, 1.0);
    REQUIRE_THAT(par[0], WithinRel(3.21921322281178, 1e-14));
    REQUIRE_THAT(par[1], WithinAbs(0.63283595179672, 1e-14));
    REQUIRE_THAT(par[2], WithinRel(162.82902986040048 / 180.0 * pi,
                                   1e-14)); // incl
    REQUIRE_THAT(par[3], WithinAbs(1.13023105373051 / 180.0 * pi,
                                   1e-14)); // Omega
    REQUIRE_THAT(par[4], WithinRel(179.22698703370386 / 180.0 * pi,
                                   1e-14)); // omega
    REQUIRE_THAT(par[5], WithinAbs(3.21031850920605 / 180.0 * pi,
                                   1e-14)); // true anomaly
  }
}

TEST_CASE("par2ic") {
  // Orbit at 90 degrees inclination
  {
    auto [r, v] =
        par2ic({1.2658227848101269, 0.21, pi / 2, 0.0, 0.0, 0.0}, 1.0);
    REQUIRE_THAT(r[0], WithinRel(1., 1e-14));
    REQUIRE_THAT(r[1], WithinAbs(0., 1e-14));
    REQUIRE_THAT(r[2], WithinAbs(0., 1e-14));
    REQUIRE_THAT(v[0], WithinAbs(0., 1e-14));
    REQUIRE_THAT(v[1], WithinAbs(0., 1e-14));
    REQUIRE_THAT(v[2], WithinRel(1.1, 1e-14));
  }
  // Orbit at 90 degrees inclination
  {
    auto [r, v] = par2ic({1.2658227848101269, 0.21, pi / 2, pi, pi, 0.0}, 1.0);
    REQUIRE_THAT(r[0], WithinRel(1., 1e-14));
    REQUIRE_THAT(r[1], WithinAbs(0., 1e-14));
    REQUIRE_THAT(r[2], WithinAbs(0., 1e-14));
    REQUIRE_THAT(v[0], WithinAbs(0., 1e-14));
    REQUIRE_THAT(v[1], WithinAbs(0., 1e-14));
    REQUIRE_THAT(v[2], WithinRel(-1.1, 1e-14));
  }
  // Retrogade orbit
  {
    auto [r, v] = par2ic(
        {1.01010101010101, 0.01, 174.289406862500 / 180.0 * pi, 1e-14, 0., 0.0},
        1.0);
    REQUIRE_THAT(r[0], WithinRel(1., 1e-14));
    REQUIRE_THAT(r[1], WithinAbs(0., 1e-14));
    REQUIRE_THAT(r[2], WithinAbs(0., 1e-14));
    REQUIRE_THAT(v[0], WithinAbs(0., 1e-14));
    REQUIRE_THAT(v[1], WithinRel(-1., 1e-14));
    REQUIRE_THAT(v[2], WithinAbs(0.1, 1e-14));
  }
  // A random orbit
  {
    auto [r, v] =
        par2ic({3.21921322281178, 0.63283595179672,
                162.82902986040048 / 180.0 * pi, 1.13023105373051 / 180.0 * pi,
                179.22698703370386 / 180.0 * pi, 3.21031850920605 / 180.0 * pi},
               1.0);

    REQUIRE_THAT(r[0], WithinRel(-1.1823467354129, 1e-14));
    REQUIRE_THAT(r[1], WithinAbs(0.0247369349235, 1e-13));
    REQUIRE_THAT(r[2], WithinAbs(-0.014848484784, 1e-14));
    REQUIRE_THAT(v[0], WithinAbs(0.00232349642367, 1e-13));
    REQUIRE_THAT(v[1], WithinRel(1.1225625625625, 1e-14));
    REQUIRE_THAT(v[2], WithinAbs(-0.34678634567, 1e-14));
  }
  // We check the convention a<0 -> e>1 is followed.
  REQUIRE_THROWS_AS(par2ic({1.3, 1.3, 0., 0., 0., 0.}, 1), std::domain_error);
  REQUIRE_THROWS_AS(par2ic({-1.3, 0.4, 0., 0., 0., 0.}, 1), std::domain_error);
  REQUIRE_THROWS_AS(par2ic({11.1, 1.4, 0., 0., 0., 5.23}, 1), std::domain_error);

}
TEST_CASE("ic2par2ic") {
  // Engines
  // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
  std::mt19937 rng_engine(122012203u);
  // Distributions for the elements
  std::uniform_real_distribution<double> sma_d(1.1, 100.);
  std::uniform_real_distribution<double> ecc_d(0, 0.99);
  std::uniform_real_distribution<double> incl_d(0., pi);
  std::uniform_real_distribution<double> Omega_d(0, 2 * pi);
  std::uniform_real_distribution<double> omega_d(0., pi);
  std::uniform_real_distribution<double> f_d(0, 2 * pi);
  {
    // Testing on N random calls on ellipses
    unsigned N = 10000;
    for (auto i = 0u; i < N; ++i) {
      double sma = sma_d(rng_engine);
      double ecc = ecc_d(rng_engine);
      double incl = incl_d(rng_engine);
      double Omega = Omega_d(rng_engine);
      double omega = omega_d(rng_engine);
      double f = f_d(rng_engine);
      auto [r, v] = par2ic({sma, ecc, incl, Omega, omega, f}, 1.);
      auto par = ic2par(r, v, 1.0);
      REQUIRE_THAT(par[0], WithinAbs(sma, 1e-10));
    }
  }
  {
    // Testing on N random calls on hyperbolas
    unsigned N = 10000;
    for (auto i = 0u; i < N; ++i) {
      double sma = -sma_d(rng_engine);
      double ecc = ecc_d(rng_engine) + 1.1;
      double incl = incl_d(rng_engine);
      double Omega = Omega_d(rng_engine);
      double omega = omega_d(rng_engine);
      double f = f_d(rng_engine);
      // Skipping if true anomaly out of asymptotes
      if (std::cos(f) < -1 / ecc) {
        continue;
      }
      auto [r, v] = par2ic({sma, ecc, incl, Omega, omega, f}, 1.);
      auto par = ic2par(r, v, 1.0);
      REQUIRE_THAT(par[0], WithinAbs(sma, 1e-10));
    }
  }
}
