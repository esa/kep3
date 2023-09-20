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

#include <kep3/core_astro/constants.hpp>
#include <kep3/planets/keplerian.hpp>

#include "catch.hpp"

TEST_CASE("constructor") {
  REQUIRE_NOTHROW(kep3::udpla::keplerian{});
  kep3::epoch ref_epoch{12.22, kep3::epoch::MJD2000};
  // From posvel
  REQUIRE_NOTHROW(kep3::udpla::keplerian{
      ref_epoch, {{{0.3, 1., 0.2}, {0.0, 1.12, 0.}}}, 1., "unknown"});
  REQUIRE_NOTHROW(kep3::udpla::keplerian{ref_epoch,
                                         {{{0.3, 1., 0.2}, {0.0, 1.12, 0.}}},
                                         1.,
                                         "unknown",
                                         {-1, -1, -1}});
  // From parameters
  std::array<double, 6> par{{1., 0., 0., 0., 0., 0.}};
  REQUIRE_NOTHROW(kep3::udpla::keplerian{ref_epoch, par, 1., "unknown"});
  REQUIRE_NOTHROW(
      kep3::udpla::keplerian{ref_epoch, par, 1., "unknown", {-1, -1, -1}});
  // Checking the data members initializations:
  kep3::udpla::keplerian udpla{ref_epoch, par, 1.1, "unknown", {1.2, 2.2, 1.9}};
  REQUIRE(udpla.get_ref_epoch() == ref_epoch);
  REQUIRE(udpla.get_name() == "unknown");
  REQUIRE(udpla.get_mu_central_body() == 1.1);
  REQUIRE(udpla.get_mu_self() == 1.2);
  REQUIRE(udpla.get_radius() == 2.2);
  REQUIRE(udpla.get_safe_radius() == 1.9);
}

TEST_CASE("eph") {
  kep3::epoch ref_epoch{0., kep3::epoch::MJD2000};
  kep3::udpla::keplerian udpla1{
      ref_epoch,
      {{{kep3::AU, 0., 0.}, {0., kep3::EARTH_VELOCITY, 0.}}},
      kep3::MU_SUN};
  kep3::udpla::keplerian udpla2{
      ref_epoch,
      {{{kep3::AU, 0., 0.}, {0., -kep3::EARTH_VELOCITY, 0.}}},
      kep3::MU_SUN};
  std::array<double, 6> par3{{kep3::AU, 0., 0., 0., 0., 0.}};
  kep3::udpla::keplerian udpla3{ref_epoch, par3};
  std::array<double, 6> par4{{kep3::AU, 0., 0., 0., 0., 0.}};
  kep3::udpla::keplerian udpla4{ref_epoch, par4};
  double period_in_days =
      (2 * kep3::pi * kep3::AU / kep3::EARTH_VELOCITY) * kep3::SEC2DAY;
  auto [r, v] = udpla1.eph(ref_epoch + 200* period_in_days);
  fmt::print("r: {},{},{}", r[0] / kep3::AU, r[1] / kep3::AU, r[2] / kep3::AU);
  fmt::print("r: {},{},{}", v[0] / kep3::EARTH_VELOCITY,
             v[1] / kep3::EARTH_VELOCITY, v[2] / kep3::EARTH_VELOCITY);
}
