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
#include <kep3/detail/exceptions.hpp>
#include <kep3/planets/keplerian.hpp>
#include <stdexcept>

#include "catch.hpp"
#include "test_helpers.hpp"

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
  // From parameters kep3::elements_type::KEP_F
  std::array<double, 6> par0{{1., 0., 0., 0., 0., 0.}};
  REQUIRE_NOTHROW(kep3::udpla::keplerian{ref_epoch, par0, 1., "unknown"});
  REQUIRE_NOTHROW(
      kep3::udpla::keplerian{ref_epoch, par0, 1., "unknown", {-1, -1, -1}});
  // Checking the data members initializations:
  kep3::udpla::keplerian udpla{
      ref_epoch, par0, 1.1, "unknown", {1.2, 2.2, 1.9}};
  REQUIRE(udpla.get_ref_epoch() == ref_epoch);
  REQUIRE(udpla.get_name() == "unknown");
  REQUIRE(udpla.get_mu_central_body() == 1.1);
  REQUIRE(udpla.get_mu_self() == 1.2);
  REQUIRE(udpla.get_radius() == 2.2);
  REQUIRE(udpla.get_safe_radius() == 1.9);
  REQUIRE(udpla.period() == 2 * kep3::pi * std::sqrt(1. / 1.1));
  // Calling constructor with different elements type
  {
    std::array<double, 6> par{{1., 0., 0., 0., 0., 0.}};
    REQUIRE_NOTHROW(kep3::udpla::keplerian{ref_epoch,
                                           par,
                                           1.,
                                           "unknown",
                                           {-1, -1, -1},
                                           kep3::elements_type::KEP_F});
  }
  {
    std::array<double, 6> par{{1., 0., 0., 0., 0., 0.}};
    REQUIRE_NOTHROW(kep3::udpla::keplerian{ref_epoch,
                                           par,
                                           1.,
                                           "unknown",
                                           {-1, -1, -1},
                                           kep3::elements_type::KEP_M});
  }
  {
    std::array<double, 6> par{{1., 0., 0., 1., 0., 0.}};
    REQUIRE_NOTHROW(kep3::udpla::keplerian{
        ref_epoch, par, 1., "unknown", {-1, -1, -1}, kep3::elements_type::MEQ});
  }
  {
    std::array<double, 6> par{{1., 0., 0., 1., 0., 0.}};
    REQUIRE_NOTHROW(kep3::udpla::keplerian{ref_epoch,
                                           par,
                                           1.,
                                           "unknown",
                                           {-1, -1, -1},
                                           kep3::elements_type::MEQ_R});
  }
  { // hyperbola and mean anomaly????
    std::array<double, 6> par{{-10., 10., 0., 1., 0., 0.}};
    REQUIRE_THROWS_AS((kep3::udpla::keplerian{ref_epoch,
                                              par,
                                              1.,
                                              "unknown",
                                              {-1, -1, -1},
                                              kep3::elements_type::KEP_M}),
                      std::logic_error);
  }
  { // posvel as 1x6 orbital parameters????
    std::array<double, 6> par{{1., 0., 0., 1., 0., 0.}};
    REQUIRE_THROWS_AS((kep3::udpla::keplerian{ref_epoch,
                                              par,
                                              1.,
                                              "unknown",
                                              {-1, -1, -1},
                                              kep3::elements_type::POSVEL}),
                      std::logic_error);
  }
  { // negative a but ecc < 1????
    std::array<double, 6> par{{-10., 0., 0., 1., 0., 0.}};
    REQUIRE_THROWS_AS((kep3::udpla::keplerian{ref_epoch,
                                              par,
                                              1.,
                                              "unknown",
                                              {-1, -1, -1},
                                              kep3::elements_type::KEP_F}),
                      std::domain_error);
  }
  { // We construct an hyperbolic planet
    std::array<double, 6> par{{-10., 10., 0., 0., 0., 0.}};
    kep3::udpla::keplerian udpla2{ref_epoch,    par,
                                 1.,           "unknown",
                                 {-1, -1, -1}, kep3::elements_type::KEP_F};
    REQUIRE(!std::isfinite(udpla2.period()));
  }
  { // We construct an hyperbolic planet
    std::array<std::array<double, 3>, 2> posvel{{{1,0,0},{0,10,0}}};
    kep3::udpla::keplerian udpla2{ref_epoch,    posvel,
                                 1.,           "unknown",
                                 {-1, -1, -1}};
    REQUIRE(!std::isfinite(udpla2.period()));
  }
}

TEST_CASE("eph") {
  // We use 2000-01-01 as a reference epoch for all these tests
  kep3::epoch ref_epoch{0., kep3::epoch::MJD2000};
  // This is a circular orbit at 1 AU.
  std::array<std::array<double, 3>, 2> pos_vel_0{
      {{kep3::AU, 0., 0.}, {0., kep3::EARTH_VELOCITY, 0.}}};
  // A keplerian planet orbiting the Sun on such a perfectly circular orbit.
  kep3::udpla::keplerian udpla1{ref_epoch, pos_vel_0, kep3::MU_SUN};
  double period_in_days =
      (2 * kep3::pi * kep3::AU / kep3::EARTH_VELOCITY) * kep3::SEC2DAY;
  auto [r, v] = udpla1.eph(ref_epoch + period_in_days);
  REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel_0[0]) < 1e-13);
  REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel_0[1]) < 1e-13);
}
