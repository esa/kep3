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
#include <stdexcept>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/core_astro/ic2eq2ic.hpp>
#include <kep3/exceptions.hpp>
#include <kep3/planets/jpl_lp.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using kep3::udpla::jpl_lp;

TEST_CASE("constructor") {
  REQUIRE_NOTHROW(jpl_lp{});
  kep3::epoch ref_epoch{12.22, kep3::epoch::MJD2000};
  // From name
  REQUIRE_NOTHROW(jpl_lp{"Mars"});
  REQUIRE_NOTHROW(jpl_lp{"mars"});
  REQUIRE_NOTHROW(kep3::planet{jpl_lp{"neptune"}});

  REQUIRE_THROWS_AS(jpl_lp{"gigi"}, std::logic_error);
  jpl_lp udpla{"earth"};
  kep3::planet earth{udpla};
  REQUIRE(kep3_tests::floating_point_error(earth.period() * kep3::SEC2DAY,
                                           365.25) < 0.01);
}

TEST_CASE("eph") {
  // We use 2030-01-01 as a reference epoch for all these tests
  kep3::epoch ref_epoch{2458849.5, kep3::epoch::JD};
  {
    // This is the Earth-Moon w.r.t. the Sun queried from JPL Horizon at
    // 2020-01-01
    std::array<std::array<double, 3>, 2> pos_vel_0{
        {{-2.488023054631234E+10, 1.449771522542222E+11,
          -6.590293144971132E+02},
         {-2.984589828430694E+04, -5.151004951052294E+03,
          3.108878527788850E-01}}};
    // The Earth in jpl_lp mode
    jpl_lp udpla{"earth"};
    auto [r, v] = udpla.eph(ref_epoch);
    REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel_0[0]) < 0.01);
    REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel_0[1]) < 0.01);
  }
  {
    // This is Mercury w.r.t. the Sun queried from JPL Horizon at
    // 2020-01-01
    std::array<std::array<double, 3>, 2> pos_vel_0{
        {{-9.474762662376745E+09, -6.894147965135109E+10,
          -4.764334842347469E+09},
         {3.848711305256677E+04, -4.155242103836629E+03,
          -3.870162659830893E+03}}};
    // Mercury in jpl_lp mode
    jpl_lp udpla{"mercury"};
    auto [r, v] = udpla.eph(ref_epoch);
    REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel_0[0]) < 0.05);
    REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel_0[1]) < 0.05);
  }
  {
    // This is Venus w.r.t. the Sun queried from JPL Horizon at
    // 2020-01-01
    std::array<std::array<double, 3>, 2> pos_vel_0{
        {{1.081892249749067E+11, 7.861125522626230E+09, -6.135421905733132E+09},
         {-2.679023504807336E+03, 3.476995213020635E+04,
          6.316923820826013E+02}}};
    // Venus in jpl_lp mode
    jpl_lp udpla{"venus"};
    auto [r, v] = udpla.eph(ref_epoch);
    REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel_0[0]) < 0.02);
    REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel_0[1]) < 0.02);
  }
}

TEST_CASE("elements") {
  kep3::epoch ref_epoch{12.22, kep3::epoch::MJD2000};
  // We use Neptune
  jpl_lp udpla{"nePTUne"}; // casing is not important
  auto pos_vel = udpla.eph(ref_epoch);
  // Test on various element types
  {
    auto par = udpla.elements(ref_epoch, kep3::elements_type::KEP_F);
    auto [r, v] = kep3::par2ic(par, kep3::MU_SUN);
    REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel[0]) < 1e-13);
    REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel[1]) < 1e-13);
  }
  {
    auto par = udpla.elements(ref_epoch, kep3::elements_type::KEP_M);
    par[5] = kep3::m2f(par[5], par[1]);
    auto [r, v] = kep3::par2ic(par, kep3::MU_SUN);
    REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel[0]) < 1e-13);
    REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel[1]) < 1e-13);
  }
  {
    auto par = udpla.elements(ref_epoch, kep3::elements_type::MEQ);
    auto [r, v] = kep3::eq2ic(par, kep3::MU_SUN);
    REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel[0]) < 1e-13);
    REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel[1]) < 1e-13);
  }
  {
    auto par = udpla.elements(ref_epoch, kep3::elements_type::MEQ_R);
    auto [r, v] = kep3::eq2ic(par, kep3::MU_SUN, true);
    REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel[0]) < 1e-13);
    REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel[1]) < 1e-13);
  }
  {
    REQUIRE_THROWS_AS(udpla.elements(ref_epoch, kep3::elements_type::POSVEL),
                      std::logic_error);
  }
}

TEST_CASE("stream_operator") {
  REQUIRE_NOTHROW((std::cout << jpl_lp{} << '\n'));
}

TEST_CASE("serialization_test") {
  // Instantiate a generic jpl_lp
  kep3::epoch ref_epoch{2423.4343, kep3::epoch::MJD2000};
  jpl_lp udpla{"neptune"};

  // Store the string representation.
  std::stringstream ss;
  auto before = boost::lexical_cast<std::string>(udpla);
  // Now serialize
  {
    boost::archive::binary_oarchive oarchive(ss);
    oarchive << udpla;
  }
  // Deserialize
  // Create a new udpla object
  jpl_lp udpla2{};
  {
    boost::archive::binary_iarchive iarchive(ss);
    iarchive >> udpla2;
  }
  auto after = boost::lexical_cast<std::string>(udpla2);
  // Compare the string represetation
  REQUIRE(before == after);
}