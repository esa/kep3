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
  REQUIRE(kep3_tests::floating_point_error(earth.period() * kep3::SEC2DAY, 365.25) < 0.01);
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