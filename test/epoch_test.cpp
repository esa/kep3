// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <chrono>
#include <ctime>
#include <iostream>
#include <random>
#include <string>

#include <kep3/epoch.hpp>

#include "catch.hpp"

using kep3::epoch;
using kep3::kep_clock;
using namespace std::literals;
namespace chr = std::chrono;

TEST_CASE("construct") {
  // test syntax

  // // > 2000
  REQUIRE_NOTHROW(epoch());
  REQUIRE_NOTHROW(epoch(123.456));
  REQUIRE_NOTHROW(epoch(123.456, epoch::julian_type::MJD2000));
  REQUIRE_NOTHROW(epoch(0.0, epoch::julian_type::JD));
  REQUIRE_NOTHROW(epoch(123.456, epoch::julian_type::JD));
  REQUIRE_NOTHROW(epoch(0.0, epoch::julian_type::MJD));
  REQUIRE_NOTHROW(epoch(123.456, epoch::julian_type::MJD));
  REQUIRE_NOTHROW(epoch(2000, 10, 17, 11, 36, 21, 121, 841));
  REQUIRE_NOTHROW(epoch(2064, 10, 17, 11, 36, 21, 121, 841).as_utc_string() == "2064-10-17T11:36:21.121841");

  // // < 2000
  REQUIRE_NOTHROW(epoch(-123.456));
  REQUIRE_NOTHROW(epoch(-123.456, epoch::julian_type::MJD2000));
  REQUIRE_NOTHROW(epoch(-0.0, epoch::julian_type::JD));
  REQUIRE_NOTHROW(epoch(-123.456, epoch::julian_type::JD));
  REQUIRE_NOTHROW(epoch(-0.0, epoch::julian_type::MJD));
  REQUIRE_NOTHROW(epoch(-123.456, epoch::julian_type::MJD));
  REQUIRE_NOTHROW(epoch(1980, 10, 17, 11, 36, 21, 121, 841));
  REQUIRE_NOTHROW(epoch(1980, 10, 17, 11, 36, 21, 121, 841).as_utc_string() == "1980-10-17T11:36:21.121841");

  // Epoch from lvalue and rvalue references
  epoch ep{2000, 1, 1};
  REQUIRE(epoch(ep) == ep);
  REQUIRE(epoch(epoch{2000,1,1}) == ep);

  // test conversions
  REQUIRE(epoch(123.456).mjd2000() ==
          epoch(123.456, epoch::julian_type::MJD2000).mjd2000());
  REQUIRE(epoch(0.).mjd() == epoch(51544., epoch::julian_type::MJD).mjd());
  REQUIRE(epoch(0.).jd() == epoch(2451544.5, epoch::julian_type::JD).jd());
}

TEST_CASE("epoch_operators") {
  REQUIRE(epoch(2034, 10, 17) == epoch(2034, 10, 17));
  REQUIRE(epoch(2034, 10, 17) != epoch(2034, 11, 17));
  // Testing us precision
  REQUIRE(epoch(2034, 10, 17) != epoch(2034, 10, 17, 0, 0, 0, 0, 1));
  // Check that ns precision is not supported
  REQUIRE(epoch(2000, 10, 17) == epoch(2000, 10, 17, 0, 0, 0, 0, 0) + chr::nanoseconds(100));

  // Conversion from double (defaults to days)
  REQUIRE(epoch(1.) > epoch(0.));
  REQUIRE(epoch(1.) >= epoch(1.));
  REQUIRE(epoch(1.) >= epoch(0.));
  REQUIRE(epoch(0.) < epoch(1.));
  REQUIRE(epoch(1.) <= epoch(1.));
  epoch today(0.);
  auto offset{chr::days(10963)};
  today += offset;
  std::cout << "TODAY: " << today << "\n";
  REQUIRE(today == epoch(2030, 1, 6));
  today -= chr::duration_cast<kep_clock::duration>(offset);
  REQUIRE(today == epoch());
  auto oneday{chr::days(1)};
  auto yesterday{today - chr::duration_cast<kep_clock::duration>(oneday)};
  auto yesterday1{today - oneday};
  REQUIRE(yesterday == yesterday1);

  REQUIRE(yesterday == epoch(1999, 12, 31));
  today = yesterday + chr::duration_cast<kep_clock::duration>(chr::days(1));
  REQUIRE(today == epoch());
  auto diff{today - yesterday};
  REQUIRE(diff == chr::duration_cast<kep_clock::duration>(chr::days(1)));
  REQUIRE_NOTHROW((std::cout << epoch()));
}
