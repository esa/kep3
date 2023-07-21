// Copyright 2023 Francesco Biscani (bluescarni@gmail.com), Dario Izzo
// (dario.izzo@gmail.com)
//
// This file is part of the cascade library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/date_time/gregorian/gregorian.hpp>
#include <iostream>
#include <random>

#include <epoch.hpp>

#include "catch.hpp"

using namespace kep3;
using namespace boost::gregorian;
using namespace boost::posix_time;

TEST_CASE("construct") {
  // test syntax
  REQUIRE_NOTHROW(epoch());
  REQUIRE_NOTHROW(epoch(123.456));
  REQUIRE_NOTHROW(epoch(123.456, epoch::MJD2000));
  REQUIRE_NOTHROW(epoch(123.456, epoch::JD));
  REQUIRE_NOTHROW(epoch(123.456, epoch::MJD));
  REQUIRE_NOTHROW(epoch(31, 12, 2034));
  boost::posix_time::ptime posix_time_test(date(2034, 12, 31));
  REQUIRE_NOTHROW(epoch(posix_time_test));
  // test conversions
  REQUIRE(epoch(123.456).mjd2000() == epoch(123.456, epoch::MJD2000).mjd2000());
  REQUIRE(epoch(0.).mjd2000() == epoch(51544, epoch::MJD).mjd2000());
  REQUIRE(epoch(0.).mjd2000() == epoch(2451544.5, epoch::JD).mjd2000());
  REQUIRE(epoch(31, 12, 2034).jd() == epoch(posix_time_test).jd());
}

TEST_CASE("epoch_operators") {
  REQUIRE(epoch(0.) == epoch(0.));
  REQUIRE(epoch(0.) != epoch(1.));
  REQUIRE(epoch(1.) > epoch(0.));
  REQUIRE(epoch(1.) >= epoch(1.));
  REQUIRE(epoch(1.) >= epoch(0.));
  REQUIRE(epoch(0.) < epoch(1.));
  REQUIRE(epoch(1.) <= epoch(1.));
  epoch today(0.);
  today += 100.;
  REQUIRE(today == epoch(100.));
  today -= 100.;
  REQUIRE(today == epoch(0.));
  auto yesterday = today - 1.;
  REQUIRE(yesterday == epoch(-1.));
  today = yesterday + 1;
  REQUIRE(today == epoch(0.));
  REQUIRE(today - yesterday == 1.);
}
