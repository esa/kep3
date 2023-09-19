// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/date_time/gregorian/gregorian.hpp>
#include <chrono>
#include <ctime>
#include <iostream>
#include <random>
#include <string>

#include "kep3/epoch.hpp"

#include "catch.hpp"

using kep3::epoch;
// using kep3::epoch_from_iso_string;
// using kep3::epoch_from_string;
using boost::gregorian::date;
using boost::posix_time::ptime;

using namespace std::literals;

TEST_CASE("construct")
{
  // test syntax
  REQUIRE_NOTHROW(epoch());
  REQUIRE_NOTHROW(epoch(123.456));
  REQUIRE_NOTHROW(epoch(123.456, epoch::julian_type::MJD2000));
  REQUIRE_NOTHROW(epoch(0.0, epoch::julian_type::JD));
  REQUIRE_NOTHROW(epoch(123.456, epoch::julian_type::JD));
  REQUIRE_NOTHROW(epoch(0.0, epoch::julian_type::MJD));
  REQUIRE_NOTHROW(epoch(123.456, epoch::julian_type::MJD));
  REQUIRE_NOTHROW(epoch(2034, 12, 31));
  // // > 2000
  boost::posix_time::ptime posix_time_test(date(2034, 12, 31));
  REQUIRE_NOTHROW(epoch(posix_time_test));
  REQUIRE(epoch(posix_time_test).get_posix_time() == posix_time_test);
  std::cout << "Chrono1: " << epoch(posix_time_test).get_posix_time() << "\n";
  std::cout << "Boost1: " << posix_time_test << "\n";

  // // < 2000
  boost::posix_time::ptime posix_time_test2(date(1980, 12, 31));
  REQUIRE_NOTHROW(epoch(posix_time_test2));
  REQUIRE(epoch(posix_time_test2).get_posix_time() == posix_time_test2);
  std::cout << "Chrono2: " << epoch(posix_time_test2).get_posix_time() << "\n";
  std::cout << "Boost2: " << posix_time_test2 << "\n";

  // // test conversions
  REQUIRE(epoch(123.456).mjd2000() == epoch(123.456, epoch::julian_type::MJD2000).mjd2000());
  REQUIRE(epoch(0.).mjd() == epoch(51544, epoch::julian_type::MJD).mjd());
  REQUIRE(epoch(0.).jd() == epoch(2451544.5, epoch::julian_type::JD).jd());
  std::cout << "Chrono3: " << kep3::epoch::as_utc_string(epoch(posix_time_test).jd()) << "\n";
  std::cout << "Chrono3-1: " << kep3::epoch::as_utc_string(epoch(2034, 12, 31).jd()) << "\n";
  REQUIRE(epoch(2034, 12, 31).jd().time_since_epoch().count() == epoch(posix_time_test).jd().time_since_epoch().count());

}

// TEST_CASE("epoch_operators")
// {
//   REQUIRE(epoch(0.) == epoch(0.));
//   REQUIRE(epoch(0.) != epoch(1.));
//   REQUIRE(epoch(1.) > epoch(0.));
//   REQUIRE(epoch(1.) >= epoch(1.));
//   REQUIRE(epoch(1.) >= epoch(0.));
//   REQUIRE(epoch(0.) < epoch(1.));
//   REQUIRE(epoch(1.) <= epoch(1.));
//   epoch today(0.);
//   today += 100.;
//   REQUIRE(today == epoch(100.));
//   today -= 100.;
//   REQUIRE(today == epoch(0.));
//   auto yesterday = today - 1.;
//   REQUIRE(yesterday == epoch(-1.));
//   today = yesterday + 1;
//   REQUIRE(today == epoch(0.));
//   REQUIRE(today - yesterday == 1.);
//   REQUIRE_NOTHROW((std::cout << epoch()));
// }

// TEST_CASE("conversions_from_string") {
//   {
//     std::string ts("20020131T000000");
//     epoch e(epoch_from_iso_string(ts));
//     REQUIRE(e == epoch(boost::posix_time::ptime(date(2002, 01, 31))));
//   }
//   {
//     std::string ts("2002-01-20 00:00:00.000");
//     epoch e(epoch_from_string(ts));
//     REQUIRE(e == epoch(boost::posix_time::ptime(date(2002, 01, 20))));
//   }
// }
