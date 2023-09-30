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
using kep3::kep_clock;
using namespace std::literals;
namespace chr = std::chrono;

TEST_CASE( "construct" )
{
    // test syntax
    // // > 2000
    REQUIRE_NOTHROW( epoch() );
    REQUIRE_NOTHROW( epoch( 123.456 ) );
    REQUIRE_NOTHROW( epoch( 123.456, epoch::julian_type::MJD2000 ) );
    REQUIRE_NOTHROW( epoch( 0.0, epoch::julian_type::JD ) );
    REQUIRE_NOTHROW( epoch( 123.456, epoch::julian_type::JD ) );
    REQUIRE_NOTHROW( epoch( 0.0, epoch::julian_type::MJD ) );
    REQUIRE_NOTHROW( epoch( 123.456, epoch::julian_type::MJD ) );
    REQUIRE_NOTHROW( epoch( 34, 364, 11, 36, 21, 121, 841 ) );
    // // > 2000
    // boost::posix_time::ptime posix_time_test( date( 2034, 12, 31 ) );
    // REQUIRE_NOTHROW( epoch( posix_time_test ) );
    // REQUIRE( epoch( posix_time_test ).get_posix_time() == posix_time_test );
    //   std::cout << "Chrono1: " << epoch(posix_time_test).get_posix_time() <<
    //   "\n"; std::cout << "Boost1: " << posix_time_test << "\n";

    // // < 2000
    REQUIRE_NOTHROW( epoch( -123.456 ) );
    REQUIRE_NOTHROW( epoch( -123.456, epoch::julian_type::MJD2000 ) );
    REQUIRE_NOTHROW( epoch( -0.0, epoch::julian_type::JD ) );
    REQUIRE_NOTHROW( epoch( -123.456, epoch::julian_type::JD ) );
    REQUIRE_NOTHROW( epoch( -0.0, epoch::julian_type::MJD ) );
    REQUIRE_NOTHROW( epoch( -123.456, epoch::julian_type::MJD ) );
    REQUIRE_NOTHROW( epoch( -34, 364, 11, 36, 21, 121, 841 ) );
    //   boost::posix_time::ptime posix_time_test2(date(1980, 12, 31));
    //   REQUIRE_NOTHROW(epoch(posix_time_test2));
    //   REQUIRE(epoch(posix_time_test2).get_posix_time() == posix_time_test2);
    //   std::cout << "Chrono2: " << epoch(posix_time_test2).get_posix_time() <<
    //   "\n"; std::cout << "Boost2: " << posix_time_test2 << "\n";

    // // test conversions
      REQUIRE(epoch(123.456).mjd2000() ==
              epoch(123.456, epoch::julian_type::MJD2000).mjd2000());
      REQUIRE(epoch(0.).mjd() == epoch(51544, epoch::julian_type::MJD).mjd());
      REQUIRE(epoch(0.).jd() == epoch(2451544.5, epoch::julian_type::JD).jd());
    //   REQUIRE(epoch(2034, 12, 31) == epoch(posix_time_test));
}

TEST_CASE( "epoch_operators" )
{
    REQUIRE( epoch( 34, 10 ) == epoch( 34, 10 ) );
    REQUIRE( epoch( 34, 10 ) != epoch( 34, 11 ) );
    // Testing us precision
    REQUIRE( epoch( 34,  10 ) != epoch( 34, 10, 0, 0, 0, 0, 1 ) );
    // Check that ns precision is not supported
    REQUIRE( epoch( 0, 10 ) ==
             epoch( 0, 10, 0, 0, 0, 0, 0 ) + chr::nanoseconds(100) );

    //    Conversion from double (defaults to days)
    REQUIRE( epoch( 1. ) > epoch( 0. ) );
    REQUIRE( epoch( 1. ) >= epoch( 1. ) );
    REQUIRE( epoch( 1. ) >= epoch( 0. ) );
    REQUIRE( epoch( 0. ) < epoch( 1. ) );
    REQUIRE( epoch( 1. ) <= epoch( 1. ) );
    epoch today( 0. );
    today += chr::days( 100 );
    REQUIRE( today == epoch( 0, 100 ) );
    today -= chr::duration_cast<kep_clock::duration>( chr::days( 100 ) );
    REQUIRE( today == epoch() );
    auto oneday{chr::days( 1 )};
    auto yesterday{today - chr::duration_cast<kep_clock::duration>( oneday )};
    auto yesterday1{today - oneday};
    REQUIRE(yesterday == yesterday1);

    REQUIRE( yesterday == epoch( 0, -1 ) );
    today = yesterday + chr::duration_cast<kep_clock::duration>( chr::days( 1 ) );
    REQUIRE( today == epoch() );
    auto diff{ today - yesterday };
    REQUIRE( diff == chr::duration_cast<kep_clock::duration>( chr::days( 1 ) ) );
    REQUIRE_NOTHROW( ( std::cout << epoch() ) );
}

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
