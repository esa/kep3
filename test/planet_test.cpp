// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <kep3/planet.hpp>
#include <kep3/epoch.hpp>


#include "catch.hpp"

using kep3::planet;
using kep3::epoch;


TEST_CASE("construction") {
  {
    // Default constructor (a null planet)
    REQUIRE_NOTHROW(planet());
    auto pla = planet();
    REQUIRE_NOTHROW(pla.eph(epoch(0.)));
    auto pos_vel = pla.eph(epoch(0.));
    REQUIRE(pos_vel[0] == std::array<double,3>{1.,0.,0.});
    REQUIRE(pos_vel[1] == std::array<double,3>{0.,1.,0.});
    REQUIRE(pla.get_name() == kep3::detail::type_name<kep3::detail::null_udpla>());
    REQUIRE(pla.get_extra_info() == std::string(""));
  }
}
