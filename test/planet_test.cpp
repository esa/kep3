// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <string>

#include <kep3/epoch.hpp>
#include <kep3/planet.hpp>

#include "catch.hpp"

using kep3::epoch;
using kep3::planet;

struct simple_udpla {
  simple_udpla() = default;
  static std::array<std::array<double, 3>, 2> eph(const epoch &) {
    std::array<double, 3> pos = {1., 0., 0.};
    std::array<double, 3> vel = {0., 1., 0.};
    return {pos, vel};
  };
  static std::string get_name() { return "A simple planet"; }
  static std::string get_extra_info() { return "The simplest planet ever!"; }
};

TEST_CASE("construction") {
  {
    // Default constructor (a null planet)
    REQUIRE_NOTHROW(planet());
    auto pla = planet();
    REQUIRE_NOTHROW(pla.eph(epoch(0.)));
    auto pos_vel = pla.eph(epoch(0.));
    REQUIRE(pos_vel[0] == std::array<double, 3>{1., 0., 0.});
    REQUIRE(pos_vel[1] == std::array<double, 3>{0., 1., 0.});
    REQUIRE(pla.get_name() ==
            kep3::detail::type_name<kep3::detail::null_udpla>());
    REQUIRE(pla.get_extra_info() == std::string(""));
  }
  {
    // Constructor from udpla
    simple_udpla udpla{};
    REQUIRE_NOTHROW(planet(udpla));
    planet pla(udpla);
    auto pos_vel = pla.eph(epoch(0.));
    REQUIRE(pos_vel[0] == std::array<double, 3>{1., 0., 0.});
    REQUIRE(pos_vel[1] == std::array<double, 3>{0., 1., 0.});
    REQUIRE(pla.get_name() == "A simple udpla");
    REQUIRE(pla.get_extra_info() == "The simplest udpla ever!");
  }
}
