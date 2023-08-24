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

#include <fmt/core.h>

#include <kep3/detail/exceptions.hpp>
#include <kep3/epoch.hpp>
#include <kep3/planet.hpp>

#include "catch.hpp"

using kep3::epoch;
using kep3::planet;
using kep3::detail::is_udpla;
using kep3::detail::null_udpla;


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

struct complete_udpla {
  explicit complete_udpla(std::array<double, 4> physical_properties = {-1., -1.,
                                                                       -1.,
                                                                       -1.})
      : m_name("A complete, albeit simple Planet"),
        m_mu_central_body(physical_properties[0]),
        m_mu_self(physical_properties[1]), m_radius(physical_properties[2]),
        m_safe_radius(physical_properties[3]){};
  static std::array<std::array<double, 3>, 2> eph(const epoch &) {
    std::array<double, 3> pos = {1., 0., 0.};
    std::array<double, 3> vel = {0., 1., 0.};
    return {pos, vel};
  };

  [[nodiscard]] std::string get_name() const { return m_name; }
  [[nodiscard]] double get_mu_central_body() const { return m_mu_central_body; }
  [[nodiscard]] double get_mu_self() const { return m_mu_self; }
  [[nodiscard]] double get_radius() const { return m_radius; }
  [[nodiscard]] double get_safe_radius() const { return m_safe_radius; }

  [[nodiscard]] std::string get_extra_info() const {
    return fmt::format(
        "Gravitational parameter of main attracting body: {}\nGravitational "
        "parameter of body: {}\nBody radius: {}\nBody safe radius: {}\n",
        get_mu_central_body(), get_mu_self(), get_radius(), get_safe_radius());
  }

  std::string m_name;
  double m_mu_central_body;
  double m_mu_self;
  double m_radius;
  double m_safe_radius;
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
            kep3::detail::type_name<null_udpla>());
    REQUIRE(pla.get_extra_info() == std::string(""));
    double dummy = 0.;
    REQUIRE_THROWS_AS(dummy = pla.get_mu_central_body(),
                      kep3::not_implemented_error);
    REQUIRE_THROWS_AS(dummy = pla.get_mu_self(), kep3::not_implemented_error);
    REQUIRE_THROWS_AS(dummy = pla.get_radius(), kep3::not_implemented_error);
    REQUIRE_THROWS_AS(dummy = pla.get_safe_radius(),
                      kep3::not_implemented_error);
    REQUIRE(pla.extract<null_udpla>() != nullptr);
  }
  {
    // Constructor from a simple udpla
    simple_udpla udpla{};
    REQUIRE_NOTHROW(planet(udpla));
    planet pla(udpla);
    auto pos_vel = pla.eph(epoch(0.));
    REQUIRE(pos_vel[0] == std::array<double, 3>{1., 0., 0.});
    REQUIRE(pos_vel[1] == std::array<double, 3>{0., 1., 0.});
    REQUIRE(pla.get_name() == "A simple planet");
    REQUIRE(pla.get_extra_info() == "The simplest planet ever!");
    double dummy = 0.;
    REQUIRE_THROWS_AS(dummy = pla.get_mu_central_body(),
                      kep3::not_implemented_error);
    REQUIRE_THROWS_AS(dummy = pla.get_mu_self(), kep3::not_implemented_error);
    REQUIRE_THROWS_AS(dummy = pla.get_radius(), kep3::not_implemented_error);
    REQUIRE_THROWS_AS(dummy = pla.get_safe_radius(),
                      kep3::not_implemented_error);
  }
  {
    // Constructor from a more complete udpla
    complete_udpla udpla({1., 2., -1., 4.});
    REQUIRE_NOTHROW(planet(udpla));
    planet pla(udpla);
    auto pos_vel = pla.eph(epoch(0.));
    REQUIRE(pos_vel[0] == std::array<double, 3>{1., 0., 0.});
    REQUIRE(pos_vel[1] == std::array<double, 3>{0., 1., 0.});
    REQUIRE(pla.get_name() == "A complete, albeit simple Planet");
    REQUIRE(pla.get_mu_central_body() == 1.);
    REQUIRE(pla.get_mu_self() == 2.);
    double dummy = 0.;
    REQUIRE_THROWS_AS(dummy = pla.get_radius(), kep3::not_implemented_error);
    REQUIRE(pla.get_safe_radius() == 4.);
  }
  // Check copy semantics.
  auto p0 = planet();
  planet p1{p0};
  REQUIRE(p0.extract<null_udpla>() != nullptr);
  REQUIRE(p1.extract<null_udpla>() != nullptr);
  planet p2{simple_udpla{}};
  p2 = p1;
  REQUIRE(p2.extract<null_udpla>() != nullptr);
  // Move semantics.
  planet p3{std::move(p0)};
  REQUIRE(p3.extract<null_udpla>() != nullptr);
  planet p4{simple_udpla{}};
  p4 = std::move(p2);
  REQUIRE((p4.extract<null_udpla>() != nullptr));
  // Check we can revive moved-from objects.
  p0 = p4;
  REQUIRE((p0.extract<null_udpla>() != nullptr));
  p2 = std::move(p4);
  REQUIRE((p2.extract<null_udpla>() != nullptr));

  // Check the is_udpla type trait.
  REQUIRE(is_udpla<simple_udpla>::value);
  REQUIRE(is_udpla<null_udpla>::value);
  REQUIRE(!is_udpla<simple_udpla &>::value);
  REQUIRE(!is_udpla<const simple_udpla>::value);
  REQUIRE(!is_udpla<int>::value);
  REQUIRE(!is_udpla<void>::value);
  REQUIRE(!is_udpla<std::string>::value);
  REQUIRE(std::is_constructible<planet, simple_udpla>::value);
  REQUIRE(std::is_constructible<planet, null_udpla>::value);
  REQUIRE(std::is_constructible<planet, simple_udpla &>::value);
  REQUIRE(
      std::is_constructible<planet, const null_udpla &>::value);
  REQUIRE(std::is_constructible<planet, simple_udpla &&>::value);
  REQUIRE(!std::is_constructible<planet, int>::value);
  REQUIRE(!std::is_constructible<planet, std::string>::value);
}
