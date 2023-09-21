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

#include <kep3/epoch.hpp>
#include <kep3/exceptions.hpp>
#include <kep3/planet.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using kep3::epoch;
using kep3::planet;
using kep3::detail::is_udpla;
using kep3::detail::null_udpla;
using kep3::detail::udpla_has_eph;
using kep3::detail::udpla_has_get_extra_info_v;
using kep3::detail::udpla_has_get_mu_central_body_v;
using kep3::detail::udpla_has_get_mu_self_v;
using kep3::detail::udpla_has_get_name_v;
using kep3::detail::udpla_has_get_radius_v;
using kep3::detail::udpla_has_get_safe_radius_v;
using kep3::detail::udpla_has_period_v;

struct simple_udpla {
  simple_udpla() = default;
  static std::array<std::array<double, 3>, 2> eph(const epoch &) {
    std::array<double, 3> pos = {1., 0., 0.};
    std::array<double, 3> vel = {0., 1., 0.};
    return {pos, vel};
  };
  static std::string get_name() { return "A simple planet"; }
  static std::string get_extra_info() { return "The simplest planet ever!"; }

private:
  friend class boost::serialization::access;
  template <typename Archive> void serialize(Archive &, unsigned) {}
};
kep3_S11N_PLANET_EXPORT(simple_udpla);

struct simple_udpla_mu {
  simple_udpla_mu() = default;
  static std::array<std::array<double, 3>, 2> eph(const epoch &) {
    std::array<double, 3> pos = {1., 0., 0.};
    std::array<double, 3> vel = {0., 1., 0.};
    return {pos, vel};
  };
  static std::string get_name() { return "A simple planet with mu"; }
  static std::string get_extra_info() {
    return "The simplest planet ever but with mu!";
  }
  static double get_mu_central_body() { return 1.; }

private:
  friend class boost::serialization::access;
  template <typename Archive> void serialize(Archive &, unsigned) {}
};
kep3_S11N_PLANET_EXPORT(simple_udpla_mu);

struct simple_udpla_mu_h {
  simple_udpla_mu_h() = default;
  static std::array<std::array<double, 3>, 2> eph(const epoch &) {
    std::array<double, 3> pos = {1., 0., 0.};
    std::array<double, 3> vel = {0., 10., 0.};
    return {pos, vel};
  };
  static std::string get_name() { return "A simple planet with mu"; }
  static std::string get_extra_info() {
    return "The simplest planet ever but with mu!";
  }
  static double get_mu_central_body() { return 1.; }

private:
  friend class boost::serialization::access;
  template <typename Archive> void serialize(Archive &, unsigned) {}
};
kep3_S11N_PLANET_EXPORT(simple_udpla_mu_h);

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
  [[nodiscard]] double period(const epoch &) const {
    return m_radius - m_radius;
  }

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

private:
  friend class boost::serialization::access;
  template <typename Archive> void serialize(Archive &ar, unsigned) {
    ar &m_name;
    ar &m_mu_central_body;
    ar &m_mu_self;
    ar &m_radius;
    ar &m_safe_radius;
  }
};
kep3_S11N_PLANET_EXPORT(complete_udpla);

TEST_CASE("construction") {
  {
    kep3::detail::null_udpla udpla{};
    udpla.eph(epoch(0.));
    // Default constructor (a null planet)
    REQUIRE_NOTHROW(planet());
    auto pla = planet();
    REQUIRE_NOTHROW(pla.eph(epoch(0.)));
    auto pos_vel = pla.eph(epoch(0.));
    REQUIRE(pos_vel[0] == std::array<double, 3>{1., 0., 0.});
    REQUIRE(pos_vel[1] == std::array<double, 3>{0., 1., 0.});
    REQUIRE(pla.get_name() == kep3::detail::type_name<null_udpla>());
    REQUIRE(pla.get_extra_info() == std::string(""));
    REQUIRE_THROWS_AS((pla.get_mu_central_body()), kep3::not_implemented_error);
    REQUIRE_THROWS_AS((pla.get_mu_self()), kep3::not_implemented_error);
    REQUIRE_THROWS_AS((pla.get_radius()), kep3::not_implemented_error);
    REQUIRE_THROWS_AS((pla.get_safe_radius()), kep3::not_implemented_error);
    REQUIRE_THROWS_AS((pla.period()), kep3::not_implemented_error);
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
    REQUIRE_THROWS_AS((pla.get_mu_central_body()), kep3::not_implemented_error);
    REQUIRE_THROWS_AS((pla.get_mu_self()), kep3::not_implemented_error);
    REQUIRE_THROWS_AS((pla.get_radius()), kep3::not_implemented_error);
    REQUIRE_THROWS_AS((pla.get_safe_radius()), kep3::not_implemented_error);
    REQUIRE_THROWS_AS((pla.period()), kep3::not_implemented_error);
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
    REQUIRE_THROWS_AS((pla.get_radius()), kep3::not_implemented_error);
    REQUIRE(pla.get_safe_radius() == 4.);
  }
  {
    // Constructor from a simple udpla with mu and hyperbolic orbit
    simple_udpla_mu_h udpla{};
    REQUIRE_NOTHROW(planet(udpla));
    planet pla(udpla);
    REQUIRE(!std::isfinite(pla.period()));
  }
  {
    // Constructor from a simple udpla with mu and elliptical orbit
    simple_udpla_mu udpla{};
    REQUIRE_NOTHROW(planet(udpla));
    planet pla(udpla);
    REQUIRE(kep3_tests::floating_point_error(pla.period(), kep3::pi * 2.) <
            1e-14);
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
    REQUIRE_THROWS_AS((pla.get_radius()), kep3::not_implemented_error);
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
}

TEST_CASE("type_traits") {
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
  REQUIRE(std::is_constructible<planet, const null_udpla &>::value);
  REQUIRE(std::is_constructible<planet, simple_udpla &&>::value);
  REQUIRE(!std::is_constructible<planet, int>::value);
  REQUIRE(!std::is_constructible<planet, std::string>::value);
  // check the udpla_has_eph type trait.
  REQUIRE(udpla_has_eph<null_udpla>::value);
  REQUIRE(udpla_has_eph<simple_udpla>::value);
  REQUIRE(!udpla_has_eph<int>::value);
  // check the udpla_has_get_name_v type trait.
  REQUIRE(!udpla_has_get_name_v<null_udpla>);
  REQUIRE(udpla_has_get_name_v<simple_udpla>);
  REQUIRE(!udpla_has_get_name_v<double>);
  // check the udpla_has_get_extra_info_v type trait.
  REQUIRE(!udpla_has_get_extra_info_v<null_udpla>);
  REQUIRE(udpla_has_get_extra_info_v<simple_udpla>);
  REQUIRE(!udpla_has_get_extra_info_v<double>);
  // check the udpla_has_get_mu_central_body_v type trait.
  REQUIRE(!udpla_has_get_mu_central_body_v<null_udpla>);
  REQUIRE(!udpla_has_get_mu_central_body_v<simple_udpla>);
  REQUIRE(!udpla_has_get_mu_central_body_v<double>);
  REQUIRE(udpla_has_get_mu_central_body_v<complete_udpla>);
  // check the udpla_has_get_mu_self_v type trait.
  REQUIRE(!udpla_has_get_mu_self_v<null_udpla>);
  REQUIRE(!udpla_has_get_mu_self_v<simple_udpla>);
  REQUIRE(!udpla_has_get_mu_self_v<double>);
  REQUIRE(udpla_has_get_mu_self_v<complete_udpla>);
  // check the udpla_has_get_radius_v type trait.
  REQUIRE(!udpla_has_get_radius_v<null_udpla>);
  REQUIRE(!udpla_has_get_radius_v<simple_udpla>);
  REQUIRE(!udpla_has_get_radius_v<double>);
  REQUIRE(udpla_has_get_radius_v<complete_udpla>);
  // check the udpla_has_get_safe_radius_v type trait.
  REQUIRE(!udpla_has_get_safe_radius_v<null_udpla>);
  REQUIRE(!udpla_has_get_safe_radius_v<simple_udpla>);
  REQUIRE(!udpla_has_get_safe_radius_v<double>);
  REQUIRE(udpla_has_get_safe_radius_v<complete_udpla>);
  // check the udpla_has_period_v type trait.
  REQUIRE(!udpla_has_period_v<null_udpla>);
  REQUIRE(!udpla_has_period_v<simple_udpla>);
  REQUIRE(!udpla_has_period_v<double>);
  REQUIRE(udpla_has_period_v<complete_udpla>);
}

TEST_CASE("copy_constructor_test") {
  // We instantiate a planet
  complete_udpla udpla({1., 2., -1., 4.});
  planet pla{udpla};

  // We call the copy constructor
  planet pla_copy(pla);
  // We extract the user planet
  auto *p1 = pla.extract<complete_udpla>();
  auto *p2 = pla_copy.extract<complete_udpla>();

  // 1 - We check the resources pointed to by m_ptr have different address
  REQUIRE(p1 != 0);
  REQUIRE(p2 != 0);
  REQUIRE(p1 != p2);
  // 2 - We check that the other members are copied
  REQUIRE(pla.get_name() == pla_copy.get_name());
  REQUIRE(pla.get_mu_central_body() == pla_copy.get_mu_central_body());
  REQUIRE(pla.get_mu_self() == pla_copy.get_mu_self());
  REQUIRE(pla.get_safe_radius() == pla_copy.get_safe_radius());
}

TEST_CASE("planet_move_constructor_test") {
  // We instantiate a planet
  complete_udpla udpla({1., 2., -1., 4.});
  planet pla{udpla};

  // We store a streaming representation of the object
  auto pla_string = boost::lexical_cast<std::string>(pla);
  // We get the memory address where the user algo is stored
  auto *p1 = pla.extract<complete_udpla>();
  // We call the move constructor
  planet moved_pla(std::move(pla));
  // We get the memory address where the user algo is stored
  auto *p2 = moved_pla.extract<complete_udpla>();
  // And the string representation of the moved algo
  auto moved_pla_string = boost::lexical_cast<std::string>(moved_pla);
  // 1 - We check the resource pointed by m_ptr has been moved from algo to
  // moved_algo
  REQUIRE(p1 == p2);
  // 2 - We check that the two string representations are identical
  REQUIRE(pla_string == moved_pla_string);
}

TEST_CASE("planet_move_assignment_test") {
  // We instantiate a planet
  complete_udpla udpla({1., 2., -1., 4.});
  planet pla{udpla};

  // We store a streaming representation of the object
  auto pla_string = boost::lexical_cast<std::string>(pla);
  // We get the memory address where the user algo is stored
  auto *p1 = pla.extract<complete_udpla>();
  // We call the move assignment
  planet moved_pla{};
  moved_pla = std::move(pla);
  // We get the memory address where the user algo is stored
  auto *p2 = moved_pla.extract<complete_udpla>();
  // And the string representation of the moved algo
  auto moved_pla_string = boost::lexical_cast<std::string>(moved_pla);
  // 1 - We check the resource pointed by m_ptr has been moved from algo to
  // moved_algo
  REQUIRE(p1 == p2);
  // 2 - We check that the two string representations are identical
  REQUIRE(pla_string == moved_pla_string);
}

TEST_CASE("copy_assignment_test") {
  // We instantiate a planet
  complete_udpla udpla({1., 2., -1., 4.});
  planet pla{udpla};

  // We call the copy assignment
  planet pla_copy{};
  pla_copy = pla;
  // We extract the user planet
  auto *p1 = pla.extract<complete_udpla>();
  auto *p2 = pla_copy.extract<complete_udpla>();

  // 1 - We check the resources pointed to by m_ptr have different address
  REQUIRE(p1 != 0);
  REQUIRE(p2 != 0);
  REQUIRE(p1 != p2);
  // 2 - We check that the other members are copied
  REQUIRE(pla.get_name() == pla_copy.get_name());
  REQUIRE(pla.get_mu_central_body() == pla_copy.get_mu_central_body());
  REQUIRE(pla.get_mu_self() == pla_copy.get_mu_self());
  REQUIRE(pla.get_safe_radius() == pla_copy.get_safe_radius());
}

TEST_CASE("planet_extract_is_test") {
  // We instantiate a planet
  planet pla{complete_udpla({1., 2., -1., 4.})};

  auto *p0 = pla.extract<complete_udpla>();

  // We check thet we can access to public data members
  REQUIRE(p0->m_mu_central_body == 1.);
  REQUIRE(p0->m_mu_self == 2.);
  REQUIRE(p0->m_radius == -1.);
  REQUIRE(p0->m_safe_radius == 4.);

  // We check that a non successful cast returns a null pointer
  REQUIRE(!pla.extract<simple_udpla>());

  // We check the is method
  REQUIRE(pla.is<complete_udpla>());
  REQUIRE(!pla.is<simple_udpla>());
}

TEST_CASE("is_valid") {
  planet p0;
  REQUIRE(p0.is_valid());
  planet p1(std::move(p0));
  REQUIRE(!p0.is_valid());
  p0 = planet{simple_udpla{}};
  REQUIRE(p0.is_valid());
  p1 = std::move(p0);
  REQUIRE(!p0.is_valid());
}

TEST_CASE("generic_assignment") {
  planet p0;
  REQUIRE(p0.is<null_udpla>());
  REQUIRE(&(p0 = simple_udpla{}) == &p0);
  REQUIRE(p0.is_valid());
  REQUIRE(p0.is<simple_udpla>());
  REQUIRE((!std::is_assignable<planet, void>::value));
  REQUIRE((!std::is_assignable<planet, int &>::value));
  REQUIRE((!std::is_assignable<planet, const int &>::value));
  REQUIRE((!std::is_assignable<planet, int &&>::value));
}

TEST_CASE("type_index") {
  planet p0 = planet{null_udpla{}};
  // REQUIRE(p0.get_type_index() == std::type_index(typeid(null_udpla)));
  // p0 = planet{simple_udpla{}};
  // REQUIRE(p0.get_type_index() == std::type_index(typeid(simple_udpla)));
}

TEST_CASE("get_ptr") {
  planet p0;
  REQUIRE(p0.get_ptr() == p0.extract<null_udpla>());
  REQUIRE(static_cast<const planet &>(p0).get_ptr() ==
          static_cast<const planet &>(p0).extract<null_udpla>());
  p0 = planet{simple_udpla{}};
  REQUIRE(p0.get_ptr() == p0.extract<simple_udpla>());
  REQUIRE(static_cast<const planet &>(p0).get_ptr() ==
          static_cast<const planet &>(p0).extract<simple_udpla>());
}

TEST_CASE("stream_operator") {
  REQUIRE_NOTHROW((std::cout << planet{} << '\n'));
}

TEST_CASE("planet_astro_methods_test") {
  // We instantiate a planet
  planet pla{complete_udpla({1.1, 2.3, 4.02, 4.5})};
  // Test eph
  auto pos_vel = pla.eph(epoch(0.));
  REQUIRE(pos_vel[0] == std::array<double, 3>{1., 0., 0.});
  REQUIRE(pos_vel[1] == std::array<double, 3>{0., 1., 0.});
  REQUIRE(pla.get_name() == "A complete, albeit simple Planet");
  REQUIRE(pla.get_mu_central_body() == 1.1);
  REQUIRE(pla.get_mu_self() == 2.3);
  REQUIRE(pla.get_radius() == 4.02);
  REQUIRE(pla.get_safe_radius() == 4.5);
  REQUIRE(pla.period(kep3::epoch(0.)) == 0.);
}

TEST_CASE("serialization_test") {
  // Instantiate a planet
  planet pla{complete_udpla({1.1, 2.3, 4.02, 4.5})};

  // Store the string representation.
  std::stringstream ss;
  auto before = boost::lexical_cast<std::string>(pla);
  // Now serialize, deserialize and compare the result.
  {
    boost::archive::binary_oarchive oarchive(ss);
    oarchive << pla;
  }
  // Create a new planet object
  auto pla2 = planet{simple_udpla{}};
  boost::lexical_cast<std::string>(pla2); // triggers the streaming operator
  {
    boost::archive::binary_iarchive iarchive(ss);
    iarchive >> pla2;
  }
  auto after = boost::lexical_cast<std::string>(pla2);
  REQUIRE(before == after);
  // Check explicitly that the properties of base_p where restored as well.
  REQUIRE(pla.extract<complete_udpla>()->m_mu_central_body ==
          pla.extract<complete_udpla>()->m_mu_central_body);
  REQUIRE(pla.extract<complete_udpla>()->m_mu_self ==
          pla.extract<complete_udpla>()->m_mu_self);
  REQUIRE(pla.extract<complete_udpla>()->m_radius ==
          pla.extract<complete_udpla>()->m_radius);
  REQUIRE(pla.extract<complete_udpla>()->m_safe_radius ==
          pla.extract<complete_udpla>()->m_safe_radius);
}