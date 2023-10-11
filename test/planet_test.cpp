// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <stdexcept>
#include <string>

#include <fmt/core.h>

#include "catch.hpp"
#include <boost/core/demangle.hpp>
#include <boost/lexical_cast.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/epoch.hpp>
#include <kep3/exceptions.hpp>
#include <kep3/planet.hpp>

#include "test_helpers.hpp"

using kep3::epoch;
using kep3::planet;
using kep3::detail::null_udpla;

struct simple_udpla {
    simple_udpla() = default;
    static std::array<std::array<double, 3>, 2> eph(const epoch &)
    {
        std::array<double, 3> pos = {1., 0., 0.};
        std::array<double, 3> vel = {0., 1., 0.};
        return {pos, vel};
    };
    static std::string get_name()
    {
        return "A simple planet";
    }
    static std::string get_extra_info()
    {
        return "The simplest planet ever!";
    }

private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &, unsigned)
    {
    }
};

TANUKI_S11N_WRAP_EXPORT(simple_udpla, kep3::detail::planet_iface)

struct simple_udpla_mu {
    simple_udpla_mu() = default;
    static std::array<std::array<double, 3>, 2> eph(const epoch &)
    {
        std::array<double, 3> pos = {1., 0., 0.};
        std::array<double, 3> vel = {0., 1., 0.};
        return {pos, vel};
    };
    static std::string get_name()
    {
        return "A simple planet with mu";
    }
    static std::string get_extra_info()
    {
        return "The simplest planet ever but with mu!";
    }
    static double get_mu_central_body()
    {
        return 1.;
    }

private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &, unsigned)
    {
    }
};

TANUKI_S11N_WRAP_EXPORT(simple_udpla_mu, kep3::detail::planet_iface)

struct simple_udpla_mu_h {
    simple_udpla_mu_h() = default;
    static std::array<std::array<double, 3>, 2> eph(const epoch &)
    {
        std::array<double, 3> pos = {1., 0., 0.};
        std::array<double, 3> vel = {0., 10., 0.};
        return {pos, vel};
    };
    static std::string get_name()
    {
        return "A simple planet with mu";
    }
    static std::string get_extra_info()
    {
        return "The simplest planet ever but with mu!";
    }
    static double get_mu_central_body()
    {
        return 1.;
    }

private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &, unsigned)
    {
    }
};

TANUKI_S11N_WRAP_EXPORT(simple_udpla_mu_h, kep3::detail::planet_iface)

struct complete_udpla {
    explicit complete_udpla(std::array<double, 4> physical_properties = {-1., -1., -1., -1.})
        : m_name("A complete, albeit simple Planet"), m_mu_central_body(physical_properties[0]),
          m_mu_self(physical_properties[1]), m_radius(physical_properties[2]), m_safe_radius(physical_properties[3]){};
    static std::array<std::array<double, 3>, 2> eph(const epoch &)
    {
        std::array<double, 3> pos = {1., 0., 0.};
        std::array<double, 3> vel = {0., 1., 0.};
        return {pos, vel};
    };

    [[nodiscard]] std::string get_name() const
    {
        return m_name;
    }
    [[nodiscard]] double get_mu_central_body() const
    {
        return m_mu_central_body;
    }
    [[nodiscard]] double get_mu_self() const
    {
        return m_mu_self;
    }
    [[nodiscard]] double get_radius() const
    {
        return m_radius;
    }
    [[nodiscard]] double get_safe_radius() const
    {
        return m_safe_radius;
    }
    [[nodiscard]] double period(const epoch &) const
    {
        return m_radius - m_radius;
    }

    [[nodiscard]] std::array<double, 6> elements(kep3::epoch, kep3::elements_type) const
    {
        return {m_radius, m_radius, m_radius, m_radius, m_radius, m_radius};
    }

    [[nodiscard]] std::string get_extra_info() const
    {
        return fmt::format("Gravitational parameter of main attracting body: {}\nGravitational "
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
    template <typename Archive>
    void serialize(Archive &ar, unsigned)
    {
        ar &m_name;
        ar &m_mu_central_body;
        ar &m_mu_self;
        ar &m_radius;
        ar &m_safe_radius;
    }
};

TANUKI_S11N_WRAP_EXPORT(complete_udpla, kep3::detail::planet_iface)

TEST_CASE("construction")
{
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
        REQUIRE(pla.get_name() == boost::core::demangle(typeid(null_udpla).name()));
        REQUIRE(pla.get_extra_info() == std::string(""));
        REQUIRE(pla.get_mu_central_body() == -1);
        REQUIRE(pla.get_mu_self() == -1);
        REQUIRE(pla.get_radius() == -1);
        REQUIRE(pla.get_safe_radius() == -1);
        REQUIRE_THROWS_AS((pla.period()), kep3::not_implemented_error);
        REQUIRE_THROWS_AS((pla.elements()), kep3::not_implemented_error);
        REQUIRE(value_isa<null_udpla>(pla));
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
        REQUIRE(pla.get_mu_central_body() == -1);
        REQUIRE(pla.get_mu_self() == -1);
        REQUIRE(pla.get_radius() == -1);
        REQUIRE(pla.get_safe_radius() == -1);
        REQUIRE_THROWS_AS((pla.period()), kep3::not_implemented_error);
        REQUIRE_THROWS_AS((pla.elements()), kep3::not_implemented_error);

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
        REQUIRE(pla.get_radius() == -1);
        REQUIRE(pla.get_safe_radius() == 4.);
        REQUIRE(std::isfinite(pla.period()));
        REQUIRE(pla.elements() == std::array<double, 6>{-1,-1,-1,-1,-1,-1});

    }
    {
        // Constructor from a simple udpla with mu and hyperbolic orbit
        simple_udpla_mu_h udpla{};
        REQUIRE_NOTHROW(planet(udpla));
        planet pla(udpla);
        REQUIRE(!std::isfinite(pla.period()));
        REQUIRE(pla.elements()[0] < 0);
        REQUIRE_THROWS_AS(pla.elements(kep3::epoch(), kep3::elements_type::KEP_M), std::logic_error);
        REQUIRE(pla.elements(kep3::epoch(), kep3::elements_type::MEQ)[0] > 0);
        REQUIRE(pla.elements(kep3::epoch(), kep3::elements_type::MEQ_R)[0] > 0);
        REQUIRE(pla.elements(kep3::epoch(), kep3::elements_type::KEP_F)[0] < 0);

    }
    {
        // Constructor from a simple udpla with mu and elliptical orbit
        simple_udpla_mu udpla{};
        REQUIRE_NOTHROW(planet(udpla));
        planet pla(udpla);
        REQUIRE(kep3_tests::floating_point_error(pla.period(), kep3::pi * 2.) < 1e-14);
        REQUIRE(pla.elements()[0] > 0);
        REQUIRE(pla.elements()[1] < 1);
        REQUIRE(pla.elements(kep3::epoch(), kep3::elements_type::KEP_M)[0] > 0);
        REQUIRE(pla.elements(kep3::epoch(), kep3::elements_type::MEQ)[0] > 0);
        REQUIRE(pla.elements(kep3::epoch(), kep3::elements_type::MEQ_R)[0] > 0);
        REQUIRE(pla.elements(kep3::epoch(), kep3::elements_type::KEP_F)[0] > 0);
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
        REQUIRE(pla.get_radius() == -1);
        REQUIRE(pla.get_safe_radius() == 4.);
    }
    // Check copy semantics.
    auto p0 = planet();
    planet p1{p0};
    REQUIRE(value_isa<null_udpla>(p0));
    REQUIRE(value_isa<null_udpla>(p1));
    planet p2{simple_udpla{}};
    p2 = p1;
    REQUIRE(value_isa<null_udpla>(p2));
    //  Move semantics.
    planet p3{std::move(p0)};
    REQUIRE(value_isa<null_udpla>(p3));
    planet p4{simple_udpla{}};
    p4 = std::move(p2);
    REQUIRE(value_isa<null_udpla>(p4));
    //  Check we can revive moved-from objects.
    p0 = p4;
    REQUIRE(value_isa<null_udpla>(p0));
    p2 = std::move(p4);
    REQUIRE(value_isa<null_udpla>(p2));
}

TEST_CASE("planet_extract_is_test")
{
    // We instantiate a planet
    planet pla{complete_udpla({1., 2., -1., 4.})};

    auto *p0 = value_ptr<complete_udpla>(pla);

    // We check thet we can access to public data members
    REQUIRE(p0->m_mu_central_body == 1.);
    REQUIRE(p0->m_mu_self == 2.);
    REQUIRE(p0->m_radius == -1.);
    REQUIRE(p0->m_safe_radius == 4.);

    // We check that a non successful cast returns a null pointer
    REQUIRE(!value_ptr<simple_udpla>(pla));

    // We check the is method
    REQUIRE(value_isa<complete_udpla>(pla));
    REQUIRE(!value_isa<simple_udpla>(pla));
}

TEST_CASE("generic_assignment")
{
    REQUIRE((!std::is_assignable<planet, void>::value));
    REQUIRE((!std::is_assignable<planet, int &>::value));
    REQUIRE((!std::is_assignable<planet, const int &>::value));
    REQUIRE((!std::is_assignable<planet, int &&>::value));
}

TEST_CASE("stream_operator")
{
    REQUIRE_NOTHROW((std::cout << planet{} << '\n'));
}

TEST_CASE("planet_astro_methods_test")
{
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

TEST_CASE("serialization_test")
{
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
    REQUIRE(value_ref<complete_udpla>(pla).m_mu_central_body == value_ref<complete_udpla>(pla2).m_mu_central_body);
    REQUIRE(value_ref<complete_udpla>(pla).m_mu_self == value_ref<complete_udpla>(pla2).m_mu_self);
    REQUIRE(value_ref<complete_udpla>(pla).m_radius == value_ref<complete_udpla>(pla2).m_radius);
    REQUIRE(value_ref<complete_udpla>(pla).m_safe_radius == value_ref<complete_udpla>(pla2).m_safe_radius);
}
