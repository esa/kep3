// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <stdexcept>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/lambert_problem.hpp>
#include <kep3/leg/sims_flanagan.hpp>
#include <kep3/planet.hpp>
#include <kep3/udpla/vsop2013.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

TEST_CASE("constructor")
{
    {
        // The default constructor constructs a valid leg with no mismatches.
        kep3::leg::sims_flanagan sf{};
        auto mc = sf.compute_mismatch_constraints();
        REQUIRE(*std::max_element(mc.begin(), mc.end()) < 1e-15);
        auto tc = sf.compute_throttle_constraints();
        REQUIRE(*std::max_element(tc.begin(), tc.end()) < 0.);
    }
    {
        // The constructor fails when data are malformed
        std::array<double, 7> xs{1, 0, 0, 0, 1, 0, 1};
        std::array<double, 7> xf{0, 1, 0, -1, 0, 0, 1};
        REQUIRE_NOTHROW(kep3::leg::sims_flanagan(xs, {0, 0, 0, 0, 0, 0}, xf, kep3::pi / 2, 1., 1., 1., 0.5));
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan(xs, {0., 0., 0., 0., 0.}, xf, kep3::pi / 2, 1., 1., 1., 0.5),
                          std::logic_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan(xs, {0, 0, 0, 0, 0, 0}, xf, -0.42, 1., 1., 1., 0.5),
                          std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan(xs, {0, 0, 0, 0, 0, 0}, xf, kep3::pi / 2, -0.3, 1., 1., 0.5),
                          std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan(xs, {0, 0, 0, 0, 0, 0}, xf, kep3::pi / 2, 1., -2., 1., 0.5),
                          std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan(xs, {0, 0, 0, 0, 0, 0}, xf, kep3::pi / 2, 1., 1., -0.32, 0.5),
                          std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan(xs, {0, 0, 0, 0, 0, 0}, xf, kep3::pi / 2, 1., 1., 1., 32),
                          std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan(xs, {0, 0, 0, 0, 0, 0}, xf, kep3::pi / 2, 1., 1., 1., -0.1),
                          std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan(xs, {}, xf, kep3::pi / 2, 1., 1., 1., 0.5), std::logic_error);
    }
}

TEST_CASE("getters_and_setters")
{
    kep3::leg::sims_flanagan sf{};
    std::array<double, 7> xf = {1., 1., 1., 1., 1., 1., 1.};
    sf.set_xf(xf);
    REQUIRE(sf.get_xf() == xf);
    sf.set_xs(xf);
    REQUIRE(sf.get_xs() == xf);
    std::vector<double> throttles{1., 2., 3., 1., 2., 3.};
    sf.set_throttles(throttles);
    REQUIRE(sf.get_throttles() == throttles);
    sf.set_cut(0.333);
    REQUIRE(sf.get_cut() == 0.333);
    sf.set_max_thrust(0.333);
    REQUIRE(sf.get_max_thrust() == 0.333);
    sf.set_isp(0.333);
    REQUIRE(sf.get_isp() == 0.333);
    sf.set_mu(0.333);
    REQUIRE(sf.get_mu() == 0.333);
    sf.set_tof(0.333);
    REQUIRE(sf.get_tof() == 0.333);
}

TEST_CASE("compute_throttle_constraints_test")
{
    kep3::leg::sims_flanagan sf({1, 0, 0, 0, 1, 0, 1}, {0, 1, 0, 1, 1, 1, 0, 1, 1}, {0, 1, 0, -1, 0, 0, 1}, 1, 1, 1, 1,
                                1);
    auto tc = sf.compute_throttle_constraints();
    REQUIRE(tc[0] == 0.);
    REQUIRE(tc[1] == 2.);
    REQUIRE(tc[2] == 1.);
}

std::array<double, 7> to_x(std::array<std::array<double, 3>, 2> rv, double mass)
{
    return {rv[0][0], rv[0][1], rv[0][2], rv[1][0], rv[1][1], rv[1][2], mass};
}

std::array<double, 7> normalize_x(std::array<double, 7> x)
{
    x[0] /= kep3::AU;
    x[1] /= kep3::AU;
    x[2] /= kep3::AU;
    x[3] /= kep3::EARTH_VELOCITY;
    x[4] /= kep3::EARTH_VELOCITY;
    x[5] /= kep3::EARTH_VELOCITY;
    x[0] /= 1000.;
    return x;
}

TEST_CASE("compute_mismatch_constraints_test")
{
    {
        // We test that an engineered ballistic arc always returns no mismatch for all cuts.
        // We use (for no reason) the ephs of the Earth and Jupiter
        kep3::udpla::vsop2013 udpla_earth("earth_moon", 1e-2);
        kep3::udpla::vsop2013 udpla_jupiter("jupiter", 1e-2);
        kep3::planet earth{udpla_earth};
        kep3::planet jupiter{udpla_jupiter};
        // And some epochs / tofs.
        double dt_days = 1000.;
        double dt = dt_days * kep3::DAY2SEC;
        double t0 = 1233.3;
        double mass = 1000;
        auto rv0 = earth.eph(t0);
        auto rv1 = jupiter.eph(t0 + dt_days);
        // We create a ballistic arc matching the two.
        kep3::lambert_problem lp{rv0[0], rv1[0], dt, kep3::MU_SUN};
        rv0[1][0] = lp.get_v0()[0][0];
        rv0[1][1] = lp.get_v0()[0][1];
        rv0[1][2] = lp.get_v0()[0][2];
        rv1[1][0] = lp.get_v1()[0][0];
        rv1[1][1] = lp.get_v1()[0][1];
        rv1[1][2] = lp.get_v1()[0][2];
        // We test for 1 to 33 segments and cuts in [0,0.1,0.2, ..., 1]
        std::vector<double> cut_values{0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
        for (unsigned long N = 1u; N < 34; ++N) {
            for (auto cut : cut_values) {
                std::vector<double> throttles(N * 3, 0.);
                kep3::leg::sims_flanagan sf(to_x(rv0, mass), throttles, to_x(rv1, mass), dt, 1., 1., kep3::MU_SUN, cut);
                auto mc = sf.compute_mismatch_constraints();
                mc = normalize_x(mc);
                REQUIRE(*std::max_element(mc.begin(), mc.end()) < 1e-12);
            }
        }
    }
}
