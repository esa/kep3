// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <functional>
#include <stdexcept>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/core_astro/kepler_equations.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using Catch::Matchers::WithinAbs;
using kep3::propagate_lagrangian;
using kep3::propagate_lagrangian_u;

void test_propagate_lagrangian_tof_inversion(
    const std::function<void(std::array<std::array<double, 3>, 2> &, double, double)> &propagate,
    unsigned int N = 10000)
{
    { // Elliptical circular motion xy
        std::array<std::array<double, 3>, 2> pos_vel = {{{1., 0, 0.}, {0., 1., 0.}}};
        propagate(pos_vel, kep3::pi / 2., 1.);
        auto &[pos, vel] = pos_vel;

        REQUIRE_THAT(pos[0], WithinAbs(0., 1e-14));
        REQUIRE_THAT(pos[1], WithinAbs(1., 1e-14));
        REQUIRE_THAT(pos[2], WithinAbs(0., 1e-14));
        REQUIRE_THAT(vel[0], WithinAbs(-1., 1e-14));
        REQUIRE_THAT(vel[1], WithinAbs(0., 1e-14));
        REQUIRE_THAT(vel[2], WithinAbs(0., 1e-14));
    }
    { // Elliptical circular motion xy
        std::array<std::array<double, 3>, 2> pos_vel = {{{1., 0, 0.}, {0., 1., 0.}}};
        propagate(pos_vel, -kep3::pi / 2., 1.);
        auto &[pos, vel] = pos_vel;

        REQUIRE_THAT(pos[0], WithinAbs(0., 1e-14));
        REQUIRE_THAT(pos[1], WithinAbs(-1., 1e-14));
        REQUIRE_THAT(pos[2], WithinAbs(0., 1e-14));
        REQUIRE_THAT(vel[0], WithinAbs(1., 1e-14));
        REQUIRE_THAT(vel[1], WithinAbs(0., 1e-14));
        REQUIRE_THAT(vel[2], WithinAbs(0., 1e-14));
    }
    { // Elliptical circular motion xz
        std::array<std::array<double, 3>, 2> pos_vel = {{{1., 0, 0.}, {0., 0., 1.}}};
        propagate(pos_vel, kep3::pi / 2., 1.);
        auto &[pos, vel] = pos_vel;

        REQUIRE_THAT(pos[0], WithinAbs(0., 1e-14));
        REQUIRE_THAT(pos[1], WithinAbs(0., 1e-14));
        REQUIRE_THAT(pos[2], WithinAbs(1., 1e-14));
        REQUIRE_THAT(vel[0], WithinAbs(-1., 1e-14));
        REQUIRE_THAT(vel[1], WithinAbs(0., 1e-14));
        REQUIRE_THAT(vel[2], WithinAbs(0., 1e-14));
    }
    { // Elliptical circular motion yz
        std::array<std::array<double, 3>, 2> pos_vel = {{{0., 1, 0.}, {0., 0., 1.}}};
        propagate(pos_vel, kep3::pi / 2., 1.);
        auto &[pos, vel] = pos_vel;

        REQUIRE_THAT(pos[0], WithinAbs(0., 1e-14));
        REQUIRE_THAT(pos[1], WithinAbs(0., 1e-14));
        REQUIRE_THAT(pos[2], WithinAbs(1., 1e-14));
        REQUIRE_THAT(vel[0], WithinAbs(0., 1e-14));
        REQUIRE_THAT(vel[1], WithinAbs(-1., 1e-14));
        REQUIRE_THAT(vel[2], WithinAbs(0., 1e-14));
    }
    // Engines
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(1220202343u);

    { // Targeting Ellipses
        std::uniform_real_distribution<double> sma_d(1.1, 10.);
        std::uniform_real_distribution<double> ecc_d(0, 0.9);
        std::uniform_real_distribution<double> incl_d(0., kep3::pi);
        std::uniform_real_distribution<double> Omega_d(0, 2 * kep3::pi);
        std::uniform_real_distribution<double> omega_d(0., kep3::pi);
        std::uniform_real_distribution<double> f_d(0, 2 * kep3::pi);
        std::uniform_real_distribution<double> time_d(-2. * kep3::pi, 2. * kep3::pi);

        // Testing on N random calls on ellipses
        for (auto i = 0u; i < N; ++i) {
            double sma = sma_d(rng_engine);
            double ecc = ecc_d(rng_engine);
            double incl = incl_d(rng_engine);
            double Omega = Omega_d(rng_engine);
            double omega = omega_d(rng_engine);
            double f = f_d(rng_engine);

            std::array<double, 6> par = {sma, ecc, incl, Omega, omega, f};
            auto pos_vel = kep3::par2ic(par, 1.);
            double tof = time_d(rng_engine);
            auto pos_vel_after = pos_vel;
            propagate(pos_vel_after, tof, 1.);
            propagate(pos_vel_after, -tof, 1.);
            REQUIRE(kep3_tests::floating_point_error_vector(pos_vel[0], pos_vel_after[0]) < 1e-13);
        }
    }

    { // Targeting Hyperbolas
        std::uniform_real_distribution<double> sma_d(1.1, 100.);
        std::uniform_real_distribution<double> ecc_d(2., 20.);
        std::uniform_real_distribution<double> incl_d(0., kep3::pi);
        std::uniform_real_distribution<double> Omega_d(0, 2 * kep3::pi);
        std::uniform_real_distribution<double> omega_d(0., kep3::pi);
        std::uniform_real_distribution<double> f_d(0, 2 * kep3::pi);
        std::uniform_real_distribution<double> time_d(0.1, 20.);
        // Testing on N random calls on hyperbolas
        for (auto i = 0u; i < N; ++i) {
            double sma = sma_d(rng_engine);
            double ecc = ecc_d(rng_engine);
            double incl = incl_d(rng_engine);
            double Omega = Omega_d(rng_engine);
            double omega = omega_d(rng_engine);
            double f = f_d(rng_engine);
            if (std::cos(f) > -1 / ecc) {
                std::array<double, 6> par = {-sma, ecc, incl, Omega, omega, f};
                auto pos_vel = kep3::par2ic(par, 1.);
                double tof = time_d(rng_engine);
                auto pos_vel_after = pos_vel;
                propagate(pos_vel_after, tof, 1.);
                propagate(pos_vel_after, -tof, 1.);
                // REQUIRE(kep3_tests::floating_point_error_vector(pos_vel[0], pos_vel_after[0]) < 1e-13);
            }
        }
    }
}

TEST_CASE("propagate_lagrangian")
{
    // We test both Normal and Universal variables versions.
    test_propagate_lagrangian_tof_inversion(&propagate_lagrangian, 10000u);
    test_propagate_lagrangian_tof_inversion(&propagate_lagrangian_u, 10000u);
}

TEST_CASE("correctness")
{
    // We test propagate_lagrangian and propagate_lagrangian_u give the same result.
    {
        // Some ellipse
        std::array<std::array<double, 3>, 2> pos_vel = {{{1.223, 0.3123, -0.432}, {0.06345, 0.43234, -0.874634}}};
        auto pos_vel_u = pos_vel;
        propagate_lagrangian(pos_vel, 3.56, 1.24);
        propagate_lagrangian_u(pos_vel_u, 3.56, 1.24);

        auto &[r, v] = pos_vel;
        auto &[ru, vu] = pos_vel_u;
        REQUIRE(kep3_tests::floating_point_error_vector(r, ru) < 1e-13);
        REQUIRE(kep3_tests::floating_point_error_vector(v, vu) < 1e-13);
    }
    {
        // Some hyperbolae
        std::array<std::array<double, 3>, 2> pos_vel = {{{1.223, 0.3123, -0.432}, {-3.06345, 4.43234, -0.874634}}};
        auto pos_vel_u = pos_vel;
        propagate_lagrangian(pos_vel, 3.56, 1.24);
        propagate_lagrangian_u(pos_vel_u, 3.56, 1.24);

        auto &[r, v] = pos_vel;
        auto &[ru, vu] = pos_vel_u;
        REQUIRE(kep3_tests::floating_point_error_vector(r, ru) < 1e-13);
        REQUIRE(kep3_tests::floating_point_error_vector(v, vu) < 1e-13);
    }
}
