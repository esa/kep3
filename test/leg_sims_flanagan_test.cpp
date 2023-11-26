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
#include <utility>
#include <vector>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/ipopt.hpp>
#include <pagmo/algorithms/nlopt.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/utils/gradients_and_hessians.hpp>

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
        std::array<std::array<double, 3>, 2> rvs{{{1, 0, 0}, {0, 1, 0}}};
        std::array<std::array<double, 3>, 2> rvf{{{0, 1, 0}, {-1, 0, 0}}};
        double ms = 1.;
        double mf = 1.;
        REQUIRE_NOTHROW(kep3::leg::sims_flanagan(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, kep3::pi / 2, 1., 1., 1., 0.5));
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan(rvs, ms, {0., 0., 0., 0., 0.}, rvf, mf, kep3::pi / 2, 1., 1., 1., 0.5),
            std::logic_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, -0.42, 1., 1., 1., 0.5),
                          std::domain_error);
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, kep3::pi / 2, -0.3, 1., 1., 0.5),
            std::domain_error);
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, kep3::pi / 2, 1., -2., 1., 0.5),
            std::domain_error);
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, kep3::pi / 2, 1., 1., -0.32, 0.5),
            std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, kep3::pi / 2, 1., 1., 1., 32),
                          std::domain_error);
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, kep3::pi / 2, 1., 1., 1., -0.1),
            std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan(rvs, ms, {}, rvf, mf, kep3::pi / 2, 1., 1., 1., 0.5),
                          std::logic_error);
    }
}

TEST_CASE("getters_and_setters")
{
    kep3::leg::sims_flanagan sf{};
    std::array<std::array<double, 3>, 2> rvf{{{1, 1, 1}, {1, 1, 1}}};
    double mass = 123.;
    sf.set_rvf(rvf);
    REQUIRE(sf.get_rvf() == rvf);
    sf.set_ms(mass);
    REQUIRE(sf.get_ms() == mass);
    sf.set_rvs(rvf);
    REQUIRE(sf.get_rvs() == rvf);
    sf.set_mf(mass);
    REQUIRE(sf.get_mf() == mass);
    std::vector<double> throttles{1., 2., 3., 1., 2., 3.};
    std::vector<double> throttles2{1.1, 2.1, 3.1, 1.1, 2.1, 3.1};
    sf.set_throttles(throttles);
    REQUIRE(sf.get_throttles() == throttles);
    sf.set_throttles(throttles2.begin(), throttles2.end());
    REQUIRE(sf.get_throttles() == throttles2);
    REQUIRE_THROWS_AS(sf.set_throttles(throttles2.begin(), throttles2.end() - 1), std::logic_error);
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
    std::array<std::array<double, 3>, 2> rvs{{{1, 0, 0}, {0, 1, 0}}};
    std::array<std::array<double, 3>, 2> rvf{{{0, 1, 0}, {-1, 0, 0}}};
    kep3::leg::sims_flanagan sf(rvs, 1., {0, 1, 0, 1, 1, 1, 0, 1, 1}, rvf, 1, 1, 1, 1, 1, 1);
    auto tc = sf.compute_throttle_constraints();
    REQUIRE(tc[0] == 0.);
    REQUIRE(tc[1] == 2.);
    REQUIRE(tc[2] == 1.);
}

std::array<double, 7> normalize_con(std::array<double, 7> con)
{
    con[0] /= kep3::AU;
    con[1] /= kep3::AU;
    con[2] /= kep3::AU;
    con[3] /= kep3::EARTH_VELOCITY;
    con[4] /= kep3::EARTH_VELOCITY;
    con[5] /= kep3::EARTH_VELOCITY;
    con[6] /= 1.;
    return con;
}

struct sf_test_udp {
    sf_test_udp() = default;
    sf_test_udp(std::array<std::array<double, 3>, 2> rvs, double ms, std::array<std::array<double, 3>, 2> rvf,
                // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                double max_thrust = 0.05, double isp = 2000, unsigned nseg = 10u)
        : m_rvs(rvs), m_rvf(rvf), m_ms(ms), m_max_thrust(max_thrust), m_isp(isp), m_nseg(nseg)
    {
    }

    [[nodiscard]] std::vector<double> fitness(const std::vector<double> &x) const
    {
        // x = [tof (in days), mf (in kg), throttles]
        // We set the leg (avoiding the allocation for the throttles is possible but requires mutable data members.)
        kep3::leg::sims_flanagan leg(m_rvs, m_ms, std::vector<double>(m_nseg * 3, 0.), m_rvf, m_ms,
                                     x[0] * kep3::DAY2SEC, m_max_thrust, m_isp, kep3::MU_SUN);
        // We compute segments and dt
        std::size_t nseg = leg.get_throttles().size() / 3u;
        double dt = leg.get_tof() / static_cast<double>(nseg);
        // We set the throttles
        leg.set_throttles(x.begin() + 1, x.end());
        // We compute the mass schedule
        double mass = m_ms;

        for (decltype(leg.get_throttles().size()) i = 0u; i < nseg; ++i) {
            // We compute the the dv
            double dv0 = leg.get_max_thrust() / mass * dt * leg.get_throttles()[3 * i];
            double dv1 = leg.get_max_thrust() / mass * dt * leg.get_throttles()[3 * i + 1];
            double dv2 = leg.get_max_thrust() / mass * dt * leg.get_throttles()[3 * i + 2];
            // Update the mass
            double norm_dv = std::sqrt(dv0 * dv0 + dv1 * dv1 + dv2 * dv2);
            mass *= std::exp(-norm_dv / leg.get_isp() / kep3::G0);
        }
        leg.set_mf(mass);

        std::vector<double> retval(7 + m_nseg, 0.);

        // Fitness
        retval[0] = -mass;
        // Equality Constraints
        auto eq_con = leg.compute_mismatch_constraints();
        retval[1] = eq_con[0] / kep3::AU;
        retval[2] = eq_con[1] / kep3::AU;
        retval[3] = eq_con[2] / kep3::AU;
        retval[4] = eq_con[3] / kep3::EARTH_VELOCITY;
        retval[5] = eq_con[4] / kep3::EARTH_VELOCITY;
        retval[6] = eq_con[5] / kep3::EARTH_VELOCITY;
        //  Inequality Constraints
        auto ineq_con = leg.compute_throttle_constraints();
        std::copy(ineq_con.begin(), ineq_con.end(), retval.begin() + 7);
        return retval;
    }

    [[nodiscard]] std::vector<double> gradient(const std::vector<double> &x) const
    {
        return pagmo::estimate_gradient([this](const std::vector<double> &x) { return this->fitness(x); }, x);
    }

    [[nodiscard]] std::pair<std::vector<double>, std::vector<double>> get_bounds() const
    {
        std::vector<double> lb(m_nseg * 3 + 1, -1.);
        std::vector<double> ub(m_nseg * 3 + 1, +1.);
        lb[0] = 1.;    // days
        ub[0] = 2500.; // days
        return {lb, ub};
    }

    [[nodiscard]] static std::vector<double>::size_type get_nec()
    {
        return 6u;
    }

    [[nodiscard]] std::vector<double>::size_type get_nic() const
    {
        return m_nseg;
    }

    std::array<std::array<double, 3>, 2> m_rvs{};
    std::array<std::array<double, 3>, 2> m_rvf{};
    double m_ms{};
    double m_max_thrust{};
    double m_isp{};
    std::size_t m_nseg{};
};

TEST_CASE("compute_mismatch_constraints_test")
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
            kep3::leg::sims_flanagan sf(rv0, 1., throttles, rv1, 1., dt, 1., 1., kep3::MU_SUN, cut);
            auto mc = sf.compute_mismatch_constraints();
            mc = normalize_con(mc);
            REQUIRE(*std::max_element(mc.begin(), mc.end()) < 1e-12);
        }
    }

    {
        pagmo::problem prob{sf_test_udp{rv0, mass, rv1}};
        std::cout << prob << std::endl;
        pagmo::population pop{prob, 1u};
        pagmo::nlopt uda{"slsqp"};
        // algo.set_integer_option("print_level", 5);
        // algo.set_integer_option("max_iter", 3000);
        uda.set_xtol_abs(0.);
        uda.set_xtol_rel(0.);
        uda.set_ftol_abs(1e-8);
        pagmo::algorithm algo{uda};

        algo.set_verbosity(1u);
        pop = algo.evolve(pop);
        fmt::print("{}", pop.champion_f());
    }
}
