// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <chrono>
#include <fmt/base.h>
#include <stdexcept>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/io/xio.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/ipopt.hpp>
#include <pagmo/algorithms/nlopt.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/lambert_problem.hpp>
#include <kep3/leg/sf_checks.hpp>
#include <kep3/leg/sims_flanagan_alpha.hpp>
#include <kep3/leg/sims_flanagan_hf.hpp>
#include <kep3/leg/sims_flanagan_hf_alpha.hpp>
#include <kep3/planet.hpp>
#include <kep3/ta/zero_hold_kep.hpp>
#include <kep3/udpla/vsop2013.hpp>

#include "catch.hpp"
#include "leg_sims_flanagan_hf_alpha_helpers.hpp"
#include "test_helpers.hpp"
#include <pagmo/utils/gradients_and_hessians.hpp>
#include <xtensor/views/xview.hpp>

#include <heyoka/config.hpp>
#include <heyoka/expression.hpp>
#include <heyoka/math/pow.hpp>
#include <heyoka/math/relational.hpp>
#include <heyoka/math/select.hpp>
#include <heyoka/math/sqrt.hpp>
#include <heyoka/math/sum.hpp>
#include <heyoka/taylor.hpp>

#include "leg_sims_flanagan_hf_alpha_udp.hpp"

TEST_CASE("constructor")
{
    {
        // The default constructor constructs a valid leg with no mismatches.
        kep3::leg::sims_flanagan_hf_alpha sf{};
        auto mc = sf.compute_mismatch_constraints();
        REQUIRE(*std::max_element(mc.begin(), mc.end()) < 1e-13);
        auto tc = sf.compute_throttle_constraints();
        REQUIRE(*std::max_element(tc.begin(), tc.end()) < 0.);
    }
    {
        // The constructor fails when data are malformed
        std::array<std::array<double, 3>, 2> rvs{{{1, 0, 0}, {0, 1, 0}}};
        std::array<std::array<double, 3>, 2> rvf{{{0, 1, 0}, {-1, 0, 0}}};
        double ms = 1.;
        double mf = 1.;
        REQUIRE_NOTHROW(kep3::leg::sims_flanagan_hf_alpha(rvs, ms, {0., 0., 0., 0., 0., 0.}, {0., 0.}, rvf, mf,
                                                          kep3::pi / 2, 1., 1., 1., 0.5));
        const std::array<double, 7> rvms{1, 0, 0, 0, 1, 0, 1};
        const std::array<double, 7> rvmf{0, 1, 0, -1, 0, 0, 1};
        REQUIRE_NOTHROW(kep3::leg::sims_flanagan_hf_alpha(rvms, {0., 0., 0., 0., 0., 0.}, {0., 0.}, rvmf, kep3::pi / 2,
                                                          1., 1., 1., 0.5));
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan_hf_alpha(rvs, ms, {0., 0., 0., 0., 0.}, {0., 0.}, rvf, mf,
                                                            kep3::pi / 2, 1., 1., 1., 0.5),
                          std::logic_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan_hf_alpha(rvs, ms, {0., 0., 0., 0., 0., 0.}, {0.}, rvf, mf,
                                                            kep3::pi / 2, 1., 1., 1., 0.5),
                          std::logic_error);
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan_hf_alpha(rvs, ms, {0, 0, 0, 0, 0, 0}, {0., 0.}, rvf, mf, -0.42, 1., 1., 1., 0.5),
            std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan_hf_alpha(rvs, ms, {0, 0, 0, 0, 0, 0}, {0., 0.}, rvf, mf,
                                                            kep3::pi / 2, -0.3, 1., 1., 0.5),
                          std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan_hf_alpha(rvs, ms, {0, 0, 0, 0, 0, 0}, {0., 0.}, rvf, mf,
                                                            kep3::pi / 2, 1., -2., 1., 0.5),
                          std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan_hf_alpha(rvs, ms, {0, 0, 0, 0, 0, 0}, {0., 0.}, rvf, mf,
                                                            kep3::pi / 2, 1., 1., -0.32, 0.5),
                          std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan_hf_alpha(rvs, ms, {0, 0, 0, 0, 0, 0}, {0., 0.}, rvf, mf,
                                                            kep3::pi / 2, 1., 1., 1., 32),
                          std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan_hf_alpha(rvs, ms, {0, 0, 0, 0, 0, 0}, {0., 0.}, rvf, mf,
                                                            kep3::pi / 2, 1., 1., 1., -0.1),
                          std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan_hf_alpha(rvs, ms, {}, {}, rvf, mf, kep3::pi / 2, 1., 1., 1., 0.5),
                          std::logic_error);
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan_hf_alpha(rvs, ms, {0, 0, 0, 0, 0, 0}, {}, rvf, mf, kep3::pi / 2, 1., 1., 1., 0.5),
            std::logic_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan_hf_alpha(rvs, ms, {0, 0, 0, 0, 0, 0}, {0., 0.}, rvf, mf,
                                                            kep3::pi / 2, 1., 1., 1., 0.5, -1e-2),
                          std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan_hf_alpha(rvs, ms, {0, 0, 0, 0, 0, 0}, {0., 0.}, rvf, mf,
                                                            kep3::pi / 2, 1., 1., 1., 0.5, 1.2),
                          std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::_check_nseg(2, 1, 2), std::logic_error);
        // Create two Taylor adaptive integrators
        auto ta1 = kep3::ta::get_ta_zero_hold_kep(1e-13);
        auto ta2 = kep3::ta::get_ta_zero_hold_kep_var(1e-13);
        // Wrap them in an optional pair of references
        auto tas_opt = std::optional{
            std::pair<const heyoka::taylor_adaptive<double>&,
                    const heyoka::taylor_adaptive<double>&>(ta1, ta2)
        };
        REQUIRE_NOTHROW(kep3::leg::sims_flanagan_hf_alpha(rvs, ms, {0., 0., 0., 0., 0., 0.}, {0., 0.}, rvf, mf,
                                                          kep3::pi / 2, 1., 1., 1., 0.5, 1e-13,
                                                          tas_opt));
    }
}

TEST_CASE("getters_and_setters")
{
    {
        kep3::leg::sims_flanagan_hf_alpha sf{};
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

        std::vector<double> talphas{1., 2., 3.};
        std::vector<double> talphas2{1.1, 2.1, 3.1};
        sf.set_talphas(talphas);
        REQUIRE(sf.get_talphas() == talphas);
        sf.set_talphas(talphas2.begin(), talphas2.end());
        REQUIRE(sf.get_talphas() == talphas2);
        REQUIRE_THROWS_AS(sf.set_throttles(talphas2.begin(), talphas2.end() - 1), std::logic_error);

        sf.set_cut(0.333);
        REQUIRE(sf.get_cut() == 0.333);
        sf.set_max_thrust(0.333);
        REQUIRE(sf.get_max_thrust() == 0.333);
        sf.set_veff(0.333);
        REQUIRE(sf.get_veff() == 0.333);
        REQUIRE(sf.get_tas().get_pars()[1] == 0.333);
        REQUIRE(sf.get_tas_var().get_pars()[1] == 0.333);
        sf.set_mu(0.333);
        REQUIRE(sf.get_mu() == 0.333);
        REQUIRE(sf.get_tas().get_pars()[0] == 0.333);
        REQUIRE(sf.get_tas_var().get_pars()[0] == 0.333);
        sf.set_tof(0.333);
        REQUIRE(sf.get_tof() == 0.333);
        sf.set_tol(1e-4);
        REQUIRE(sf.get_tol() == 1e-4);
    }
    {
        kep3::leg::sims_flanagan_hf_alpha sf{};
        std::array<std::array<double, 3>, 2> rvf{{{1, 1, 1}, {1, 1, 1}}};
        std::vector<double> throttles{1., 2., 3., 1., 2., 3.};
        std::vector<double> talphas{1., 2.};

        sf.set(rvf, 12, throttles, talphas, rvf, 12, 4, 4, 4, 4, 0.333, 2e-5);
        REQUIRE(sf.get_rvs() == rvf);
        REQUIRE(sf.get_ms() == 12);
        REQUIRE(sf.get_rvf() == rvf);
        REQUIRE(sf.get_mf() == 12);
        REQUIRE(sf.get_throttles() == throttles);
        REQUIRE(sf.get_max_thrust() == 4);
        REQUIRE(sf.get_veff() == 4);
        REQUIRE(sf.get_mu() == 4);
        REQUIRE(sf.get_tof() == 4);
        REQUIRE(sf.get_cut() == 0.333);
        REQUIRE(sf.get_tol() == 2e-5);
    }
    {
        kep3::leg::sims_flanagan_hf_alpha sf{};
        std::array<double, 7> rvms{1, 1, 1, 1, 1, 1, 1};
        std::vector<double> throttles{1., 2., 3., 1., 2., 3.};
        std::vector<double> talphas{1., 2.};

        sf.set(rvms, throttles, talphas, rvms, 4, 4, 4, 4, 0.333, 2e-5);
        REQUIRE(sf.get_rvms() == rvms);
        REQUIRE(sf.get_rvmf() == rvms);
        REQUIRE(sf.get_throttles() == throttles);
        REQUIRE(sf.get_max_thrust() == 4);
        REQUIRE(sf.get_veff() == 4);
        REQUIRE(sf.get_mu() == 4);
        REQUIRE(sf.get_tof() == 4);
        REQUIRE(sf.get_cut() == 0.333);
        REQUIRE(sf.get_tol() == 2e-5);
        REQUIRE(typeid(sf.get_tas()) == typeid(kep3::ta::get_ta_zero_hold_kep(sf.get_tol())));
        REQUIRE(typeid(sf.get_tas_var()) == typeid(kep3::ta::get_ta_zero_hold_kep_var(sf.get_tol())));
    }
    {
        kep3::leg::sims_flanagan_hf_alpha sf{};
        std::array<double, 7> rvms{1, 1, 1, 1, 1, 1, 1};
        std::vector<double> throttles{1., 2., 3., 1., 2., 3.};
        std::vector<double> talphas{1., 2.};

        sf.set(rvms, throttles, talphas, rvms, 4);
        REQUIRE(sf.get_rvms() == rvms);
        REQUIRE(sf.get_rvmf() == rvms);
        REQUIRE(sf.get_tof() == 4);
        REQUIRE(typeid(sf.get_tas()) == typeid(kep3::ta::get_ta_zero_hold_kep(sf.get_tol())));
        REQUIRE(typeid(sf.get_tas_var()) == typeid(kep3::ta::get_ta_zero_hold_kep_var(sf.get_tol())));
    }
}

TEST_CASE("compute_throttle_constraints_test")
{
    std::array<std::array<double, 3>, 2> rvs{{{1, 0, 0}, {0, 1, 0}}};
    std::array<std::array<double, 3>, 2> rvf{{{0, 1, 0}, {-1, 0, 0}}};
    kep3::leg::sims_flanagan_hf_alpha sf(rvs, 1., {0, 1, 0, 1, 1, 1, 0, 1, 1}, {0.1, 0.1, 0.1}, rvf, 1, 1, 1, 1, 1, 1);
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
    con[6] /= 1000;
    return con;
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
        // double mass = 1000;
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
                std::vector<double> talphas(N, dt / static_cast<double>(N));
                kep3::leg::sims_flanagan_hf_alpha sf(rv0, 1., throttles, talphas, rv1, 1., dt, 1., 1., kep3::MU_SUN,
                                                     cut);
                auto mc = sf.compute_mismatch_constraints();
                auto mc_before = mc;
                mc = normalize_con(mc);
                auto tmp = *std::max_element(mc.begin(), mc.end());
                if (!std::isfinite(tmp)) {
                    fmt::print("\n\n");
                    fmt::print("Before Norm Mismatch constraint {}\n", mc_before);
                    fmt::print("Norm Mismatch constraint {}\n", mc);
                    fmt::print("Cut {}\n", cut);
                    fmt::print("rv0 {}\n", rv0);
                    fmt::print("Throttles {}\n", throttles);
                    fmt::print("talphas {}\n", talphas);
                    fmt::print("dt {}\n", dt);
                    fmt::print("rv1 {}\n", rv1);
                    fmt::print("\n\n");
                }
                REQUIRE(*std::max_element(mc.begin(), mc.end()) < 1e-8);
            }
        }
    }
    {
        // We test that some random thrusted arc computes correctly the mismatches
        std::array<std::array<double, 3>, 2> rv0{{{-25216645728.283768, 144924279081.32498, -38276.915766745136},
                                                  {-29833.034155296387, -5217.946770284042, 0.0013781466450451985}}};
        std::array<std::array<double, 3>, 2> rv1{{{207987766344.9237, -3139291734.542421, -5177822135.395065},
                                                  {1296.9096025329445, 26295.415668645317, 518.960634127031}}};
        std::vector<double> throttles{0.1, -0.2, 0.3, -0.4, 0.1, 0.2, -0.3, 0.2, 0.2, -0.1, -0.3, 0.1};
        // std::vector<double> throttles{0,0,0,0,0,0,0,0,0,0,0,0};
        double m0 = 4500.;
        double m1 = 4500.;
        double tof = 340. * kep3::DAY2SEC;
        std::vector<double> talphas{tof / 4, tof / 4, tof / 4, tof / 4};
        double max_thrust = 0.05;
        double veff = 2500. * kep3::G0;
        kep3::leg::sims_flanagan_hf_alpha sf(rv0, m0, throttles, talphas, rv1, m1, tof, max_thrust, veff,
                                             1.32712440018e+20, 0.5);
        auto mc = sf.compute_mismatch_constraints();
        // This was computed using directly an independent method (manually via Taylor integration)
        std::array<double, 7> ground_truth
            = {49962343234.42602, 63860682492.61005, 3188074669.721971,  4846.696643712443,
               2752.267007482824, 696.3982365414013, -23.610620941327397};
        REQUIRE(kep3_tests::L_infinity_norm_rel(mc, ground_truth) < 1e-13);
    }
}

TEST_CASE("UDP Fitness function timing")
{
    // We test that an engineered ballistic arc always returns no mismatch for all cuts.
    // We use (for no reason) the ephs of the Earth and Jupiter
    // Define the vectors
    std::array<std::array<double, 3>, 2> rv0 = {{{-125036811000.422, -83670919168.87277, 2610252.8064399767},
                                                 {16081.829029183446, -24868.923007449284, 0.7758272135425942}}};

    std::array<std::array<double, 3>, 2> rv1 = {{{-169327023332.1986, -161931354587.78766, 763967345.9733696},
                                                 {17656.297796509956, -15438.116653052988, -756.9165272457421}}};
    // And some epochs / tofs.
    double dt_days = 550.;
    double dt = dt_days * kep3::DAY2SEC;
    double mass = 1500;
    double max_thrust = 0.6;
    double veff = 3000.0;

    // Define optimization parameters
    size_t nseg = 8u;
    std::vector<double> lb(nseg * 4 + 2, 0.0);
    std::fill(lb.begin() + static_cast<long>(nseg) * 3, lb.begin() + static_cast<long>(nseg) * 4, 0.7);
    lb[nseg * 4] = 500.0;          // days
    lb[nseg * 4 + 1] = mass / 2.0; // kg

    std::vector<double> throttles(nseg * 3, 0.0);
    std::vector<double> talphas(nseg, dt / static_cast<double>(nseg));

    pagmo::problem prob{sf_hf_alpha_test_udp{rv0, mass, rv1, max_thrust, veff, 8u}};

    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();
    unsigned long Ntotal = 10000;
    for (unsigned long N = 1; N <= Ntotal; ++N) { // Start from 1 to avoid division by zero
        auto fitness_result = prob.fitness(lb);
    }
    // End timing
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::milliseconds>(end_time - start_time);

    // Print execution time
    fmt::print("Time taken for {} fitness function calls: {} ms\n", Ntotal, duration.count());
}

TEST_CASE("compute_mismatch_constraints_test_SLSQP")
{
    // We test that an engineered ballistic arc always returns no mismatch for all cuts.
    // We use (for no reason) the ephs of the Earth and Jupiter
    // Define the vectors
    std::array<std::array<double, 3>, 2> rv0 = {{{-125036811000.422, -83670919168.87277, 2610252.8064399767},
                                                 {16081.829029183446, -24868.923007449284, 0.7758272135425942}}};

    std::array<std::array<double, 3>, 2> rv1 = {{{-169327023332.1986, -161931354587.78766, 763967345.9733696},
                                                 {17656.297796509956, -15438.116653052988, -756.9165272457421}}};
    // And some epochs / tofs.
    double dt_days = 550.;
    double dt = dt_days * kep3::DAY2SEC;
    double mass = 1500;
    double max_thrust = 0.6;
    double veff = 3000.0 * kep3::G0;
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
            std::vector<double> talphas(N, dt / static_cast<double>(N));
            kep3::leg::sims_flanagan_hf_alpha sf(rv0, 1., throttles, talphas, rv1, 1., dt, 1., 1., kep3::MU_SUN, cut);
            auto mc = sf.compute_mismatch_constraints();
            mc = normalize_con(mc);
            REQUIRE(*std::max_element(mc.begin(), mc.end()) < 1e-8);
        }
    }

    {
        // Here we reuse the ballitic arc as a ground truth for an optimization.
        // We check that, when feasible, the optimal mass solution is indeed ballistic.
        // pagmo::problem prob{sf_hf_test_udp{rv0, mass, rv1, max_thrust, 2000, 10u}};
        pagmo::problem prob{sf_hf_alpha_test_udp{rv0, mass, rv1, max_thrust, veff, 8u}};
        prob.set_c_tol(1e-6);
        bool found = false;
        // pagmo::ipopt uda{};
        pagmo::nlopt uda{"slsqp"};
        uda.set_xtol_abs(1e-8);
        uda.set_xtol_rel(1e-8);
        uda.set_ftol_abs(0);
        uda.set_maxeval(1000);
        pagmo::algorithm algo{uda};
        std::vector<double> x0 = {4.4974400664853576e-06,
                                  -2.228606694039324e-05,
                                  -4.127503984025221e-06,
                                  -4.4974400664853576e-05,
                                  -3.962503037691836e-05,
                                  -5.53484490582199e-05,
                                  -4.1704137150877467e-05,
                                  -3.061691084379221e-05,
                                  0.00026180272039488274,
                                  1.4072230620899963e-05,
                                  2.1704137150877467e-05,
                                  -1.704323879748168e-05,
                                  8.099867069119767e-05,
                                  0.00022733814598852931,
                                  -0.0001560602131207774,
                                  -4.863705957848785e-05,
                                  4.9172235359169806e-05,
                                  3.125332265206107e-06,
                                  1.9893423744072716e-06,
                                  -8.678891965163936e-05,
                                  3.704323879748168e-05,
                                  -3.9172235359169806e-05,
                                  -2.5492703724246742e-05,
                                  -8.794501922065229e-05,
                                  0.72,
                                  0.8955821741388365,
                                  0.8871475866789136,
                                  0.8958681122690402,
                                  0.7955821741388365,
                                  0.8862111905777488,
                                  0.7327544903143346,
                                  0.8899466274787871,
                                  500.0084959099231,
                                  1230};
        pagmo::population pop{prob};
        pop.push_back(x0);
        algo.set_verbosity(10u);
        pop = algo.evolve(pop);
        auto champ = pop.champion_f();
        found = prob.feasibility_f(champ);

        fmt::print("Champ {}\n", champ);
        fmt::print("Best x{}\n", pop.champion_x());

        REQUIRE_FALSE(!found); // If this does not pass, then the optimization above never found a ballistic arc ...
                               // theres a problem somewhere.
    }
}

// Compare low-fidelity and high-fidelity methods with zero thrust (ought to be the same)
TEST_CASE("compare_low_and_high_fidelity_with_alpha")
{

    // Initialise unique test quantities
    double cut = 0.6;
    auto sf_helper_object = sf_hf_test_alpha_object(cut);

    kep3::leg::sims_flanagan_hf_alpha sf(sf_helper_object.m_rvs, sf_helper_object.m_ms, sf_helper_object.m_throttles,
                                         sf_helper_object.m_talphas, sf_helper_object.m_rvf, sf_helper_object.m_mf,
                                         sf_helper_object.m_tof, sf_helper_object.m_max_thrust, sf_helper_object.m_veff,
                                         sf_helper_object.m_mu, sf_helper_object.m_cut, 1e-16);
    kep3::leg::sims_flanagan_alpha sf_lf(sf_helper_object.m_rvs, sf_helper_object.m_ms, sf_helper_object.m_throttles,
                                         sf_helper_object.m_talphas, sf_helper_object.m_rvf, sf_helper_object.m_mf,
                                         sf_helper_object.m_tof, sf_helper_object.m_max_thrust, sf_helper_object.m_veff,
                                         sf_helper_object.m_mu, sf_helper_object.m_cut);

    auto retval = sf.compute_mismatch_constraints();
    auto retval_lf = sf_lf.compute_mismatch_constraints();

    std::array<double, 3> r1 = {retval[0], retval[1], retval[2]};
    std::array<double, 3> r2 = {retval_lf[0], retval_lf[1], retval_lf[2]};
    std::array<double, 3> v1 = {retval[3], retval[4], retval[5]};
    std::array<double, 3> v2 = {retval_lf[3], retval_lf[4], retval_lf[5]};

    REQUIRE(kep3_tests::floating_point_error_vector(r1, r2) < 1e-14);
    REQUIRE(kep3_tests::floating_point_error_vector(v1, v2) < 1e-14);
    REQUIRE(std::abs((retval[6] - retval_lf[6]) / retval[6]) < 1e-14);
}

// Compare low-fidelity and high-fidelity methods with zero thrust (ought to be the same)
TEST_CASE("compare_withandwithout_alpha")
{

    // Initialise unique test quantities
    double cut = 0.6;
    auto sf_helper_object = sf_hf_test_alpha_object(cut);

    kep3::leg::sims_flanagan_hf sf(sf_helper_object.m_rvs, sf_helper_object.m_ms, sf_helper_object.m_throttles,
                                   sf_helper_object.m_rvf, sf_helper_object.m_mf, sf_helper_object.m_tof,
                                   sf_helper_object.m_max_thrust, sf_helper_object.m_veff, sf_helper_object.m_mu,
                                   sf_helper_object.m_cut, 1e-16);
    kep3::leg::sims_flanagan_hf_alpha sf_alpha(
        sf_helper_object.m_rvs, sf_helper_object.m_ms, sf_helper_object.m_throttles, sf_helper_object.m_talphas,
        sf_helper_object.m_rvf, sf_helper_object.m_mf, sf_helper_object.m_tof, sf_helper_object.m_max_thrust,
        sf_helper_object.m_veff, sf_helper_object.m_mu, sf_helper_object.m_cut);

    auto retval = sf.compute_mismatch_constraints();
    auto retval_alpha = sf_alpha.compute_mismatch_constraints();

    std::array<double, 3> r1 = {retval[0], retval[1], retval[2]};
    std::array<double, 3> r2 = {retval_alpha[0], retval_alpha[1], retval_alpha[2]};
    std::array<double, 3> v1 = {retval[3], retval[4], retval[5]};
    std::array<double, 3> v2 = {retval_alpha[3], retval_alpha[4], retval_alpha[5]};

    REQUIRE(kep3_tests::floating_point_error_vector(r1, r2) < 1e-14);
    REQUIRE(kep3_tests::floating_point_error_vector(v1, v2) < 1e-14);
    REQUIRE(std::abs((retval[6] - retval_alpha[6]) / retval[6]) < 1e-14);
}

TEST_CASE("compute_state_history")
{
    // Get state history
    kep3::leg::sims_flanagan_hf_alpha sf{};
    auto mc = sf.compute_mismatch_constraints();
    unsigned grid_points_per_segment = 4;
    auto state_history = sf.get_state_history(grid_points_per_segment);

    // Get fwd final state
    std::vector<double> fwd_seg_sh = state_history.at(sf.get_nseg_fwd() - 1);
    std::array<double, 7> final_fwd_state{};
    std::copy(fwd_seg_sh.begin() + (grid_points_per_segment - 1) * 7l,
              fwd_seg_sh.begin() + grid_points_per_segment * 7l, final_fwd_state.begin());

    // Get bck final state
    std::vector<double> bck_seg_sh = state_history.at(sf.get_nseg_fwd());
    std::array<double, 7> final_bck_state{};
    std::copy(bck_seg_sh.begin() + (grid_points_per_segment - 1) * 7l,
              bck_seg_sh.begin() + grid_points_per_segment * 7l, final_bck_state.begin());

    // Get mismatch and calculate Linfty norm
    std::transform(final_fwd_state.begin(), final_fwd_state.end(), final_bck_state.begin(), final_fwd_state.begin(),
                   std::minus<double>());
    std::array<double, 7> manual_mismatch = final_fwd_state; // final_fwd_state is overridden with the subtracted values
    REQUIRE(kep3_tests::L_infinity_norm(manual_mismatch, mc) < 1e-15);
}

TEST_CASE("compute_state_history_2")
{
    // Initialise unique test quantities
    std::vector<double> throttles
        = {0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24,
           0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34};
    double cut = 0.6;
    // Initialise helper quantities
    auto sf_test_object = sf_hf_test_alpha_object(throttles, cut);

    kep3::leg::sims_flanagan_hf_alpha sf(
        sf_test_object.m_rvs, sf_test_object.m_ms, sf_test_object.m_throttles,
        {sf_test_object.m_tof / 10, sf_test_object.m_tof / 10, sf_test_object.m_tof / 10, sf_test_object.m_tof / 10,
         sf_test_object.m_tof / 10, sf_test_object.m_tof / 10, sf_test_object.m_tof / 10, sf_test_object.m_tof / 10,
         sf_test_object.m_tof / 10, sf_test_object.m_tof / 10},
        sf_test_object.m_rvf, sf_test_object.m_mf, sf_test_object.m_tof, sf_test_object.m_max_thrust,
        sf_test_object.m_veff, sf_test_object.m_mu, sf_test_object.m_cut, 1e-16);

    // Get state history
    auto mc = sf.compute_mismatch_constraints();
    unsigned grid_points_per_segment = 4;
    auto state_history = sf.get_state_history(grid_points_per_segment);

    // Get fwd final state
    std::vector<double> fwd_seg_sh = state_history.at(sf.get_nseg_fwd() - 1);
    std::array<double, 7> final_fwd_state{};
    std::copy(fwd_seg_sh.begin() + (grid_points_per_segment - 1) * 7l,
              fwd_seg_sh.begin() + grid_points_per_segment * 7l, final_fwd_state.begin());

    // Get bck final state
    std::vector<double> bck_seg_sh = state_history.at(sf.get_nseg_fwd());
    std::array<double, 7> final_bck_state{};
    std::copy(bck_seg_sh.begin() + (grid_points_per_segment - 1) * 7l,
              bck_seg_sh.begin() + grid_points_per_segment * 7l, final_bck_state.begin());

    // Get mismatch and calculate Linfty norm
    std::transform(final_fwd_state.begin(), final_fwd_state.end(), final_bck_state.begin(), final_fwd_state.begin(),
                   std::minus<double>());
    std::array<double, 7> manual_mismatch = final_fwd_state; // final_fwd_state is overridden with the subtracted values
    REQUIRE(kep3_tests::L_infinity_norm(manual_mismatch, mc) < 1e-15);
}

TEST_CASE("serialization_test")
{
    // Instantiate a generic lambert problem
    std::array<std::array<double, 3>, 2> rvs{{{-1, -1, -1}, {-1, -1, -1}}};
    std::array<std::array<double, 3>, 2> rvf{{{0.1, 1.1, 0.1}, {-1.1, 0.1, 0.1}}};
    kep3::leg::sims_flanagan_hf_alpha sf1{rvs, 12., {1, 2, 3, 4, 5, 6}, {0.1, 0.1}, rvf, 10, 2.3, 2.3, 2.3, 1.1, 0.2};

    // Store the string representation.
    std::stringstream ss;
    auto before = boost::lexical_cast<std::string>(sf1);
    // Now serialize
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << sf1;
    }
    // Deserialize
    // Create a new lambert problem object
    kep3::leg::sims_flanagan_hf_alpha sf_a{};
    {
        boost::archive::binary_iarchive iarchive(ss);
        iarchive >> sf_a;
    }
    auto after = boost::lexical_cast<std::string>(sf_a);
    // Compare the string represetation
    REQUIRE(before == after);
}