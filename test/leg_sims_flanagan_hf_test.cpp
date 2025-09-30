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
#include <stdexcept>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/io/xio.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/nlopt.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/lambert_problem.hpp>
#include <kep3/leg/sf_checks.hpp>
#include <kep3/leg/sims_flanagan.hpp>
#include <kep3/leg/sims_flanagan_hf.hpp>
#include <kep3/planet.hpp>
#include <kep3/ta/zero_hold_kep.hpp>
#include <kep3/udpla/vsop2013.hpp>

#include "catch.hpp"
#include "leg_sims_flanagan_hf_helpers.hpp"
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

#include "leg_sims_flanagan_hf_udp.hpp"

TEST_CASE("constructor")
{
    {
        // The default constructor constructs a valid leg with no mismatches.
        kep3::leg::sims_flanagan_hf sf{};
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
        REQUIRE_NOTHROW(
            kep3::leg::sims_flanagan_hf(rvs, ms, {0., 0., 0., 0., 0., 0.}, rvf, mf, kep3::pi / 2, 1., 1., 1., 0.5));
        const std::array<double, 7> rvms{1, 0, 0, 0, 1, 0, 1};
        const std::array<double, 7> rvmf{0, 1, 0, -1, 0, 0, 1};
        REQUIRE_NOTHROW(
            kep3::leg::sims_flanagan_hf(rvms, {0., 0., 0., 0., 0., 0.}, rvmf, kep3::pi / 2, 1., 1., 1., 0.5));
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan_hf(rvs, ms, {0., 0., 0., 0., 0.}, rvf, mf, kep3::pi / 2, 1., 1., 1., 0.5),
            std::logic_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan_hf(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, -0.42, 1., 1., 1., 0.5),
                          std::domain_error);
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan_hf(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, kep3::pi / 2, -0.3, 1., 1., 0.5),
            std::domain_error);
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan_hf(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, kep3::pi / 2, 1., -2., 1., 0.5),
            std::domain_error);
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan_hf(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, kep3::pi / 2, 1., 1., -0.32, 0.5),
            std::domain_error);
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan_hf(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, kep3::pi / 2, 1., 1., 1., 32),
            std::domain_error);
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan_hf(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, kep3::pi / 2, 1., 1., 1., -0.1),
            std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan_hf(rvs, ms, {}, rvf, mf, kep3::pi / 2, 1., 1., 1., 0.5),
                          std::logic_error);
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan_hf(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, kep3::pi / 2, 1., 1., 1., 0.5, -1e-2),
            std::domain_error);
        REQUIRE_THROWS_AS(
            kep3::leg::sims_flanagan_hf(rvs, ms, {0, 0, 0, 0, 0, 0}, rvf, mf, kep3::pi / 2, 1., 1., 1., 0.5, 1.2),
            std::domain_error);
        REQUIRE_THROWS_AS(kep3::leg::_check_nseg(2, 1, 2), std::logic_error);
    }
}

TEST_CASE("getters_and_setters")
{
    {
        kep3::leg::sims_flanagan_hf sf{};
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
        kep3::leg::sims_flanagan_hf sf{};
        std::array<std::array<double, 3>, 2> rvf{{{1, 1, 1}, {1, 1, 1}}};
        std::vector<double> throttles{1., 2., 3., 1., 2., 3.};

        sf.set(rvf, 12, throttles, rvf, 12, 4, 4, 4, 4, 0.333, 2e-5);
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
        kep3::leg::sims_flanagan_hf sf{};
        std::array<double, 7> rvms{1, 1, 1, 1, 1, 1, 1};
        std::vector<double> throttles{1., 2., 3., 1., 2., 3.};

        sf.set(rvms, throttles, rvms, 4, 4, 4, 4, 0.333, 2e-5);
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
}

TEST_CASE("compute_throttle_constraints_test")
{
    std::array<std::array<double, 3>, 2> rvs{{{1, 0, 0}, {0, 1, 0}}};
    std::array<std::array<double, 3>, 2> rvf{{{0, 1, 0}, {-1, 0, 0}}};
    kep3::leg::sims_flanagan_hf sf(rvs, 1., {0, 1, 0, 1, 1, 1, 0, 1, 1}, rvf, 1, 1, 1, 1, 1, 1);
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
                kep3::leg::sims_flanagan_hf sf(rv0, 1., throttles, rv1, 1., dt, 1., 1., kep3::MU_SUN, cut);
                auto mc = sf.compute_mismatch_constraints();
                mc = normalize_con(mc);
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
        double max_thrust = 0.05;
        double veff = 2500.*kep3::G0;
        kep3::leg::sims_flanagan_hf sf(rv0, m0, throttles, rv1, m1, tof, max_thrust, veff, 1.32712440018e+20, 0.5);
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
    double veff = 3000.0 * kep3::G0;

    // Define optimization parameters
    size_t nseg = 8u;
    std::vector<double> lb(nseg * 3 + 2, 0.0);
    lb[nseg * 3] = 500.0;          // days
    lb[nseg * 3 + 1] = mass / 2.0; // kg

    std::vector<double> throttles(nseg * 3, 0.0);
    std::vector<double> talphas(nseg, dt / static_cast<double>(nseg));

    pagmo::problem prob{sf_hf_test_udp{rv0, mass, rv1, max_thrust, veff, 8u}};

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
            kep3::leg::sims_flanagan_hf sf(rv0, 1., throttles, rv1, 1., dt, 1., 1., kep3::MU_SUN, cut);
            auto mc = sf.compute_mismatch_constraints();
            mc = normalize_con(mc);
            REQUIRE(*std::max_element(mc.begin(), mc.end()) < 1e-8);
        }
    }

    {
        // Here we reuse the ballitic arc as a ground truth for an optimization.
        // We check that, when feasible, the optimal mass solution is indeed ballistic.
        pagmo::problem prob{sf_hf_test_udp{rv0, mass, rv1, max_thrust, veff, 8u}};
        prob.set_c_tol(1e-6);
        bool found = false;
        pagmo::nlopt uda{"slsqp"};
        uda.set_xtol_abs(1e-8);
        uda.set_xtol_rel(1e-8);
        uda.set_ftol_abs(0);
        uda.set_maxeval(1000);
        pagmo::algorithm algo{uda};

        std::vector<double> x0 = {2.981309716881921e-05,   4.386247528314636e-06,
                                  -8.386247528314636e-06,  -1.0489297915311248e-05,
                                  -2.75863662269476e-06,   -2.0057791800205834e-05,
                                  -2.0168670164500432e-05, -1.0168670164500432e-05,
                                  2.6796831660229594e-05,  -1.981309716881921e-05,
                                  1.3984567008747744e-06,  4.485693307454289e-06,
                                  1.323815701173753e-05,   1.7142463462947175e-06,
                                  3.70831124044131e-08,    -3.6176089862227655e-06,
                                  5.9303303176247626e-06,  -1.0678752064652906e-05,
                                  1.5370329134036994e-05,  1.7881716332279772e-06,
                                  1.183633361556691e-05,   -9.323815701173753e-07,
                                  -5.046909298335364e-06,  -1.3506155624149922e-05,
                                  649.9992117409238,       1230};
        pagmo::population pop{prob};
        pop.push_back(x0);
        algo.set_verbosity(10u);
        pop = algo.evolve(pop);
        auto champ = pop.champion_f();
        found = prob.feasibility_f(champ);
        fmt::print("{}\n", champ);
        fmt::print("{}\n", pop.champion_x());
        REQUIRE_FALSE(!found); // If this does not pass, then the optimization above never found a ballistic arc ...
                               // theres a problem somewhere.
    }
}

// Compare low-fidelity and high-fidelity methods with zero thrust (ought to be the same)
TEST_CASE("compare_low_and_high_fidelity")
{

    // Initialise unique test quantities
    double cut = 0.6;
    auto sf_helper_object = sf_hf_test_object(cut);

    kep3::leg::sims_flanagan_hf sf(sf_helper_object.m_rvs, sf_helper_object.m_ms, sf_helper_object.m_throttles,
                                   sf_helper_object.m_rvf, sf_helper_object.m_mf, sf_helper_object.m_tof,
                                   sf_helper_object.m_max_thrust, sf_helper_object.m_veff, sf_helper_object.m_mu,
                                   sf_helper_object.m_cut, 1e-16);
    kep3::leg::sims_flanagan sf_lf(sf_helper_object.m_rvs, sf_helper_object.m_ms, sf_helper_object.m_throttles,
                                   sf_helper_object.m_rvf, sf_helper_object.m_mf, sf_helper_object.m_tof,
                                   sf_helper_object.m_max_thrust, sf_helper_object.m_veff, sf_helper_object.m_mu,
                                   sf_helper_object.m_cut);

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

// Compare high-fidelity method with manually calculated (direct heyoka interfacing) Taylor integration.
TEST_CASE("compute_mismatch_constraints_test3")
{

    // Initialise unique test quantities
    double cut = 0.5;
    std::vector<double> throttles = {0.10, 0.11, 0.12, 0.13, 0.14, 0.15};
    auto sf_test_object = sf_hf_test_object(throttles, cut);
    std::array<double, 7> mc_manual = sf_test_object.compute_manual_mc();

    // Calculate equivalent with hf leg.
    kep3::leg::sims_flanagan_hf sf(sf_test_object.m_rvs, sf_test_object.m_ms, sf_test_object.m_throttles,
                                   sf_test_object.m_rvf, sf_test_object.m_mf, sf_test_object.m_tof,
                                   sf_test_object.m_max_thrust, sf_test_object.m_veff, sf_test_object.m_mu,
                                   sf_test_object.m_cut, 1e-16);
    auto mc_sf_hf = sf.compute_mismatch_constraints();

    std::array<double, 3> r1 = {mc_sf_hf[0], mc_sf_hf[1], mc_sf_hf[2]};
    std::array<double, 3> r2 = {mc_manual[0], mc_manual[1], mc_manual[2]};
    REQUIRE(kep3_tests::floating_point_error_vector(r1, r2) < 1e-16);
    std::array<double, 3> v1 = {mc_sf_hf[3], mc_sf_hf[4], mc_sf_hf[5]};
    std::array<double, 3> v2 = {mc_manual[3], mc_manual[4], mc_manual[5]};
    REQUIRE(kep3_tests::floating_point_error_vector(v1, v2) < 1e-16);
    REQUIRE(std::abs((mc_sf_hf[6] - mc_manual[6]) / mc_sf_hf[6]) < 1e-16);
}

TEST_CASE("compute_mc_grad_test")
{
    // Initialise unique test quantities
    std::vector<double> throttles_full
        = {0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24,
           0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34};
    std::array<unsigned long, 4> nseg_array = {1, 2, 5, 10};
    std::array<double, 7> cut_array = {0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
    std::vector<double> throttles;
    for (unsigned long nseg : nseg_array) {
        for (double cut : cut_array) {
            throttles
                = std::vector<double>(throttles_full.begin(), throttles_full.begin() + static_cast<long>(nseg) * 3);
            auto sf_test_object = sf_hf_test_object(throttles, cut);

            // Numerical gradient
            std::vector<double> num_grad = sf_test_object.compute_numerical_gradient();
            auto xt_num_gradients = xt::adapt(num_grad, {7u + nseg, 7u + 3u * nseg + 1u + 7u});
            auto xt_num_mc_gradients = xt::view(xt_num_gradients, xt::range(0, 7), xt::all());

            // Analytical gradient
            std::vector<double> a_gradients = sf_test_object.compute_analytical_gradient();
            auto xt_a_gradients = xt::adapt(a_gradients, {7u, 7u + 3u * static_cast<unsigned>(nseg) + 1u + 7u});

            REQUIRE(xt::linalg::norm(xt_num_mc_gradients - xt_a_gradients)
                    < 1e-8); // With the high fidelity gradient this is still the best we can achieve. The difference is
                             // like 4.56e-8
        }
    }
}

TEST_CASE("compute_tc_grad_test")
{

    // Initialise unique test quantities
    std::vector<double> throttles
        = {0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24};
    unsigned int nseg = static_cast<unsigned int>(throttles.size()) / 3;
    double cut = 0.6;
    // Initialise helper quantities
    auto sf_test_object = sf_hf_test_object(throttles, cut);

    // Numerical gradient
    std::vector<double> num_grad = sf_test_object.compute_numerical_gradient();
    auto xt_num_gradients = xt::adapt(num_grad, {7u + nseg, 30u});
    auto xt_num_tc_gradients = xt::view(xt_num_gradients, xt::range(7, 12), xt::range(7, 22));

    // Calculate throttle constraint gradients
    kep3::leg::sims_flanagan_hf sf(sf_test_object.m_rvs, sf_test_object.m_ms, sf_test_object.m_throttles,
                                   sf_test_object.m_rvf, sf_test_object.m_mf, sf_test_object.m_tof,
                                   sf_test_object.m_max_thrust, sf_test_object.m_veff, sf_test_object.m_mu,
                                   sf_test_object.m_cut, 1e-16);
    std::vector<double> tc_a_grad = sf.compute_tc_grad();
    auto xt_tc_a_grad = xt::adapt(tc_a_grad, {nseg, 3u * nseg});

    REQUIRE(xt::linalg::norm(xt_num_tc_gradients - xt_tc_a_grad) < 1e-13); // 1e-14 fails
}

TEST_CASE("compute_state_history")
{
    // Get state history
    kep3::leg::sims_flanagan_hf sf{};
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
    auto sf_test_object = sf_hf_test_object(throttles, cut);

    kep3::leg::sims_flanagan_hf sf(sf_test_object.m_rvs, sf_test_object.m_ms, sf_test_object.m_throttles,
                                   sf_test_object.m_rvf, sf_test_object.m_mf, sf_test_object.m_tof,
                                   sf_test_object.m_max_thrust, sf_test_object.m_veff, sf_test_object.m_mu,
                                   sf_test_object.m_cut, 1e-16);

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
    kep3::leg::sims_flanagan_hf sf1{rvs, 12., {1, 2, 3, 4, 5, 6}, rvf, 10, 2.3, 2.3, 2.3, 1.1, 0.2};

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
    kep3::leg::sims_flanagan_hf sf_a{};
    {
        boost::archive::binary_iarchive iarchive(ss);
        iarchive >> sf_a;
    }
    auto after = boost::lexical_cast<std::string>(sf_a);
    // Compare the string represetation
    REQUIRE(before == after);
}