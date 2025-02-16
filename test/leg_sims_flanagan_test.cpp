// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <algorithm>
#include <stdexcept>
#include <vector>

#include <xtensor/xadapt.hpp>
#include <xtensor/xio.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/nlopt.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/lambert_problem.hpp>
#include <kep3/leg/sims_flanagan.hpp>
#include <kep3/planet.hpp>
#include <kep3/udpla/vsop2013.hpp>

#include "catch.hpp"
#include "leg_sims_flanagan_udp.hpp"

TEST_CASE("constructor")
{
    {
        // The default constructor constructs a valid leg with no mismatches.
        kep3::leg::sims_flanagan sf{};
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
            kep3::leg::sims_flanagan(rvs, ms, {0., 0., 0., 0., 0., 0.}, rvf, mf, kep3::pi / 2, 1., 1., 1., 0.5));
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
        REQUIRE_THROWS_AS(kep3::leg::sims_flanagan(rvs, ms, {}, rvf, mf, kep3::pi / 2, 1., 1., 1., 0.5),
                          std::logic_error);
    }
}

TEST_CASE("getters_and_setters")
{
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
    {
        kep3::leg::sims_flanagan sf{};
        std::array<std::array<double, 3>, 2> rvf{{{1, 1, 1}, {1, 1, 1}}};
        std::vector<double> throttles{1., 2., 3., 1., 2., 3.};

        sf.set(rvf, 12, throttles, rvf, 12, 4, 4, 4, 4, 0.333);
        REQUIRE(sf.get_rvs() == rvf);
        REQUIRE(sf.get_ms() == 12);
        REQUIRE(sf.get_rvf() == rvf);
        REQUIRE(sf.get_mf() == 12);
        REQUIRE(sf.get_throttles() == throttles);
        REQUIRE(sf.get_max_thrust() == 4);
        REQUIRE(sf.get_isp() == 4);
        REQUIRE(sf.get_mu() == 4);
        REQUIRE(sf.get_tof() == 4);
        REQUIRE(sf.get_cut() == 0.333);
    }
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
//
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
            REQUIRE(*std::max_element(mc.begin(), mc.end()) < 1e-8);
        }
    }

    {
        // Here we reuse the ballitic arc as a ground truth for an optimization.
        // We check that, when feasible, the optimal mass solution is indeed ballistic.
        pagmo::problem prob{sf_test_udp{rv0, mass, rv1, 0.05, 2000, 10u}};
        prob.set_c_tol(1e-8);
        bool found = false;
        unsigned trial = 0u;
        pagmo::nlopt uda{"slsqp"};
        uda.set_xtol_abs(1e-10);
        uda.set_xtol_rel(1e-10);
        uda.set_ftol_abs(0);
        uda.set_maxeval(1000);
        pagmo::algorithm algo{uda};
        while ((!found) && (trial < 20u)) {
            pagmo::population pop{prob, 1u};
            algo.set_verbosity(10u);
            pop = algo.evolve(pop);
            auto champ = pop.champion_f();
            found = prob.feasibility_f(champ);
            if (found) {
                fmt::print("{}\n", champ);
                found = *std::min_element(champ.begin() + 7, champ.end()) < -0.99999;
            }
            trial++;
        }
        REQUIRE_FALSE(!found); // If this does not pass, then the optimization above never found a ballistic arc ...
                               // theres a problem somewhere.
    }
    {
        // Here we create an ALMOST ballistic arc as a ground truth for an optimization.
        // We check that, when feasible, the optimal mass solution is indeed ballistic.
        auto rv1_modified = rv1;
        rv1_modified[1][0] += 1000; // Adding 1km/s along x
        pagmo::problem prob{sf_test_udp{rv0, mass, rv1_modified, 0.05, 2000, 10u}};
        prob.set_c_tol(1e-8);
        bool found = false;
        unsigned trial = 0u;
        pagmo::nlopt uda{"slsqp"};
        uda.set_xtol_abs(1e-10);
        uda.set_xtol_rel(1e-10);
        uda.set_ftol_abs(0);
        uda.set_maxeval(1000);
        pagmo::algorithm algo{uda};
        while ((!found) && (trial < 20u)) {
            pagmo::population pop{prob, 1u};
            algo.set_verbosity(10u);
            pop = algo.evolve(pop);
            auto champ = pop.champion_f();
            found = prob.feasibility_f(champ);
            if (found) {
                fmt::print("{}\n", champ);
            }
            trial++;
        }
        // If this does not pass, then the optimization above never converged to a feasible solution.
        REQUIRE_FALSE(!found);
    }
}

TEST_CASE("mismatch_constraints_test2")
{
    // The ground truth were computed with these values of the astro constants, hence we cannot use here the pykep ones.
    double AU_OLD = 149597870700.0;
    double EV_OLD = 29784.691831696804;
    double MU_OLD = 1.32712440018e20;
    // We test the correctness of the compute_mismatch_constraints computations against a ground truth (computed with a
    // different program)
    std::array<std::array<double, 3>, 2> rvs{
        {{1 * AU_OLD, 0.1 * AU_OLD, -0.1 * AU_OLD}, {0.2 * EV_OLD, 1 * EV_OLD, -0.2 * EV_OLD}}};

    std::array<std::array<double, 3>, 2> rvf{
        {{1.2 * AU_OLD, -0.1 * AU_OLD, 0.1 * AU_OLD}, {-0.2 * EV_OLD, 1.023 * EV_OLD, -0.44 * EV_OLD}}};

    double ms = 1500.;
    double mf = 1300.;
    std::vector<double> throttles
        = {0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24};
    kep3::leg::sims_flanagan sf(rvs, ms, throttles, rvf, mf, 324.0 * kep3::DAY2SEC, 0.12, 100, MU_OLD, 0.6);
    auto retval = sf.compute_mismatch_constraints();
    std::vector<double> ground_truth
        = {-1.9701274809621304e+11, 4.6965044246848071e+11, -1.5007523306033661e+11, -2.9975151466948650e+04,
           -2.8264916164742666e+04, 1.0264806797549732e+04, -8.2807673427721159e+02};
    REQUIRE(std::abs((retval[0] - ground_truth[0]) / retval[0]) < 1e-13);
    REQUIRE(std::abs((retval[1] - ground_truth[1]) / retval[1]) < 1e-13);
    REQUIRE(std::abs((retval[2] - ground_truth[2]) / retval[2]) < 1e-13);
    REQUIRE(std::abs((retval[3] - ground_truth[3]) / retval[3]) < 1e-13);
    REQUIRE(std::abs((retval[4] - ground_truth[4]) / retval[4]) < 1e-13);
    REQUIRE(std::abs((retval[5] - ground_truth[5]) / retval[5]) < 1e-13);
}

TEST_CASE("grad_test")
{
    // Here we test the analytical gradient against an equivalent numerical one. We do so through the udp "sf_test_udp"
    std::array<std::array<double, 3>, 2> rvs{
        {{1 * kep3::AU, 0.1 * kep3::AU, -0.1 * kep3::AU},
         {0.2 * kep3::EARTH_VELOCITY, 1 * kep3::EARTH_VELOCITY, -0.2 * kep3::EARTH_VELOCITY}}};
    //
    std::array<std::array<double, 3>, 2> rvf{
        {{1.2 * kep3::AU, -0.1 * kep3::AU, 0.1 * kep3::AU},
         {-0.2 * kep3::EARTH_VELOCITY, 1.023 * kep3::EARTH_VELOCITY, -0.44 * kep3::EARTH_VELOCITY}}};
    //
    double ms = 1500.;

    sf_test_udp udp{rvs, ms, rvf, 0.12, 100, 5};
    std::vector<double> x
        = {0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 324., 1300.};
    auto grad = udp.gradient_numerical(x);
    auto grad_a = udp.gradient(x);
    auto xgrad = xt::adapt(grad, {1u + 7u + 5u, 17u});
    auto xgrad_a = xt::adapt(grad_a, {1u + 7u + 5u, 17u});
    REQUIRE(xt::linalg::norm(xgrad - xgrad_a)
            < 1e-8); // With the high fidelity gradient this is still the best we can achieve
}

TEST_CASE("serialization_test")
{
    // Instantiate a generic lambert problem
    std::array<std::array<double, 3>, 2> rvs{{{-1, -1, -1}, {-1, -1, -1}}};
    std::array<std::array<double, 3>, 2> rvf{{{0.1, 1.1, 0.1}, {-1.1, 0.1, 0.1}}};
    kep3::leg::sims_flanagan sf1{rvs, 12., {1, 2, 3, 4, 5, 6}, rvf, 10, 2.3, 2.3, 2.3, 1.1, 0.2};

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
    kep3::leg::sims_flanagan sf2{};
    {
        boost::archive::binary_iarchive iarchive(ss);
        iarchive >> sf2;
    }
    auto after = boost::lexical_cast<std::string>(sf2);
    // Compare the string represetation
    REQUIRE(before == after);
}
