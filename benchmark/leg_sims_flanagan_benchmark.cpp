// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <chrono>
#include <iostream>
#include <random>

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

#include "leg_sims_flanagan_udp_bench.hpp"

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void perform_test_speed(unsigned N, unsigned nseg, unsigned pop_size)
{
    //
    // Engines
    //
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(122012203u);
    //
    // Distributions
    //
    std::uniform_real_distribution<double> dv_pert_d(0., 1000.);
    std::uniform_real_distribution<double> mass_d(500, 1500);
    std::uniform_real_distribution<double> tof_d(1000, 1500);
    std::uniform_real_distribution<double> ts_d(1100, 1300);

    // We construct the solver
    pagmo::nlopt uda{"slsqp"};
    uda.set_xtol_abs(1e-8);
    uda.set_xtol_rel(1e-8);
    uda.set_ftol_abs(1e-8);
    uda.set_maxeval(1000);
    pagmo::algorithm algo{uda};
    algo.set_verbosity(0u);

    // The initial positions
    kep3::udpla::vsop2013 udpla_earth("earth_moon", 1e-2);
    kep3::udpla::vsop2013 udpla_jupiter("jupiter", 1e-2);
    kep3::planet earth{udpla_earth};
    kep3::planet jupiter{udpla_jupiter};
    double count_a = 0;
    double count_n = 0;
    std::cout << "\n";
    for (auto i = 0u; i < N; ++i) {
        // And some epochs / tofs.
        const double tof_days = tof_d(rng_engine);
        const double tof = tof_days * kep3::DAY2SEC;
        const double ts = ts_d(rng_engine);
        const double mass = mass_d(rng_engine);
        auto rvs = earth.eph(ts);
        auto rvf = jupiter.eph(ts + tof_days);
        // We create a ballistic arc matching the two.
        const kep3::lambert_problem lp{rvs[0], rvf[0], tof, kep3::MU_SUN};
        rvs[1][0] = lp.get_v0()[0][0];
        rvs[1][1] = lp.get_v0()[0][1];
        rvs[1][2] = lp.get_v0()[0][2];
        rvf[1][0] = lp.get_v1()[0][0] + dv_pert_d(rng_engine);
        rvf[1][1] = lp.get_v1()[0][1] + dv_pert_d(rng_engine);
        rvf[1][2] = lp.get_v1()[0][2] + dv_pert_d(rng_engine);

        // We construct two problems (analytical gradient and numerical gradient)
        pagmo::problem prob_a{sf_bench_udp{rvs, mass, rvf, 0.05, 2000, nseg, true}};
        pagmo::problem prob_n{sf_bench_udp{rvs, mass, rvf, 0.05, 2000, nseg, false}};
        prob_a.set_c_tol(1e-8);
        prob_n.set_c_tol(1e-8);

        // We construct the random chromosmes 
        const pagmo::population pop{prob_a, pop_size};

        // First we time the analytical gradients
        auto start = high_resolution_clock::now();
        for (decltype(pop_size) j = 0u; j < pop_size; ++j) {
            prob_a.gradient(pop.get_x()[j]);
        }
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        count_a+=static_cast<double>(duration.count()) / 1e6;

        // then the numerical ones
        start = high_resolution_clock::now();
        for (decltype(pop_size) j = 0u; j < pop_size; ++j) {
            prob_n.gradient(pop.get_x()[j]);
        }
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        count_n+=static_cast<double>(duration.count()) / 1e6;
    }
    fmt::print("{} nseg - timing: analytical {} - numerical {}", nseg, count_a, count_n);
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void perform_test_convergence(unsigned N, unsigned nseg)
{
    //
    // Engines
    //
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(122012203u);
    //
    // Distributions
    //
    std::uniform_real_distribution<double> dv_pert_d(0., 1000.);
    std::uniform_real_distribution<double> mass_d(500, 1500);
    std::uniform_real_distribution<double> tof_d(1000, 1500);
    std::uniform_real_distribution<double> ts_d(1100, 1300);

    // We construct the solver
    pagmo::nlopt uda{"slsqp"};
    uda.set_xtol_abs(0);
    uda.set_xtol_rel(0);
    uda.set_ftol_abs(0);
    uda.set_maxeval(1000);
    pagmo::algorithm algo{uda};
    algo.set_verbosity(0u);

    // The initial positions
    kep3::udpla::vsop2013 udpla_earth("earth_moon", 1e-2);
    kep3::udpla::vsop2013 udpla_jupiter("jupiter", 1e-2);
    kep3::planet earth{udpla_earth};
    kep3::planet jupiter{udpla_jupiter};
    unsigned count_a = 0;
    unsigned count_n = 0;
    std::cout << "\n";
    for (auto i = 0u; i < N; ++i) {
        // And some epochs / tofs.
        const double tof_days = tof_d(rng_engine);
        const double tof = tof_days * kep3::DAY2SEC;
        double ts = ts_d(rng_engine);
        const double mass = mass_d(rng_engine);
        auto rvs = earth.eph(ts);
        auto rvf = jupiter.eph(ts + tof_days);
        // We create a ballistic arc matching the two.
        const kep3::lambert_problem lp{rvs[0], rvf[0], tof, kep3::MU_SUN};
        rvs[1][0] = lp.get_v0()[0][0];
        rvs[1][1] = lp.get_v0()[0][1];
        rvs[1][2] = lp.get_v0()[0][2];
        rvf[1][0] = lp.get_v1()[0][0] + dv_pert_d(rng_engine);
        rvf[1][1] = lp.get_v1()[0][1] + dv_pert_d(rng_engine);
        rvf[1][2] = lp.get_v1()[0][2] + dv_pert_d(rng_engine);

        // We construct two problems (analytical gradient and numerical gradient)
        pagmo::problem prob_a{sf_bench_udp{rvs, mass, rvf, 0.05, 2000, nseg, true}};
        pagmo::problem prob_n{sf_bench_udp{rvs, mass, rvf, 0.05, 2000, nseg, false}};
        prob_a.set_c_tol(1e-8);
        prob_n.set_c_tol(1e-8);

        // We construct the same random population
        pagmo::population pop_a{prob_a, 1u};
        pagmo::population pop_n{prob_n};
        pop_n.push_back(pop_a.get_x()[0]);

        // We solve first a
        pop_a = algo.evolve(pop_a);
        if (prob_a.feasibility_f(pop_a.get_f()[0])) {
            count_a++;
            std::cout << "." << std::flush;
        } else {
            std::cout << "x" << std::flush;
        }
        // then n
        pop_n = algo.evolve(pop_n);
        if (prob_n.feasibility_f(pop_n.get_f()[0])) {
            count_n++;
            std::cout << "." << std::flush;
        } else {
            std::cout << "x" << std::flush;
        }
    }
    fmt::print("\n{} nseg - success rates: analytical {}/{} - numerical {}/{}", nseg, count_a, N, count_n, N);
}

int main()
{
    // performing tests
    fmt::print("\nSolves the same optimization problems with and without analytical gradients:");
    perform_test_convergence(200, 5);
    perform_test_convergence(200, 10);
    perform_test_convergence(200, 15);

    fmt::print("\nComputes the same analytical and numerical gradients and tests for speed:");
    perform_test_speed(100, 5, 10);
    perform_test_speed(100, 10, 10);
    perform_test_speed(100, 15, 10);
    perform_test_speed(100, 20, 10);
    perform_test_speed(100, 70, 10);
    fmt::print("\n");

}