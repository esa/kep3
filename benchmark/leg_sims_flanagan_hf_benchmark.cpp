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
#include <kep3/leg/sims_flanagan_hf.hpp>
#include <kep3/planet.hpp>
#include <kep3/udpla/vsop2013.hpp>

#include "leg_sims_flanagan_hf_udp_bench.hpp"

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;

void perform_single_nogradient_speed_test()
{
    std::array<std::array<double, 3>, 2> m_rvs{{{1, 0.1, -0.1}, {0.2, 1, -0.2}}};
    std::array<std::array<double, 3>, 2> m_rvf{{{1.2, -0.1, 0.1}, {-0.2, 1.023, -0.44}}};
    double m_ms = 1;
    double m_mf = m_ms * 13 / 15;
    double m_isp = 1;
    double m_max_thrust = 1;
    double m_cut = 0.5;
    double m_mu = 1;
    double m_tof = 1;
    std::vector<double> m_throttles = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

    auto start_con = high_resolution_clock::now();
    auto sf_leg
        = kep3::leg::sims_flanagan(m_rvs, m_ms, m_throttles, m_rvf, m_mf, m_tof, m_max_thrust, m_isp, m_mu, m_cut);
    auto stop_con = high_resolution_clock::now();
    auto duration_con = duration_cast<microseconds>(stop_con - start_con);
    fmt::print("\nLow-fidelity leg construction: {} nseg - timing: {}", m_throttles.size() / 3,
               static_cast<double>(duration_con.count()) / 1e6);
    auto start = high_resolution_clock::now();
    auto mc = sf_leg.compute_mismatch_constraints();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    fmt::print("\nLow-fidelity leg mc: {} nseg - timing: {}", m_throttles.size() / 3,
               static_cast<double>(duration.count()) / 1e6);
    auto two_start = high_resolution_clock::now();
    auto two_mc = sf_leg.compute_mc_grad();
    auto two_stop = high_resolution_clock::now();
    auto two_duration = duration_cast<microseconds>(two_stop - two_start);
    fmt::print("\nLow-fidelity leg mc_grad: {} nseg - timing: {}", m_throttles.size() / 3,
               static_cast<double>(two_duration.count()) / 1e6);

    auto start_hf_con = high_resolution_clock::now();
    auto sf_hf_leg = kep3::leg::sims_flanagan_hf(m_rvs, m_ms, m_throttles, m_rvf, m_mf, m_tof, m_max_thrust, m_isp,
                                                 m_mu, m_cut, 1e-16);
    auto stop_hf_con = high_resolution_clock::now();
    auto duration_hf_con = duration_cast<microseconds>(stop_hf_con - start_hf_con);
    fmt::print("\nHigh-fidelity leg construction: {} nseg - timing: {}", m_throttles.size() / 3,
               static_cast<double>(duration_hf_con.count()) / 1e6);
    auto hf_start = high_resolution_clock::now();
    auto hf_mc = sf_hf_leg.compute_mismatch_constraints();
    auto hf_stop = high_resolution_clock::now();
    auto hf_duration = duration_cast<microseconds>(hf_stop - hf_start);
    fmt::print("\nHigh-fidelity leg mc: {} nseg - timing: {}", m_throttles.size() / 3,
               static_cast<double>(hf_duration.count()) / 1e6);
    auto hf_two_start = high_resolution_clock::now();
    auto hf_two_mc = sf_hf_leg.compute_mc_grad();
    auto hf_two_stop = high_resolution_clock::now();
    auto hf_two_duration = duration_cast<microseconds>(hf_two_stop - hf_two_start);
    fmt::print("\nHigh-fidelity leg mc_grad: {} nseg - timing: {}", m_throttles.size() / 3,
               static_cast<double>(hf_two_duration.count()) / 1e6);
}

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
    std::uniform_real_distribution<double> dv_pert_d(0., 0.1);
    std::uniform_real_distribution<double> mass_d(0.9, 1.1);
    std::uniform_real_distribution<double> tof_d(0.9, 1.1);
    std::uniform_real_distribution<double> ts_d(1100, 1300);

    // We construct the solver
    pagmo::nlopt uda{"slsqp"};
    uda.set_xtol_abs(1e-8);
    uda.set_xtol_rel(1e-8);
    uda.set_ftol_abs(1e-8);
    uda.set_maxeval(1000);
    pagmo::algorithm algo{uda};
    // algo.set_verbosity(5u);

    // The initial positions
    kep3::udpla::vsop2013 udpla_earth("earth_moon", 1e-2);
    kep3::udpla::vsop2013 udpla_jupiter("jupiter", 1e-2);
    kep3::planet earth{udpla_earth};
    kep3::planet jupiter{udpla_jupiter};
    double count_a = 0;
    double count_n = 0;
    std::cout << "\n";
    auto rvs = earth.eph(1000);
    auto rvf = jupiter.eph(1000);
    double mass = 1;
    auto bench_udp_a = sf_hf_bench_udp{rvs, mass, rvf, 1, 1, nseg, true};
    auto bench_udp_n = sf_hf_bench_udp{rvs, mass, rvf, 1, 1, nseg, false};

    for (auto i = 0u; i < N; ++i) {
        // And some epochs / tofs.
        double ts = ts_d(rng_engine);
        rvs = earth.eph(ts);
        rvf = jupiter.eph(ts);
        mass = mass_d(rng_engine);
        // We create a ballistic arc matching the two.
        rvs[0][0] /= kep3::AU;
        rvs[0][1] /= kep3::AU;
        rvs[0][2] /= kep3::AU;
        rvf[0][0] /= kep3::AU;
        rvf[0][1] /= kep3::AU;
        rvf[0][2] /= kep3::AU;
        const double tof_days = tof_d(rng_engine);
        const kep3::lambert_problem lp{rvs[0], rvf[0], tof_days, 1.0};
        rvs[1][0] = lp.get_v0()[0][0] + dv_pert_d(rng_engine);
        rvs[1][1] = lp.get_v0()[0][1] + dv_pert_d(rng_engine);
        rvs[1][2] = lp.get_v0()[0][2] + dv_pert_d(rng_engine);
        rvf[1][0] = lp.get_v1()[0][0] + dv_pert_d(rng_engine);
        rvf[1][1] = lp.get_v1()[0][1] + dv_pert_d(rng_engine);
        rvf[1][2] = lp.get_v1()[0][2] + dv_pert_d(rng_engine);

        // We construct two problems (analytical gradient and numerical gradient)
        bench_udp_a.set_leg(rvs, mass, rvf, 1, 1, nseg, true);
        bench_udp_n.set_leg(rvs, mass, rvf, 1, 1, nseg, false);
        pagmo::problem prob_a{bench_udp_a};
        pagmo::problem prob_n{bench_udp_n};
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
        count_a += static_cast<double>(duration.count()) / 1e6;

        // then the numerical ones
        start = high_resolution_clock::now();
        for (decltype(pop_size) j = 0u; j < pop_size; ++j) {
            prob_n.gradient(pop.get_x()[j]);
        }
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        count_n += static_cast<double>(duration.count()) / 1e6;
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
    std::uniform_real_distribution<double> dv_pert_d(0., 0.1);
    std::uniform_real_distribution<double> mass_d(0.9, 1.1);
    std::uniform_real_distribution<double> tof_d(0.9, 1.1);
    std::uniform_real_distribution<double> ts_d(1100, 1300);

    // We construct the solver
    pagmo::nlopt uda{"slsqp"};
    uda.set_xtol_abs(0);
    uda.set_xtol_rel(0);
    uda.set_ftol_abs(0);
    uda.set_maxeval(1000);
    pagmo::algorithm algo{uda};
    // algo.set_verbosity(5u);

    // The initial positions
    kep3::udpla::vsop2013 udpla_earth("earth_moon", 1e-2);
    kep3::udpla::vsop2013 udpla_jupiter("jupiter", 1e-2);
    kep3::planet earth{udpla_earth};
    kep3::planet jupiter{udpla_jupiter};
    unsigned count_a = 0;
    unsigned count_n = 0;
    std::cout << "\n";
    // Create sf_hf_leg outside loop to save time
    auto rvs = earth.eph(1000);
    auto rvf = jupiter.eph(1000);
    double mass = 1;
    auto bench_udp_a = sf_hf_bench_udp{rvs, mass, rvf, 1, 1, nseg, true};
    auto bench_udp_n = sf_hf_bench_udp{rvs, mass, rvf, 1, 1, nseg, false};

    for (auto i = 0u; i < N; ++i) {
        // And some epochs / tofs.
        double ts = ts_d(rng_engine);
        rvs = earth.eph(ts);
        rvf = jupiter.eph(ts);
        mass = mass_d(rng_engine);
        // We create a ballistic arc matching the two.
        rvs[0][0] /= kep3::AU;
        rvs[0][1] /= kep3::AU;
        rvs[0][2] /= kep3::AU;
        rvf[0][0] /= kep3::AU;
        rvf[0][1] /= kep3::AU;
        rvf[0][2] /= kep3::AU;
        const double tof_days = tof_d(rng_engine);
        const kep3::lambert_problem lp{rvs[0], rvf[0], tof_days, 1.0};
        rvs[1][0] = lp.get_v0()[0][0] + dv_pert_d(rng_engine);
        rvs[1][1] = lp.get_v0()[0][1] + dv_pert_d(rng_engine);
        rvs[1][2] = lp.get_v0()[0][2] + dv_pert_d(rng_engine);
        rvf[1][0] = lp.get_v1()[0][0] + dv_pert_d(rng_engine);
        rvf[1][1] = lp.get_v1()[0][1] + dv_pert_d(rng_engine);
        rvf[1][2] = lp.get_v1()[0][2] + dv_pert_d(rng_engine);

        // We construct two problems (analytical gradient and numerical gradient)
        bench_udp_a.set_leg(rvs, mass, rvf, 1, 1, nseg, true);
        bench_udp_n.set_leg(rvs, mass, rvf, 1, 1, nseg, false);
        pagmo::problem prob_a{bench_udp_a};
        pagmo::problem prob_n{bench_udp_n};
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
    fmt::print("\nComputes the speed of a single compute_mismatch_constraints() run for a lf and hf leg.");
    perform_single_nogradient_speed_test();

    fmt::print("\nComputes the same analytical and numerical gradients and tests for speed:");
    perform_test_speed(100, 5, 10);
    perform_test_speed(100, 10, 10);
    perform_test_speed(100, 15, 10);
    perform_test_speed(100, 20, 10);
    perform_test_speed(100, 70, 10);

    // performing tests
    fmt::print("\nSolves the same optimization problems with and without analytical gradients:");
    perform_test_convergence(200, 5);
    perform_test_convergence(200, 10);
    perform_test_convergence(200, 15);

    fmt::print("\n");
}