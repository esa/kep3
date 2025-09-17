// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the term_ms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_TEST_LEG_SIMS_FLANAGAN_HF_HELPERS_H
#define kep3_TEST_LEG_SIMS_FLANAGAN_HF_HELPERS_H

#include <cstddef>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/io/xio.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <pagmo/algorithm.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/lambert_problem.hpp>
#include <kep3/leg/sims_flanagan.hpp>
#include <kep3/leg/sims_flanagan_hf.hpp>
#include <kep3/planet.hpp>
#include <kep3/ta/zero_hold_kep.hpp>
#include <kep3/udpla/vsop2013.hpp>

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

struct sf_hf_test_object {

    explicit sf_hf_test_object(double cut) : m_cut(cut) {}

    sf_hf_test_object(std::vector<double> &throttles, double cut) : m_throttles(throttles), m_cut(cut)
    {
        for (double m_throttle : m_throttles) {
            m_thrusts.push_back(m_throttle * m_max_thrust);
        }
    }

    // Retrieve mismatch constraints from manual heyoka Taylor adaptive integrator
    [[nodiscard]] std::array<double, 7> compute_manual_mc()
    {
        for (double m_throttle : m_throttles) {
            m_thrusts.push_back(m_throttle * m_max_thrust);
        }

        m_new_ta = heyoka::taylor_adaptive<double>{kep3::ta::zero_hold_kep_dyn(), m_rvms, heyoka::kw::tol = m_tol};
        *(m_new_ta.get_pars_data()) = m_mu;
        *(m_new_ta.get_pars_data() + 1) = m_isp * kep3::G0;

        // Fwd leg
        std::copy(m_thrusts.begin(), std::next(m_thrusts.begin(), 3), m_new_ta.get_pars_data() + 2);
        // Set the Taylor Integration initial conditions
        m_new_ta.set_time(0.);
        std::copy(m_rvms.begin(), m_rvms.end(), m_new_ta.get_state_data());
        // ... and integrate
        auto out = m_new_ta.propagate_until(m_tof / 2);
        std::copy(m_new_ta.get_state().begin(), m_new_ta.get_state().end(), m_fwd_final_state.begin());

        // Bck leg
        std::copy(std::next(m_thrusts.begin(), 3), std::next(m_thrusts.begin(), 6), m_new_ta.get_pars_data() + 2);
        // Set the Taylor Integration initial conditions
        m_new_ta.set_time(m_tof);
        std::copy(m_rvmf.begin(), m_rvmf.end(), m_new_ta.get_state_data());
        // ... and integrate
        auto out2 = m_new_ta.propagate_until(m_tof / 2);
        std::copy(m_new_ta.get_state().begin(), m_new_ta.get_state().end(), m_bck_final_state.begin());

        for (unsigned int i(0); i < m_mc_manual.size(); ++i) {
            m_mc_manual[i] = m_fwd_final_state[i] - m_bck_final_state[i];
        }
        return m_mc_manual;
    };

    void set_cut(double cut)
    {
        m_cut = cut;
    }

    // We assemble all constraints (equality and inequality) in a single vector from a decision vector dv
    // so we can use this in estimation of numerical gradients.
    static std::vector<double> compute_constraints(kep3::leg::sims_flanagan_hf &leg, const std::vector<double> &dv)
    {
        // First we set the leg from data
        auto nseg = leg.get_nseg();
        std::array<double, 7> rvms{};
        std::copy(dv.begin(), dv.begin() + 7, rvms.begin());
        std::vector<double> throttles(nseg * 3l);
        std::copy(dv.begin() + 7, dv.begin() + 7 + nseg * 3l, throttles.begin());
        std::array<double, 7> rvmf{};
        std::copy(dv.begin() + 7 + nseg * 3l, dv.begin() + 7 + nseg * 3l + 7, rvmf.begin());
        double time_of_flight = dv[(7 + nseg * 3 + 7 + 1) - 1];
        leg.set(rvms, throttles, rvmf, time_of_flight);
        // Then we vvaluate and return all constraints
        std::vector<double> retval(7 + nseg, 0.);
        // Equality Constraints
        auto eq_con = leg.compute_mismatch_constraints();
        retval[0] = eq_con[0];
        retval[1] = eq_con[1];
        retval[2] = eq_con[2];
        retval[3] = eq_con[3];
        retval[4] = eq_con[4];
        retval[5] = eq_con[5];
        retval[6] = eq_con[6];
        // Inequality Constraints
        auto ineq_con = leg.compute_throttle_constraints();
        std::copy(ineq_con.begin(), ineq_con.end(), retval.begin() + 7);
        return retval;
    }

    [[nodiscard]] std::vector<double> compute_numerical_gradient()
    {
        // Create SF leg.
        kep3::leg::sims_flanagan_hf sf_num(m_rvs, m_ms, m_throttles, m_rvf, m_mf, m_tof, m_max_thrust, m_isp, m_mu,
                                           m_cut, 1e-16);
        // Assemble a flattened decision vector dv
        std::vector<double> dv;
        dv.insert(dv.end(), m_rvms.begin(), m_rvms.end());
        dv.insert(dv.end(), m_throttles.begin(), m_throttles.end());
        dv.insert(dv.end(), m_rvmf.begin(), m_rvmf.end());
        dv.push_back(m_tof);

        // Calculate numerical gradient
        return pagmo::estimate_gradient_h(
            [&sf_num](const std::vector<double> &x) { return compute_constraints(sf_num, x); }, dv);
    }

    [[nodiscard]] std::vector<double> compute_analytical_gradient() const
    {
        // Initialise
        kep3::leg::sims_flanagan_hf sf_a(m_rvs, m_ms, m_throttles, m_rvf, m_mf, m_tof, m_max_thrust, m_isp, m_mu, m_cut,
                                         1e-16);
        std::array<double, 49> grad_rvm = {0};
        std::array<double, 49> grad_rvm_bck = {0};
        unsigned int nseg = static_cast<unsigned int>(m_throttles.size()) / 3;
        std::vector<double> grad_utof(static_cast<size_t>(7) * (nseg * 3u + 1u), 0.);
        std::tie(grad_rvm, grad_rvm_bck, grad_utof) = sf_a.compute_mc_grad();
        auto xgrad_rvm = xt::adapt(grad_rvm, {7u, 7u});
        auto xgrad_rvm_bck = xt::adapt(grad_rvm_bck, {7u, 7u});
        auto xgrad_utof = xt::adapt(grad_utof, {7u, nseg * 3u + 1u});

        // Cast gradients into a single vector
        std::vector<double> gradient(static_cast<size_t>(7u * (7u + static_cast<unsigned>(nseg) * 3u + 1u + 7u)), 0);
        auto xgradient = xt::adapt(gradient, {7u, 7u + static_cast<unsigned>(nseg) * 3u + 1u + 7u});
        xt::view(xgradient, xt::all(), xt::range(0u, 7u)) = xt::view(xgrad_rvm, xt::all(), xt::all()); // dmc_dxs
        xt::view(xgradient, xt::all(), xt::range(7u, 7u + nseg * 3u))
            = xt::view(xgrad_utof, xt::all(), xt::range(0, nseg * 3u)); // throttles
        xt::view(xgradient, xt::all(), xt::range(7u + nseg * 3u, 7u + nseg * 3u + 7u))
            = xt::view(xgrad_rvm_bck, xt::all(), xt::all()); // dmc_dxf
        xt::view(xgradient, xt::all(), xt::range(7u + nseg * 3u + 7u, 7u + nseg * 3u + 7u + 1u))
            = xt::view(xgrad_utof, xt::all(), xt::range(nseg * 3u, nseg * 3u + 1)); // tof

        return gradient;
    }

    // Member attributes
    std::vector<double> m_num_grad;
    heyoka::taylor_adaptive<double> m_new_ta;
    std::array<double, 7> m_fwd_final_state{};
    std::array<double, 7> m_bck_final_state{};
    std::array<double, 7> m_mc_manual{};
    std::array<std::array<double, 3>, 2> m_rvs{{{1, 0.1, -0.1}, {0.2, 1, -0.2}}};
    double m_ms = 1;
    std::vector<double> m_thrusts;
    std::vector<double> m_throttles = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::array<std::array<double, 3>, 2> m_rvf{{{1.2, -0.1, 0.1}, {-0.2, 1.023, -0.44}}};
    double m_mf = m_ms * 13 / 15;
    double m_tof = 1.234;
    double m_max_thrust = 1.6;
    double m_isp = 1.234;
    double m_mu = 0.9234;
    double m_cut = 0.5;
    double m_tol = 1e-16;
    std::vector<double> m_rvms = {m_rvs[0][0], m_rvs[0][1], m_rvs[0][2], m_rvs[1][0], m_rvs[1][1], m_rvs[1][2], m_ms};
    std::vector<double> m_rvmf = {m_rvf[0][0], m_rvf[0][1], m_rvf[0][2], m_rvf[1][0], m_rvf[1][1], m_rvf[1][2], m_mf};
};

#endif