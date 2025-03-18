// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the term_ms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_TEST_LEG_SIMS_FLANAGAN_HF_ALPHA_HELPERS_H
#define kep3_TEST_LEG_SIMS_FLANAGAN_HF_ALPHA_HELPERS_H

#include <cstddef>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xio.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <pagmo/algorithm.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/lambert_problem.hpp>
#include <kep3/leg/sims_flanagan.hpp>
#include <kep3/leg/sims_flanagan_hf_alpha.hpp>
#include <kep3/planet.hpp>
#include <kep3/ta/stark.hpp>
#include <kep3/udpla/vsop2013.hpp>

#include <pagmo/utils/gradients_and_hessians.hpp>
#include <xtensor/xview.hpp>

#include <heyoka/config.hpp>
#include <heyoka/expression.hpp>
#include <heyoka/math/pow.hpp>
#include <heyoka/math/relational.hpp>
#include <heyoka/math/select.hpp>
#include <heyoka/math/sqrt.hpp>
#include <heyoka/math/sum.hpp>
#include <heyoka/taylor.hpp>

struct sf_hf_test_alpha_object {

    // Default constructor
    sf_hf_test_alpha_object() = default;

    explicit sf_hf_test_alpha_object(std::vector<double> &throttles) : m_throttles(throttles)
    {
        for (double m_throttle : m_throttles) {
            m_thrusts.push_back(m_throttle * m_max_thrust);
        }
    }

    explicit sf_hf_test_alpha_object(double cut) : m_cut(cut) {}

    sf_hf_test_alpha_object(std::vector<double> &throttles, double cut) : m_throttles(throttles), m_cut(cut)
    {
        for (double m_throttle : m_throttles) {
            m_thrusts.push_back(m_throttle * m_max_thrust);
        }
    }

    explicit sf_hf_test_alpha_object(std::array<std::array<double, 3>, 2> rvs, double ms, std::vector<double> throttles, std::vector<double> talphas,
                               std::array<std::array<double, 3>, 2> rvf, double mf, double tof, double max_thrust,
                               double isp, double mu, double cut, double tol)
        : m_rvs(rvs), m_ms(ms), m_throttles(throttles), m_talphas(talphas), m_rvf(rvf), m_mf(mf), m_tof(tof), m_max_thrust(max_thrust),
          m_isp(isp), m_mu(mu), m_cut(cut), m_tol(tol)
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

        m_new_ta = heyoka::taylor_adaptive<double>{kep3::ta::stark_dyn(), m_rvms, heyoka::kw::tol = m_tol};
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

    [[nodiscard]] std::vector<double> compute_numerical_gradient()
    {
        // Create SF leg.
        kep3::leg::sims_flanagan_hf_alpha sf_num(m_rvs, m_ms, m_throttles, m_talphas, m_rvf, m_mf, m_tof, m_max_thrust, m_isp, m_mu,
                                           m_cut, 1e-16);
        // Create chromosome
        std::vector<double> rvms_vec = std::vector<double>(m_rvms.begin(), m_rvms.end());
        std::vector<double> rvmf_vec = std::vector<double>(m_rvmf.begin(), m_rvmf.end());
        std::vector<double> chromosome;
        chromosome.insert(chromosome.end(), rvms_vec.begin(), rvms_vec.end());
        chromosome.insert(chromosome.end(), m_throttles.begin(), m_throttles.end());
        chromosome.insert(chromosome.end(), rvmf_vec.begin(), rvmf_vec.end());
        chromosome.push_back(m_tof);

        // Calculate numerical gradient
        return pagmo::estimate_gradient_h(
            [&sf_num](const std::vector<double> &x) { return sf_num.set_and_compute_constraints(x); }, chromosome);
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
    std::vector<double> m_talphas = {0, 0, 0, 0, 0};
    std::array<std::array<double, 3>, 2> m_rvf{{{1.2, -0.1, 0.1}, {-0.2, 1.023, -0.44}}};
    double m_mf = m_ms * 13 / 15;
    double m_tof = 1;
    double m_max_thrust = 1;
    double m_isp = 1;
    double m_mu = 1;
    double m_cut = 0.5;
    double m_tol = 1e-16;
    std::vector<double> m_rvms = {m_rvs[0][0], m_rvs[0][1], m_rvs[0][2], m_rvs[1][0], m_rvs[1][1], m_rvs[1][2], m_ms};
    std::vector<double> m_rvmf = {m_rvf[0][0], m_rvf[0][1], m_rvf[0][2], m_rvf[1][0], m_rvf[1][1], m_rvf[1][2], m_mf};
};

#endif