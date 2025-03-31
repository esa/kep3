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
#include <utility>
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
        // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)                       
        std::array<std::array<double, 3>, 2> rvf, double mf, double tof, double max_thrust,
                               double isp, double mu, double cut, double tol)
        : m_rvs(rvs), m_ms(ms), m_throttles(std::move(throttles)), m_rvf(rvf), m_mf(mf), m_tof(tof), m_talphas(std::move(talphas)), m_max_thrust(max_thrust),
          m_isp(isp), m_mu(mu), m_cut(cut), m_tol(tol)
    {
        for (double m_throttle : m_throttles) {
            m_thrusts.push_back(m_throttle * m_max_thrust);
        }
    }

    void set_cut(double cut)
    {
        m_cut = cut;
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
    double m_tof = 1;
    std::vector<double> m_talphas = {0.2, 0.2, 0.2, 0.2, 0.2};
    double m_max_thrust = 1;
    double m_isp = 1;
    double m_mu = 1;
    double m_cut = 0.5;
    double m_tol = 1e-16;
    std::vector<double> m_rvms = {m_rvs[0][0], m_rvs[0][1], m_rvs[0][2], m_rvs[1][0], m_rvs[1][1], m_rvs[1][2], m_ms};
    std::vector<double> m_rvmf = {m_rvf[0][0], m_rvf[0][1], m_rvf[0][2], m_rvf[1][0], m_rvf[1][1], m_rvf[1][2], m_mf};
};

#endif