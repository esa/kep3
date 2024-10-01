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
#include <kep3/leg/sims_flanagan_hf.hpp>
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

struct sf_hf_test_object {

    // Default constructor
    sf_hf_test_object() = default;

    sf_hf_test_object(std::vector<double> throttles) : m_throttles(throttles)
    {
        for (unsigned int i(0); i < m_throttles.size(); ++i) {
            m_thrusts.push_back(m_throttles[i] * m_max_thrust);
        }
    }

    sf_hf_test_object(double cut) : m_cut(cut) {}

    sf_hf_test_object(std::vector<double> throttles, double cut) : m_cut(cut), m_throttles(throttles)
    {
        for (unsigned int i(0); i < m_throttles.size(); ++i) {
            m_thrusts.push_back(m_throttles[i] * m_max_thrust);
        }
    }

    // Retrieve mismatch constraints from manual heyoka Taylor adaptive integrator
    [[nodiscard]] std::array<double, 7> compute_manual_mc()
    {
        for (unsigned int i(0); i < m_throttles.size(); ++i) {
            m_thrusts.push_back(m_throttles[i] * m_max_thrust);
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

    [[nodiscard]] void set_cut(double cut)
    {
        m_cut = cut;
    }

    std::vector<double> compute_numerical_gradient()
    {
        // Create SF leg.
        kep3::leg::sims_flanagan_hf sf_num(m_rvs, m_ms, m_throttles, m_rvf, m_mf, m_tof, m_max_thrust, m_isp, m_mu,
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

    std::tuple<xt::xarray<double>, xt::xarray<double>, xt::xarray<double>, xt::xarray<double>, xt::xarray<double>,
               xt::xarray<double>, xt::xarray<double>, xt::xarray<double>>
    process_mc_numerical_gradient(std::vector<double> num_grad)
    {
        m_num_grad = num_grad;

        std::array<std::array<double, 7>, 7> num_dmc_dxs = {{{0}}};
        std::array<std::array<double, 7>, 7> num_dmc_dxf = {{{0}}};
        std::array<std::array<double, 15>, 7> num_dmc_du = {{{0}}};
        std::array<double, 7> num_dmc_dtof = {{0}};

        // Loop over first 7 constraints (the states) and fill in the respective matrices
        for (unsigned int i(0); i < num_dmc_dxs.size(); ++i) {
            // dmc_dxs
            std::copy(std::next(m_num_grad.begin(), 30 * i), std::next(m_num_grad.begin(), 7 + 30 * i),
                      std::next(num_dmc_dxs.begin(), i)->begin());
            // dmc_du
            std::copy(std::next(m_num_grad.begin(), 7 + 30 * i), std::next(m_num_grad.begin(), 22 + 30 * i),
                      std::next(num_dmc_du.begin(), i)->begin());
            // dmc_dxf
            std::copy(std::next(m_num_grad.begin(), 22 + 30 * i), std::next(m_num_grad.begin(), 29 + 30 * i),
                      std::next(num_dmc_dxf.begin(), i)->begin());
            // dmc_dtof
            num_dmc_dtof[i] = m_num_grad[29 + 30 * i];
        }

        xt::xarray<double> xt_num_dmc_dxs = xt::adapt(reinterpret_cast<double *>(num_dmc_dxs.data()), {7, 7});
        xt::xarray<double> xt_num_dmc_dxf = xt::adapt(reinterpret_cast<double *>(num_dmc_dxf.data()), {7, 7})
                                            * -1; // Multiple by -1 because mass correlation is -1.
        xt::xarray<double> xt_num_dmc_du = xt::adapt(reinterpret_cast<double *>(num_dmc_du.data()), {7, 15});
        xt::xarray<double> xt_num_dmc_dtof = xt::adapt(reinterpret_cast<double *>(num_dmc_dtof.data()), {7, 1});
        auto xt_num_dmc_du0 = xt::view(xt_num_dmc_du, xt::all(), xt::range(0, 3));
        auto xt_num_dmc_du1 = xt::view(xt_num_dmc_du, xt::all(), xt::range(3, 6));
        auto xt_num_dmc_du2 = xt::view(xt_num_dmc_du, xt::all(), xt::range(6, 9));
        auto xt_num_dmc_du3 = xt::view(xt_num_dmc_du, xt::all(), xt::range(9, 12))
                              * -1; // Multiple by -1 because mass correlation is -1.
        auto xt_num_dmc_du4 = xt::view(xt_num_dmc_du, xt::all(), xt::range(12, 15))
                              * -1; // Multiple by -1 because mass correlation is -1.

        return std::make_tuple(xt_num_dmc_dxs, xt_num_dmc_dxf, xt_num_dmc_du0, xt_num_dmc_du1, xt_num_dmc_du2,
                               xt_num_dmc_du3, xt_num_dmc_du4, xt_num_dmc_dtof);
    }

    std::tuple<xt::xarray<double>, xt::xarray<double>, xt::xarray<double>, xt::xarray<double>, xt::xarray<double>,
               xt::xarray<double>, xt::xarray<double>, xt::xarray<double>>
    compute_analytical_gradient()
    {
        unsigned int nseg = static_cast<unsigned int>(m_throttles.size()) / 3;
        auto nseg_fwd = nseg * 0.6;
        auto nseg_bck = nseg - nseg_fwd;
        // Initialise
        std::array<std::array<double, 7u>, 5u> xf_per_seg;
        std::array<std::array<double, 7u>, 5u> x0_per_seg;
        std::array<std::array<double, 49u>, 5u> dxdx_per_seg;
        std::array<std::array<double, 21u>, 5u> dxdu_per_seg;
        std::array<std::array<double, 7u>, 5u> dxdtof_per_seg;
        kep3::leg::sims_flanagan_hf sf_a(m_rvs, m_ms, m_throttles, m_rvf, m_mf, m_tof, m_max_thrust, m_isp, m_mu, m_cut,
                                         1e-16);
        std::tie(xf_per_seg, dxdx_per_seg, dxdu_per_seg) = sf_a.compute_mc_grad();
        // Initialize initial state matrix
        x0_per_seg[0] = sf_a.get_rvms();
        for (unsigned int i(0); i < nseg_fwd - 1; ++i) {
            x0_per_seg[i + 1] = xf_per_seg[i];
        }
        x0_per_seg[nseg - 1] = sf_a.get_rvmf();
        for (unsigned int i(0); i < nseg_bck - 1; ++i) {
            x0_per_seg[(nseg - 1) - (i + 1)] = xf_per_seg[(nseg - 1) - i];
        }

        for (unsigned int i(0); i < dxdtof_per_seg.size(); ++i) {
            std::array<double, 3> current_throttles
                = {m_throttles[i * 3], m_throttles[i * 3 + 1], m_throttles[i * 3 + 2]};
            dxdtof_per_seg[i] = sf_a.get_state_derivative(x0_per_seg[i], current_throttles);
        }

        xt::xarray<double> xt_dxdx_per_seg = xt::zeros<double>({5, 49});
        for (size_t col = 0; col < 49; ++col) {    // Iterate over columns
            for (size_t row = 0; row < 5; ++row) { // Iterate over rows
                xt_dxdx_per_seg(row, col) = dxdx_per_seg[row][col];
            }
        }

        // Create matrices from final states to previous states
        // Fwd leg
        auto M0 = xt::reshape_view(xt::view(xt_dxdx_per_seg, 0, xt::all()), {7, 7});
        auto M1 = xt::reshape_view(xt::view(xt_dxdx_per_seg, 1, xt::all()), {7, 7});
        auto M2 = xt::reshape_view(xt::view(xt_dxdx_per_seg, 2, xt::all()), {7, 7});
        // Bck leg
        auto M3 = xt::reshape_view(xt::view(xt_dxdx_per_seg, 3, xt::all()), {7, 7});
        auto M4 = xt::reshape_view(xt::view(xt_dxdx_per_seg, 4, xt::all()), {7, 7});

        // Create matrices from final states to throttles
        xt::xarray<double> xt_dxdu_per_seg = xt::adapt(reinterpret_cast<double *>(dxdu_per_seg.data()), {5, 21});
        // Fwd leg
        auto U0 = xt::reshape_view(xt::view(xt_dxdu_per_seg, 0, xt::all()), {7, 3});
        auto U1 = xt::reshape_view(xt::view(xt_dxdu_per_seg, 1, xt::all()), {7, 3});
        auto U2 = xt::reshape_view(xt::view(xt_dxdu_per_seg, 2, xt::all()), {7, 3});
        // Bck leg
        auto U3 = xt::reshape_view(xt::view(xt_dxdu_per_seg, 3, xt::all()), {7, 3});
        auto U4 = xt::reshape_view(xt::view(xt_dxdu_per_seg, 4, xt::all()), {7, 3});

        // Create matrices from final states to throttles
        xt::xarray<double> xt_dxdtof_per_seg = xt::adapt(reinterpret_cast<double *>(dxdtof_per_seg.data()), {5, 7});
        // Fwd leg
        auto f0 = xt::reshape_view(xt::view(xt_dxdtof_per_seg, 0, xt::all()), {7, 1});
        auto f1 = xt::reshape_view(xt::view(xt_dxdtof_per_seg, 1, xt::all()), {7, 1});
        auto f2 = xt::reshape_view(xt::view(xt_dxdtof_per_seg, 2, xt::all()), {7, 1});
        // Bck leg
        auto f5 = xt::reshape_view(xt::view(xt_dxdtof_per_seg, 3, xt::all()),
                                   {7, 1}); // f5 is because we don't care about f3 and f4 and the mc
        auto f6 = xt::reshape_view(xt::view(xt_dxdtof_per_seg, 4, xt::all()), {7, 1});

        // Initial and final displacements
        auto xt_a_dmc_dxs = xt::linalg::dot(xt::linalg::dot(M2, M1), M0);
        auto xt_a_dmc_dxf = xt::linalg::dot(M3, M4);
        // Throttle derivatives
        auto xt_a_dmc_du0 = xt::linalg::dot(xt::linalg::dot(M2, M1), U0);
        auto xt_a_dmc_du1 = xt::linalg::dot(M2, U1);
        auto xt_a_dmc_du2 = U2;
        auto xt_a_dmc_du3 = U3;
        auto xt_a_dmc_du4 = xt::linalg::dot(M3, U4);
        // ToF derivatives
        auto xt_a_dmc_dtof0 = xt::linalg::dot(xt::linalg::dot(xt::linalg::dot(M2, M1), M0), f0);
        auto xt_a_dmc_dtof1 = xt::linalg::dot(xt::linalg::dot(M2, M1), f1);
        auto xt_a_dmc_dtof2 = xt::linalg::dot(M2, f2);
        auto xt_a_dmc_dtof3 = xt::linalg::dot(M3, f5);
        auto xt_a_dmc_dtof4 = xt::linalg::dot(xt::linalg::dot(M3, M4), f6);
        auto xt_a_dmc_dtof
            = (xt_a_dmc_dtof0 + xt_a_dmc_dtof1 + xt_a_dmc_dtof2 + xt_a_dmc_dtof3 + xt_a_dmc_dtof4) / nseg;

        return std::make_tuple(xt_a_dmc_dxs, xt_a_dmc_dxf, xt_a_dmc_du0, xt_a_dmc_du1, xt_a_dmc_du2, xt_a_dmc_du3,
                               xt_a_dmc_du4, xt_a_dmc_dtof);
    }

    std::vector<double> m_num_grad;

    heyoka::taylor_adaptive<double> m_new_ta;
    std::array<double, 7> m_fwd_final_state;
    std::array<double, 7> m_bck_final_state;
    std::array<double, 7> m_mc_manual;

    std::array<std::array<double, 3>, 2> m_rvs{{{1, 0.1, -0.1}, {0.2, 1, -0.2}}};

    std::array<std::array<double, 3>, 2> m_rvf{{{1.2, -0.1, 0.1}, {-0.2, 1.023, -0.44}}};

    double m_ms = 1;
    double m_mf = m_ms * 13 / 15;
    double m_isp = 1;
    double m_max_thrust = 1;
    double m_cut = 0.5;
    double m_mu = 1;

    double m_tof = 1;

    std::vector<double> m_throttles = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> m_thrusts;

    double m_tol = 1e-16;
    const std::vector<double> m_rvms
        = {m_rvs[0][0], m_rvs[0][1], m_rvs[0][2], m_rvs[1][0], m_rvs[1][1], m_rvs[1][2], m_ms};
    const std::vector<double> m_rvmf
        = {m_rvf[0][0], m_rvf[0][1], m_rvf[0][2], m_rvf[1][0], m_rvf[1][1], m_rvf[1][2], m_mf};
};

#endif