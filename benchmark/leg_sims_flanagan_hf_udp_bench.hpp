// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_TEST_LEG_SIMS_FLANAGAN_HF_UDP_BENCH_H
#define kep3_TEST_LEG_SIMS_FLANAGAN_HF_UDP_BENCH_H

#include <array>
#include <vector>

#include <xtensor/xadapt.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>

#include <pagmo/utils/gradients_and_hessians.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/leg/sims_flanagan_hf.hpp>

struct sf_hf_bench_udp {
    sf_hf_bench_udp() = default;
    sf_hf_bench_udp(std::array<std::array<double, 3>, 2> rvs, double ms, std::array<std::array<double, 3>, 2> rvf,
                    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                    double max_thrust, double isp, unsigned nseg, bool analytical)
        : m_rvs(rvs), m_rvf(rvf), m_ms(ms), m_max_thrust(max_thrust), m_isp(isp), m_nseg(nseg),
          m_analytical(analytical),
          m_leg(kep3::leg::sims_flanagan_hf(m_rvs, m_ms, std::vector<double>(m_nseg * 3, 0.), m_rvf, 0.0, 0.0,
                                            m_max_thrust, m_isp, 1.0, 0.5, 1e-16))
    {
    }
    [[nodiscard]] void create_leg(std::array<std::array<double, 3>, 2> rvs, double ms,
                                  std::array<std::array<double, 3>, 2> rvf,
                                  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                                  double max_thrust, double isp, unsigned nseg, bool analytical)
    {
        m_rvs = rvs;
        m_rvf = rvf;
        m_ms = ms;
        m_max_thrust = max_thrust;
        m_isp = isp;
        m_nseg = nseg;
        m_analytical = analytical;
        m_leg = kep3::leg::sims_flanagan_hf(m_rvs, m_ms, std::vector<double>(m_nseg * 3, 0.), m_rvf, 0.0, 0.0,
                                            m_max_thrust, m_isp, 1.0, 0.5, 1e-16);
    }

    [[nodiscard]] void set_leg(std::array<std::array<double, 3>, 2> rvs, double ms,
                               std::array<std::array<double, 3>, 2> rvf,
                               // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                               double max_thrust, double isp, unsigned nseg, bool analytical)
    {
        m_rvs = rvs;
        m_rvf = rvf;
        m_ms = ms;
        m_max_thrust = max_thrust;
        m_isp = isp;
        m_nseg = nseg;
        m_analytical = analytical;
        m_leg.set(m_rvs, m_ms, std::vector<double>(m_nseg * 3, 0.), m_rvf, 0.0, 0.0, m_max_thrust, m_isp, 1.0, 0.5,
                  1e-16);
    }

    [[nodiscard]] std::vector<double> fitness(const std::vector<double> &x) const
    {
        // x = [throttles, tof (in days), mf (in kg)]
        // We set the leg (avoiding the allocation for the throttles is possible but requires mutable data members.)
        double tof = x[m_nseg * 3];    // in s
        double mf = x[m_nseg * 3 + 1]; // in kg
        m_leg.set_tof(tof);
        m_leg.set_mf(mf);

        // We set the throttles
        m_leg.set_throttles(x.begin(), x.end() - 2);

        std::vector<double> retval(1 + 7 + m_nseg, 0.);
        // Fitness
        retval[0] = -mf;
        // Equality Constraints
        auto eq_con = m_leg.compute_mismatch_constraints();
        retval[1] = eq_con[0];
        retval[2] = eq_con[1];
        retval[3] = eq_con[2];
        retval[4] = eq_con[3];
        retval[5] = eq_con[4];
        retval[6] = eq_con[5];
        retval[7] = eq_con[6];
        //  Inequality Constraints
        auto ineq_con = m_leg.compute_throttle_constraints();
        std::copy(ineq_con.begin(), ineq_con.end(), retval.begin() + 8);
        return retval;
    }

    [[nodiscard]] std::vector<double> gradient(const std::vector<double> &x) const
    {
        if (m_analytical) {
            return _gradient_analytical(x);
        } else {
            return _gradient_numerical(x);
        }
    }

    [[nodiscard]] std::vector<double> _gradient_numerical(const std::vector<double> &x) const
    {
        auto num_grad = pagmo::estimate_gradient([this](const std::vector<double> &x) { return this->fitness(x); }, x);
        return num_grad;
    }

    [[nodiscard]] std::vector<double> _gradient_analytical(const std::vector<double> &x) const
    {
        // x = [throttles, tof (in days), mf (in kg)]
        // We set the leg (avoiding the allocation for the throttles is possible but requires mutable data members.)
        double tof = x[m_nseg * 3];    // in s
        double mf = x[m_nseg * 3 + 1]; // in kg
        m_leg.set_tof(tof);
        m_leg.set_mf(mf);
        // We set the throttles
        m_leg.set_throttles(x.begin(), x.end() - 2);

        // We compute the gradients
        std::array<double, 49> grad_rvm = {0};
        std::array<double, 49> grad_rvm_bck = {0};
        std::vector<double> grad_final(static_cast<size_t>(7) * (m_nseg * 3u + 1u), 0.);
        std::tie(grad_rvm, grad_rvm_bck, grad_final) = m_leg.compute_mc_grad();
        auto xgrad_rvm = xt::adapt(grad_rvm, {7u, 7u});
        auto xgrad_rvm_bck = xt::adapt(grad_rvm_bck, {7u, 7u});
        auto xgrad_final = xt::adapt(grad_final, {7u, static_cast<unsigned int>(m_nseg) * 3u + 1u});

        std::vector<double> grad_tc = m_leg.compute_tc_grad();
        auto xt_grad_tc = xt::adapt(grad_tc, {m_nseg, 3u * m_nseg});

        // Initialise gradient
        std::vector<double> gradient((1u + 7u + m_nseg) * (m_nseg * 3u + 2u), 0);
        // Create the various xtensor objects adapting the std containers
        auto xgradient
            = xt::adapt(gradient, {1u + 7u + static_cast<unsigned>(m_nseg), static_cast<unsigned>(m_nseg) * 3u + 2u});

        xgradient(0, m_nseg * 3 + 1) = -1.; // fitness gradient - obj fun
        xt::view(xgradient, xt::range(1u, 4u), xt::range(0, m_nseg * 3u + 1u))
            = xt::view(xgrad_final, xt::range(0u, 3u), xt::all()); // dmc/du
        xt::view(xgradient, xt::range(4u, 7u), xt::range(0, m_nseg * 3u + 1u))
            = xt::view(xgrad_final, xt::range(3u, 6u), xt::all()); // dmc/du
        xt::view(xgradient, xt::range(7u, 8u), xt::range(0, m_nseg * 3u + 1u))
            = xt::view(xgrad_final, xt::range(6u, 7u), xt::all()); // dmc/du

        xt::view(xgradient, xt::range(1u, 4u), xt::range(m_nseg * 3u + 1u, m_nseg * 3u + 2u))
            = xt::view(xgrad_rvm_bck, xt::range(0u, 3u), xt::range(6u, 7u)); // dmc/dm_f
        xt::view(xgradient, xt::range(4u, 7u), xt::range(m_nseg * 3u + 1u, m_nseg * 3u + 2u))
            = xt::view(xgrad_rvm_bck, xt::range(3u, 6u), xt::range(6u, 7u)); // dmc/dm_f
        xt::view(xgradient, xt::range(7u, 8u), xt::range(m_nseg * 3u + 1u, m_nseg * 3u + 2u))
            = xt::view(xgrad_rvm_bck, xt::range(6u, 7u), xt::range(6u, 7u)); // dmc/dm_f
        xt::view(xgradient, xt::range(8u, 8u + m_nseg), xt::range(0, m_nseg * 3u))
            = xt::view(xt_grad_tc, xt::all(), xt::all()); // throttle constraints

        xt::view(xgradient, xt::all(), xt::range(m_nseg * 3u, m_nseg * 3u + 1u));

        return gradient;
    }

    [[nodiscard]] std::pair<std::vector<double>, std::vector<double>> get_bounds() const
    {
        // x = [throttles, tof (in days), mf (in kg)]
        std::vector<double> lb(m_nseg * 3 + 2, -1.);
        std::vector<double> ub(m_nseg * 3 + 2, +1.);
        lb[m_nseg * 3] = kep3::pi / 12;     // days
        ub[m_nseg * 3] = 2 * kep3::pi;     // days
        lb[m_nseg * 3 + 1] = 0.5; // kg
        ub[m_nseg * 3 + 1] = 1;   // kg
        return {lb, ub};
    }

    [[nodiscard]] static std::vector<double>::size_type get_nec()
    {
        return 7u;
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
    bool m_analytical{};
    // m_leg needs to be mutable because the heyoka integrator needs to be modifiable
    mutable kep3::leg::sims_flanagan_hf m_leg{};
};

#endif