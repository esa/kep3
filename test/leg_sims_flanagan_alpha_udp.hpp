// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_TEST_LEG_SIMS_FLANAGAN_ALPHA_UDP_H
#define kep3_TEST_LEG_SIMS_FLANAGAN_ALPHA_UDP_H

#include <array>
#include <vector>

#include <xtensor/xadapt.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>

#include <pagmo/utils/gradients_and_hessians.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/leg/sims_flanagan_alpha.hpp>

struct sf_test_udp {
    sf_test_udp() = default;
    sf_test_udp(std::array<std::array<double, 3>, 2> rvs, double ms, std::array<std::array<double, 3>, 2> rvf,
                // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                double max_thrust, double isp, unsigned nseg)
        : m_rvs(rvs), m_rvf(rvf), m_ms(ms), m_max_thrust(max_thrust), m_isp(isp), m_nseg(nseg)
    {
    }

    [[nodiscard]] std::vector<double> fitness(const std::vector<double> &x) const
    {
        // x = [throttles, tof (in days), mf (in kg)]
        // We set the leg (avoiding the allocation for the throttles is possible but requires mutable data members.)
        double tof = x[m_nseg * 3] * kep3::DAY2SEC; // in s
        double mf = x[m_nseg * 3 + 1];              // in kg
        kep3::leg::sims_flanagan_alpha leg(m_rvs, m_ms, std::vector<double>(m_nseg * 3, 0.), std::vector<double>(m_nseg, tof/m_nseg), m_rvf, mf, tof, m_max_thrust,
                                     m_isp, kep3::MU_SUN);

        // We set the throttles
        leg.set_throttles(x.begin(), x.end() - 2);

        std::vector<double> retval(1 + 7 + m_nseg, 0.);
        // Fitness
        retval[0] = -mf;
        // Equality Constraints
        auto eq_con = leg.compute_mismatch_constraints();
        retval[1] = eq_con[0] / kep3::AU;
        retval[2] = eq_con[1] / kep3::AU;
        retval[3] = eq_con[2] / kep3::AU;
        retval[4] = eq_con[3] / kep3::EARTH_VELOCITY;
        retval[5] = eq_con[4] / kep3::EARTH_VELOCITY;
        retval[6] = eq_con[5] / kep3::EARTH_VELOCITY;
        retval[7] = eq_con[6] / 1e8; //
        //  Inequality Constraints
        auto ineq_con = leg.compute_throttle_constraints();
        std::copy(ineq_con.begin(), ineq_con.end(), retval.begin() + 8);
        return retval;
    }

    [[nodiscard]] std::vector<double> gradient_numerical(const std::vector<double> &x) const
    {
        return pagmo::estimate_gradient_h([this](const std::vector<double> &x) { return this->fitness(x); }, x);
    }

    [[nodiscard]] std::vector<double> gradient(const std::vector<double> &x) const
    {
        // Currently numerical, in the future we can add analytical
        return gradient_numerical(x);
    }

    [[nodiscard]] std::pair<std::vector<double>, std::vector<double>> get_bounds() const
    {
        // x = [throttles, tof (in days), mf (in kg)]
        std::vector<double> lb(m_nseg * 3 + 2, -1.);
        std::vector<double> ub(m_nseg * 3 + 2, +1.);
        lb[m_nseg * 3] = 1.;            // days
        ub[m_nseg * 3] = 2500.;         // days
        lb[m_nseg * 3 + 1] = m_ms / 2.; // kg
        ub[m_nseg * 3 + 1] = m_ms;      // kg
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
};

#endif