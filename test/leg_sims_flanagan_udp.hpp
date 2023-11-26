// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_TEST_LEG_SIMS_FLANAGAN_UDP_H
#define kep3_TEST_LEG_SIMS_FLANAGAN_UDP_H

#include <array>
#include <vector>

#include <pagmo/utils/gradients_and_hessians.hpp>

#include <kep3/leg/sims_flanagan.hpp>

struct sf_test_udp {
    sf_test_udp() = default;
    sf_test_udp(std::array<std::array<double, 3>, 2> rvs, double ms, std::array<std::array<double, 3>, 2> rvf,
                // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                double max_thrust = 0.05, double isp = 2000, unsigned nseg = 10u)
        : m_rvs(rvs), m_rvf(rvf), m_ms(ms), m_max_thrust(max_thrust), m_isp(isp), m_nseg(nseg)
    {
    }

    [[nodiscard]] std::vector<double> fitness(const std::vector<double> &x) const
    {
        // x = [tof (in days), mf (in kg), throttles]
        // We set the leg (avoiding the allocation for the throttles is possible but requires mutable data members.)
        kep3::leg::sims_flanagan leg(m_rvs, m_ms, std::vector<double>(m_nseg * 3, 0.), m_rvf, m_ms,
                                     x[0] * kep3::DAY2SEC, m_max_thrust, m_isp, kep3::MU_SUN);
        // We compute segments and dt
        std::size_t nseg = leg.get_throttles().size() / 3u;
        double dt = leg.get_tof() / static_cast<double>(nseg);
        // We set the throttles
        leg.set_throttles(x.begin() + 1, x.end());
        // We compute the mass schedule
        double mass = m_ms;

        for (decltype(leg.get_throttles().size()) i = 0u; i < nseg; ++i) {
            // We compute the dv
            double dv0 = leg.get_max_thrust() / mass * dt * leg.get_throttles()[3 * i];
            double dv1 = leg.get_max_thrust() / mass * dt * leg.get_throttles()[3 * i + 1];
            double dv2 = leg.get_max_thrust() / mass * dt * leg.get_throttles()[3 * i + 2];
            // Update the mass
            double norm_dv = std::sqrt(dv0 * dv0 + dv1 * dv1 + dv2 * dv2);
            mass *= std::exp(-norm_dv / leg.get_isp() / kep3::G0);
        }
        leg.set_mf(mass);

        std::vector<double> retval(7 + m_nseg, 0.);

        // Fitness
        retval[0] = -mass;
        // Equality Constraints
        auto eq_con = leg.compute_mismatch_constraints();
        retval[1] = eq_con[0] / kep3::AU;
        retval[2] = eq_con[1] / kep3::AU;
        retval[3] = eq_con[2] / kep3::AU;
        retval[4] = eq_con[3] / kep3::EARTH_VELOCITY;
        retval[5] = eq_con[4] / kep3::EARTH_VELOCITY;
        retval[6] = eq_con[5] / kep3::EARTH_VELOCITY;
        //  Inequality Constraints
        auto ineq_con = leg.compute_throttle_constraints();
        std::copy(ineq_con.begin(), ineq_con.end(), retval.begin() + 7);
        return retval;
    }

    [[nodiscard]] std::vector<double> gradient(const std::vector<double> &x) const
    {
        return pagmo::estimate_gradient([this](const std::vector<double> &x) { return this->fitness(x); }, x);
    }

    [[nodiscard]] std::pair<std::vector<double>, std::vector<double>> get_bounds() const
    {
        // x = [tof (in days), mf (in kg), throttles]
        std::vector<double> lb(m_nseg * 3 + 1, -1.);
        std::vector<double> ub(m_nseg * 3 + 1, +1.);
        lb[0] = 1.;    // days
        ub[0] = 2500.; // days
        return {lb, ub};
    }

    [[nodiscard]] static std::vector<double>::size_type get_nec()
    {
        return 6u;
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