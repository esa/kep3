// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_TEST_LEG_SIMS_FLANAGAN_HF_ALPHA_UDP_H
#define kep3_TEST_LEG_SIMS_FLANAGAN_HF_ALPHA_UDP_H

#include <array>
#include <vector>

#include <xtensor/xadapt.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>

#include <pagmo/utils/gradients_and_hessians.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/leg/sims_flanagan_hf_alpha.hpp>

#include <kep3/core_astro/encodings.hpp>

struct sf_hf_alpha_test_udp {
    sf_hf_alpha_test_udp() = default;
    sf_hf_alpha_test_udp(std::array<std::array<double, 3>, 2> rvs, double ms, std::array<std::array<double, 3>, 2> rvf,
                // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                double max_thrust, double isp, unsigned nseg)
        : m_rvs(rvs), m_rvf(rvf), m_ms(ms), m_max_thrust(max_thrust), m_isp(isp), m_nseg(nseg)
    {
    }

    [[nodiscard]] std::vector<double> fitness(const std::vector<double> &x) const
    {
        // x = [throttles, alphas, tof (in days), mf (in kg)]
        // We set the leg (avoiding the allocation for the throttles is possible but requires mutable data members.)
        double tof = x[m_nseg * 4] * kep3::DAY2SEC; // in s
        double mf = x[m_nseg * 4 + 1];              // in kg
        kep3::leg::sims_flanagan_hf_alpha leg(m_rvs, m_ms, std::vector<double>(m_nseg * 3, 0.), std::vector<double>(m_nseg, tof/m_nseg), m_rvf, mf, tof, m_max_thrust,
                                     m_isp, kep3::MU_SUN);

        // We set the throttles
        int nseg = static_cast<int>(m_nseg);  // Convert safely
        leg.set_throttles(x.begin(), x.end() - 2 - nseg);

        // Transform alphas using alpha2direct
        const std::vector<double> alphas(x.begin() + 3 * nseg, x.end() - 2);
        std::vector<double> transformed_alphas = kep3::alpha2direct(alphas, tof);
        leg.set_talphas(transformed_alphas);

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
        return pagmo::estimate_gradient_h([this](const std::vector<double> &x) { return this->fitness(x); }, x, 1e-4);
    }

    [[nodiscard]] std::vector<double> gradient(const std::vector<double> &x) const
    {
        // Currently numerical, in the future we can add analytical
        return gradient_numerical(x);
    }

    [[nodiscard]] std::pair<std::vector<double>, std::vector<double>> get_bounds() const
    {
        // x = [throttles, alphas, tof (in days), mf (in kg)]
        int nseg = static_cast<int>(m_nseg);  // Convert safely
        std::vector<double> lb(nseg * 4 + 2, -1.);
        std::vector<double> ub(nseg * 4 + 2, +1.);
        // Set specific values for the range [m_nseg*3, m_nseg*4)
        std::fill(lb.begin() + nseg * 3, lb.begin() + nseg * 4, 0.7);
        std::fill(ub.begin() + nseg * 3, ub.begin() + nseg * 4, 0.9);
        lb[nseg * 4] = 1.;            // days
        ub[nseg * 4] = 2500.;         // days
        lb[nseg * 4 + 1] = m_ms / 2.; // kg
        ub[nseg * 4 + 1] = m_ms;      // kg
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