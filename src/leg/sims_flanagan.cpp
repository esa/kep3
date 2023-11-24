// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <iomanip>
#include <vector>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/epoch.hpp>
#include <kep3/leg/sims_flanagan.hpp>

namespace kep3::leg
{
sims_flanagan::sims_flanagan(const std::array<double, 7> &xs, std::vector<double> throttles,
                             // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                             const std::array<double, 7> &xf, double tof, double max_thrust, double isp, double mu,
                             double cut)
    : m_xs(xs), m_throttles(std::move(throttles)), m_xf(xf), m_tof(tof), m_max_thrust(max_thrust), m_isp(isp), m_mu(mu),
      m_cut(cut)
{
    // Sanity Checks
    if (m_tof < 0.) {
        throw std::domain_error("The time of flight of a sims_flanagan leg needs to be larger or equal to zero.");
    }
    if (m_throttles.size() % 3 != 0u) {
        throw std::domain_error("The throttles of a sims_flanagan leg are detected to be not a multiple of 3 in size "
                                "[u0x, u0y, u0z, .....].");
    }
    if (m_max_thrust < 0.) {
        throw std::domain_error(
            "The maximum allowed thrust of a sims_flanagan leg is detected to be smaller than zero.");
    }
    if (m_isp < 0.) {
        throw std::domain_error("The specific impulse of a sims_flanagan leg is detected to be smaller than zero.");
    }
    if (m_mu < 0.) {
        throw std::domain_error(
            "The gravitational parameter of a sims_flanagan leg is detected to be smaller than zero.");
    }
    if (m_cut < 1 || m_cut > 1.) {
        throw std::domain_error("The parameter cut of a sims_flanagan leg must be in [0, 1].");
    }
};

// Setters
void sims_flanagan::set_tof(double tof)
{
    m_tof = tof;
};
void sims_flanagan::set_xs(std::array<double, 7> xs)
{
    m_xs = xs;
};
void sims_flanagan::set_throttles(std::vector<double> throttles)
{
    m_throttles = std::move(throttles);
};
void sims_flanagan::set_xf(std::array<double, 7> xf)
{
    m_xf = xf;
};
void sims_flanagan::set_max_thrust(double max_thrust)
{
    m_max_thrust = max_thrust;
}
void sims_flanagan::set_isp(double isp)
{
    m_isp = isp;
}
void sims_flanagan::set_mu(double mu)
{
    m_mu = mu;
}
void sims_flanagan::set_cut(double cut)
{
    m_cut = cut;
}

// Getters
double sims_flanagan::get_tof() const
{
    return m_tof;
};
const std::array<double, 7> &sims_flanagan::get_xs() const
{
    return m_xs;
};
const std::vector<double> &sims_flanagan::get_throttles() const
{
    return m_throttles;
};
const std::array<double, 7> &sims_flanagan::get_xf() const
{
    return m_xf;
};
double sims_flanagan::get_max_thrust() const
{
    return m_max_thrust;
};
double sims_flanagan::get_isp() const
{
    return m_isp;
};
double sims_flanagan::get_mu() const
{
    return m_mu;
};
double sims_flanagan::get_cut() const
{
    return m_cut;
};

// The core routines
std::array<double, 7> sims_flanagan::compute_mismatch_constraints() const
{
    // We start defining the number of forward and backward segments.
    auto n_seg = m_throttles.size() / 3u;
    auto n_seg_fwd = static_cast<unsigned>(static_cast<double>(n_seg) * m_cut);
    auto n_seg_bck = n_seg - n_seg_fwd;

    // We introduce some convenience variables
    double max_thrust = get_max_thrust();
    double isp = get_isp();
    double mu = get_mu();
    auto xs = get_xs();
    auto xf = get_xf();
    std::array<double, 3> dv{};

    // Forward pass
    // Initial state
    std::array<std::array<double, 3>, 2> rv_fwd{{{xs[0], xs[1], xs[2]}, {xs[3], xs[4], xs[5]}}};
    double mass_fwd = xs[6];
    double dt = m_tof / static_cast<double>(n_seg);
    // We propagate for a first dt/2
    rv_fwd = propagate_lagrangian(rv_fwd, dt / 2, mu, false).first;
    // We now loop through the forward segments and 1) add a dv + 2) propagate for dt (except on the last segment, where
    // we propagate for dt/2).
    for (decltype(m_throttles.size()) i = 0u; i < n_seg_fwd; ++i) {
        // We compute the the dv
        dv[0] = max_thrust / mass_fwd * dt * m_throttles[3 * i];
        dv[1] = max_thrust / mass_fwd * dt * m_throttles[3 * i + 1];
        dv[2] = max_thrust / mass_fwd * dt * m_throttles[3 * i + 2];
        // Add it to the current spacecraft velocity
        rv_fwd[1][0] += dv[0];
        rv_fwd[1][1] += dv[1];
        rv_fwd[1][2] += dv[2];
        // Update the mass accordingly
        double norm_dv = std::sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
        mass_fwd *= std::exp(-norm_dv / isp / kep3::G0);
        // Perform the propagation
        double prop_duration = (i == n_seg_fwd - 1) ? dt / 2 : dt;
        rv_fwd = propagate_lagrangian(rv_fwd, prop_duration, mu, false).first;
    }

    // Backward pass
    // Final state
    std::array<std::array<double, 3>, 2> rv_bck{{{xf[0], xf[1], xf[2]}, {xf[3], xf[4], xf[5]}}};
    double mass_bck = xf[6];
    // We propagate for a first dt/2
    rv_bck = propagate_lagrangian(rv_bck, -dt / 2, mu, false).first;
    // We now loop through the backward segments and 1) add a dv + 2) propagate for -dt (except on the last segment,
    // where we propagate for -dt/2).
    for (decltype(m_throttles.size()) i = 0u; i < n_seg_bck; ++i) {
        // We compute the the dv
        dv[0] = max_thrust / mass_bck * dt * m_throttles[3 * (i + n_seg_fwd)];
        dv[1] = max_thrust / mass_bck * dt * m_throttles[3 * (i + n_seg_fwd) + 1];
        dv[2] = max_thrust / mass_bck * dt * m_throttles[3 * (i + n_seg_fwd) + 2];
        // Subtract it (remember we are going backward) to the current spacecraft velocity
        rv_bck[1][0] -= dv[0];
        rv_bck[1][1] -= dv[1];
        rv_bck[1][2] -= dv[2];
        // Update the mass accordingly
        double norm_dv = std::sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
        mass_bck *= std::exp(norm_dv / isp / kep3::G0);
        // Perform the propagation
        double prop_duration = (i == n_seg_bck - 1) ? -dt / 2 : -dt;
        rv_bck = propagate_lagrangian(rv_bck, prop_duration, mu, false).first;
    }
    return {rv_fwd[0][0] - rv_bck[0][0], rv_fwd[0][1] - rv_bck[0][1], rv_fwd[0][2] - rv_bck[0][2],
            rv_fwd[1][0] - rv_bck[1][0], rv_fwd[1][1] - rv_bck[1][1], rv_fwd[1][2] - rv_bck[1][2],
            mass_fwd - mass_bck};
}

std::vector<double> sims_flanagan::compute_throttle_constraints() const
{
    auto n_seg = m_throttles.size() / 3u;
    std::vector<double> retval(n_seg);
    for (decltype(m_throttles.size()) i = 0u; i < n_seg; ++i) {
        retval[i] = m_throttles[3 * i] * m_throttles[3 * i] + m_throttles[3 * i + 1] * m_throttles[3 * i + 1]
                    + m_throttles[3 * i + 2] * m_throttles[3 * i + 2] - 1.;
    }
    return retval;
}

std::ostream &operator<<(std::ostream &s, const sims_flanagan &sf)
{
    auto n_seg = sf.get_throttles().size() / 3u;
    auto n_seg_fwd = static_cast<unsigned>(static_cast<double>(n_seg) * sf.get_cut());
    auto n_seg_bck = n_seg - n_seg_fwd;
    s << fmt::format("Number of segments: {}\n", n_seg);
    s << fmt::format("Number of fwd segments: {}\n", n_seg_fwd);
    s << fmt::format("Number of bck segments: {}\n", n_seg_bck);
    s << fmt::format("Maximum thrust: {}\n", sf.get_max_thrust());
    s << fmt::format("Central body gravitational parameter: {}\n", sf.get_mu());
    s << fmt::format("Specific impulse: {}\n\n", sf.get_isp());
    s << fmt::format("Time of flight: {}\n", sf.get_tof());
    s << fmt::format("Initial mass: {}\n", sf.get_xs()[6]);
    s << fmt::format("Final mass: {}\n", sf.get_xf()[6]);
    s << fmt::format("State at departure: {}\n", sf.get_xs());
    s << fmt::format("State at arrival: {}\n", sf.get_xf());
    s << fmt::format("Throttles values: {}\n\n", sf.get_throttles());
    s << fmt::format("Mismatch constraints: {}\n", sf.compute_mismatch_constraints());
    s << fmt::format("Throttle constraints: {}\n\n", sf.compute_throttle_constraints());

    return s;
}

} // namespace kep3::leg