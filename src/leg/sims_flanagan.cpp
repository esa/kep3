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
#include <iterator>
#include <vector>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <xtensor/xadapt.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/epoch.hpp>
#include <kep3/leg/sims_flanagan.hpp>
#include <kep3/linalg.hpp>

namespace kep3::leg
{

//using kep3::linalg::_dot;

void _check_tof(double tof)
{
    if (tof < 0.) {
        throw std::domain_error("The time of flight of a sims_flanagan leg needs to be larger or equal to zero.");
    }
}
void _check_throttles(const std::vector<double> &throttles)
{
    if ((throttles.size() % 3) != 0u) {
        throw std::logic_error("The throttles of a sims_flanagan leg are detected to be not a multiple of 3 in size "
                               "[u0x, u0y, u0z, .....].");
    }
    if (throttles.empty()) {
        throw std::logic_error(
            "The throttles of a sims_flanagan leg are detected to be empty! At least one segment is necessary.");
    }
}
void _check_max_thrust(double max_thrust)
{
    if (max_thrust < 0.) {
        throw std::domain_error(
            "The maximum allowed thrust of a sims_flanagan leg is detected to be smaller than zero.");
    }
}
void _check_isp(double isp)
{
    if (isp < 0.) {
        throw std::domain_error("The specific impulse of a sims_flanagan leg is detected to be smaller than zero.");
    }
}
void _check_mu(double mu)
{
    if (mu < 0.) {
        throw std::domain_error(
            "The gravitational parameter of a sims_flanagan leg is detected to be smaller than zero.");
    }
}
void _check_cut(double cut)
{
    if (cut < 0. || cut > 1.) {
        throw std::domain_error("The parameter cut of a sims_flanagan leg must be in [0, 1].");
    }
}
void _sanity_checks(const std::array<std::array<double, 3>, 2> &, double, const std::vector<double> &throttles,
                    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                    const std::array<std::array<double, 3>, 2> &, double, double tof, double max_thrust, double isp,
                    double mu, double cut)
{
    _check_throttles(throttles);
    _check_tof(tof);
    _check_max_thrust(max_thrust);
    _check_isp(isp);
    _check_mu(mu);
    _check_cut(cut);
}

sims_flanagan::sims_flanagan(const std::array<std::array<double, 3>, 2> &rvs, double ms, std::vector<double> throttles,
                             // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                             const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust,
                             double isp, double mu, double cut)
    : m_rvs(rvs), m_ms(ms), m_throttles(std::move(throttles)), m_rvf(rvf), m_mf(mf), m_tof(tof),
      m_max_thrust(max_thrust), m_isp(isp), m_mu(mu), m_cut(cut)
{
    _sanity_checks(rvs, ms, m_throttles, rvf, mf, tof, max_thrust, isp, mu, cut);
}

// Setters
void sims_flanagan::set_tof(double tof)
{
    _check_tof(tof);
    m_tof = tof;
}
void sims_flanagan::set_rvs(std::array<std::array<double, 3>, 2> rv)
{
    m_rvs = rv;
}
void sims_flanagan::set_ms(double mass)
{
    m_ms = mass;
}
void sims_flanagan::set_throttles(std::vector<double> throttles)
{
    _check_throttles(throttles);
    m_throttles = std::move(throttles);
}
void sims_flanagan::set_throttles(std::vector<double>::const_iterator it1, std::vector<double>::const_iterator it2)
{
    if (((std::distance(it1, it2) % 3) != 0) || std::distance(it1, it2) <= 0) {
        throw std::logic_error("The throttles of a sims_flanagan leg are being set with invalid iterators.");
    }
    std::copy(it1, it2, m_throttles.begin());
}
void sims_flanagan::set_rvf(std::array<std::array<double, 3>, 2> rv)
{
    m_rvf = rv;
}
void sims_flanagan::set_mf(double mass)
{
    m_mf = mass;
}
void sims_flanagan::set_max_thrust(double max_thrust)
{
    _check_max_thrust(max_thrust);
    m_max_thrust = max_thrust;
}
void sims_flanagan::set_isp(double isp)
{
    _check_isp(isp);
    m_isp = isp;
}
void sims_flanagan::set_mu(double mu)
{
    _check_mu(mu);
    m_mu = mu;
}
void sims_flanagan::set_cut(double cut)
{
    _check_cut(cut);
    m_cut = cut;
}
void sims_flanagan::set(const std::array<std::array<double, 3>, 2> &rvs, double ms,
                        const std::vector<double> &throttles,
                        // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                        const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust,
                        double isp, double mu, double cut)
{
    _sanity_checks(rvs, ms, throttles, rvf, mf, tof, max_thrust, isp, mu, cut);
    m_rvs = rvs;
    m_ms = ms;
    m_throttles = throttles;
    m_rvf = rvf;
    m_mf = mf;
    m_tof = tof;
    m_max_thrust = max_thrust;
    m_isp = isp;
    m_mu = mu;
    m_cut = cut;
}

// Getters
double sims_flanagan::get_tof() const
{
    return m_tof;
}
const std::array<std::array<double, 3>, 2> &sims_flanagan::get_rvs() const
{
    return m_rvs;
}
double sims_flanagan::get_ms() const
{
    return m_ms;
}
const std::vector<double> &sims_flanagan::get_throttles() const
{
    return m_throttles;
}
const std::array<std::array<double, 3>, 2> &sims_flanagan::get_rvf() const
{
    return m_rvf;
}
double sims_flanagan::get_mf() const
{
    return m_mf;
}
double sims_flanagan::get_max_thrust() const
{
    return m_max_thrust;
}
double sims_flanagan::get_isp() const
{
    return m_isp;
}
double sims_flanagan::get_mu() const
{
    return m_mu;
}
double sims_flanagan::get_cut() const
{
    return m_cut;
}

// The core routines
std::array<double, 7> sims_flanagan::compute_mismatch_constraints() const
{
    // We start defining the number of forward and backward segments.
    auto nseg = m_throttles.size() / 3u;
    auto nseg_fwd = static_cast<unsigned>(static_cast<double>(nseg) * m_cut);
    auto nseg_bck = nseg - nseg_fwd;

    // We introduce some convenience variables
    double max_thrust = get_max_thrust();
    double isp = get_isp();
    double mu = get_mu();
    std::array<double, 3> dv{};

    // Forward pass
    // Initial state
    std::array<std::array<double, 3>, 2> rv_fwd(get_rvs());
    double mass_fwd = get_mf();
    double dt = m_tof / static_cast<double>(nseg);
    // We propagate for a first dt/2 (only if there is at least one forward segment)
    if (nseg_fwd > 0) {
        rv_fwd = propagate_lagrangian(rv_fwd, dt / 2, mu, false).first;
    }
    // We now loop through the forward segments and 1) add a dv + 2) propagate for dt (except on the last segment, where
    // we propagate for dt/2).
    for (decltype(m_throttles.size()) i = 0u; i < nseg_fwd; ++i) {
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
        double prop_duration = (i == nseg_fwd - 1) ? dt / 2 : dt;
        rv_fwd = propagate_lagrangian(rv_fwd, prop_duration, mu, false).first;
    }

    // Backward pass
    // Final state
    std::array<std::array<double, 3>, 2> rv_bck(get_rvf());
    double mass_bck = get_mf();
    // We propagate for a first dt/2 (only if there is at least one backward segment)
    if (nseg_bck > 0) {
        rv_bck = propagate_lagrangian(rv_bck, -dt / 2, mu, false).first;
    }
    // We now loop through the backward segments and 1) add a dv + 2) propagate for -dt (except on the last segment,
    // where we propagate for -dt/2).
    for (decltype(m_throttles.size()) i = 0u; i < nseg_bck; ++i) {
        // We compute the the dv
        dv[0] = max_thrust / mass_bck * dt * m_throttles[m_throttles.size() - 1 - 3 * i - 2];
        dv[1] = max_thrust / mass_bck * dt * m_throttles[m_throttles.size() - 1 - 3 * i - 1];
        dv[2] = max_thrust / mass_bck * dt * m_throttles[m_throttles.size() - 1 - 3 * i];
        // Subtract it (remember we are going backward) to the current spacecraft velocity
        rv_bck[1][0] -= dv[0];
        rv_bck[1][1] -= dv[1];
        rv_bck[1][2] -= dv[2];
        // Update the mass accordingly (will increase as we go backward)
        double norm_dv = std::sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
        mass_bck *= std::exp(norm_dv / isp / kep3::G0);
        // Perform the propagation
        double prop_duration = (i == nseg_bck - 1) ? -dt / 2 : -dt;
        rv_bck = propagate_lagrangian(rv_bck, prop_duration, mu, false).first;
    }
    return {rv_fwd[0][0] - rv_bck[0][0], rv_fwd[0][1] - rv_bck[0][1], rv_fwd[0][2] - rv_bck[0][2],
            rv_fwd[1][0] - rv_bck[1][0], rv_fwd[1][1] - rv_bck[1][1], rv_fwd[1][2] - rv_bck[1][2],
            mass_fwd - mass_bck};
}

std::vector<double> sims_flanagan::compute_throttle_constraints() const
{
    auto nseg = m_throttles.size() / 3u;

    std::vector<double> retval(nseg);
    for (decltype(m_throttles.size()) i = 0u; i < nseg; ++i) {
        retval[i] = m_throttles[3 * i] * m_throttles[3 * i] + m_throttles[3 * i + 1] * m_throttles[3 * i + 1]
                    + m_throttles[3 * i + 2] * m_throttles[3 * i + 2] - 1.;
    }
    return retval;
}

std::pair<std::array<double, 49>, std::vector<double>> sims_flanagan::compute_mc_grad() const
{
    // Preliminaries
    auto nseg = m_throttles.size() / 3u;
    //auto nseg_fwd = static_cast<unsigned>(static_cast<double>(nseg) * m_cut);
    //auto nseg_bck = nseg - nseg_fwd;
    //auto c = m_max_thrust * m_tof / static_cast<double>(nseg); // T*tof/nseg
    //auto a = 1. / m_isp / kep3::G0;                            // 1/veff
    //auto dt = m_tof / static_cast<double>(nseg);               // dt

    // Allocate the return values
    std::array<double, 49> grad_rvm{};          // The mismatch constraints gradient w.r.t. extended state r,v,m
    std::vector<double> grad(nseg * 3 + 1, 0.); // The mismatch constraints gradient w.r.t. throttles and tof
    return {grad_rvm, grad};
}

std::ostream &operator<<(std::ostream &s, const sims_flanagan &sf)
{
    auto nseg = sf.get_throttles().size() / 3u;
    auto nseg_fwd = static_cast<unsigned>(static_cast<double>(nseg) * sf.get_cut());
    auto nseg_bck = nseg - nseg_fwd;
    s << fmt::format("Number of segments: {}\n", nseg);
    s << fmt::format("Number of fwd segments: {}\n", nseg_fwd);
    s << fmt::format("Number of bck segments: {}\n", nseg_bck);
    s << fmt::format("Maximum thrust: {}\n", sf.get_max_thrust());
    s << fmt::format("Central body gravitational parameter: {}\n", sf.get_mu());
    s << fmt::format("Specific impulse: {}\n\n", sf.get_isp());
    s << fmt::format("Time of flight: {}\n", sf.get_tof());
    s << fmt::format("Initial mass: {}\n", sf.get_ms());
    s << fmt::format("Final mass: {}\n", sf.get_mf());
    s << fmt::format("State at departure: {}\n", sf.get_rvs());
    s << fmt::format("State at arrival: {}\n", sf.get_rvf());
    s << fmt::format("Throttles values: {}\n\n", sf.get_throttles());
    s << fmt::format("Mismatch constraints: {}\n", sf.compute_mismatch_constraints());
    s << fmt::format("Throttle constraints: {}\n\n", sf.compute_throttle_constraints());
    return s;
}

} // namespace kep3::leg