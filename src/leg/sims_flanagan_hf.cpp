// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <boost/range/algorithm.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xaxis_iterator.hpp>
#include <xtensor/xaxis_slice_iterator.hpp>
#include <xtensor/xbuilder.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/epoch.hpp>
#include <kep3/leg/sf_checks.hpp>
#include <kep3/leg/sims_flanagan_hf.hpp>
#include <kep3/linalg.hpp>
#include <kep3/ta/stark.hpp>

#include <heyoka/taylor.hpp>

namespace kep3::leg
{

using kep3::linalg::mat13;
using kep3::linalg::mat61;
using kep3::linalg::mat63;
using kep3::linalg::mat66;

// Constructors

sims_flanagan_hf::sims_flanagan_hf()
{
    // We perform some sanity checks on the user provided inputs
    _sanity_checks(m_throttles, m_tof, m_max_thrust, m_isp, m_mu, m_cut, m_tol, m_nseg, m_nseg_fwd, m_nseg_bck);
    // We set mu and veff for the non variational
    *m_tas.get_pars_data() = m_mu;
    *(m_tas.get_pars_data() + 1) = m_isp * kep3::G0;

    // ... and variational version of the integrator
    *(m_tas_var.get_pars_data()) = m_mu;
    *(m_tas_var.get_pars_data() + 1) = m_isp * kep3::G0;
    // We copy the initial conditions for the variational equations
    std::copy(m_tas_var.get_state().begin() + 7, m_tas_var.get_state().end(), m_vars.begin());

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}

sims_flanagan_hf::sims_flanagan_hf(const std::array<std::array<double, 3>, 2> &rvs, double ms,
                                   std::vector<double> throttles,
                                   // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                                   const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof,
                                   double max_thrust, double isp, double mu, double cut, double tol)
    : m_throttles(std::move(throttles)), m_tof(tof), m_max_thrust(max_thrust), m_isp(isp), m_mu(mu), m_cut(cut),
      m_tol(tol), m_nseg(static_cast<unsigned>(m_throttles.size()) / 3u),
      m_nseg_fwd(static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut)), m_nseg_bck(m_nseg - m_nseg_fwd)
{
    // We perform some sanity checks on the user provided inputs
    _sanity_checks(m_throttles, m_tof, m_max_thrust, m_isp, m_mu, m_cut, m_tol, m_nseg, m_nseg_fwd, m_nseg_bck);
    // We set mu and veff for the non variational
    *m_tas.get_pars_data() = m_mu;
    *(m_tas.get_pars_data() + 1) = m_isp * kep3::G0;

    // ... and variational version of the integrator
    *(m_tas_var.get_pars_data()) = m_mu;
    *(m_tas_var.get_pars_data() + 1) = m_isp * kep3::G0;
    // We copy the initial conditions for the variational equations
    std::copy(m_tas_var.get_state().begin() + 7, m_tas_var.get_state().end(), m_vars.begin());

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
    // Fill in m_rvm from m_rvs and m_ms
    std::copy(rvs[0].begin(), rvs[0].end(), m_rvms.begin());
    std::copy(rvs[1].begin(), rvs[1].end(), std::next(m_rvms.begin(), 3));
    set_ms(ms);
    // Fill in m_rvm from m_rvf and m_mf
    std::copy(rvf[0].begin(), rvf[0].end(), m_rvmf.begin());
    std::copy(rvf[1].begin(), rvf[1].end(), std::next(m_rvmf.begin(), 3));
    set_mf(mf);
}

sims_flanagan_hf::sims_flanagan_hf(const std::array<double, 7> &rvms, std::vector<double> throttles,
                                   const std::array<double, 7> &rvmf, double tof, double max_thrust, double isp,
                                   double mu, double cut, double tol)
    : m_rvms(rvms), m_throttles(std::move(throttles)), m_rvmf(rvmf), m_tof(tof), m_max_thrust(max_thrust), m_isp(isp),
      m_mu(mu), m_cut(cut), m_tol(tol), m_nseg(static_cast<unsigned>(m_throttles.size()) / 3u),
      m_nseg_fwd(static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut)), m_nseg_bck(m_nseg - m_nseg_fwd)
{
    // We perform some sanity checks on the user provided inputs
    _sanity_checks(m_throttles, m_tof, m_max_thrust, m_isp, m_mu, m_cut, m_tol, m_nseg, m_nseg_fwd, m_nseg_bck);
    // We set mu and veff for the non variational
    *m_tas.get_pars_data() = m_mu;
    *(m_tas.get_pars_data() + 1) = m_isp * kep3::G0;

    // ... and variational version of the integrator
    *(m_tas_var.get_pars_data()) = m_mu;
    *(m_tas_var.get_pars_data() + 1) = m_isp * kep3::G0;
    // We copy the initial conditions for the variational equations
    std::copy(m_tas_var.get_state().begin() + 7, m_tas_var.get_state().end(), m_vars.begin());

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}

// Setters
void sims_flanagan_hf::set_tof(double tof)
{
    _check_tof(tof);
    m_tof = tof;
}
void sims_flanagan_hf::set_rvs(std::array<std::array<double, 3>, 2> rv)
{
    std::copy(rv[0].begin(), rv[0].end(), m_rvms.begin());
    std::copy(rv[1].begin(), rv[1].end(), std::next(m_rvms.begin(), 3));
}
void sims_flanagan_hf::set_ms(double mass)
{
    m_rvms[6] = mass;
}
void sims_flanagan_hf::set_throttles(std::vector<double> throttles)
{
    _check_throttles(throttles, m_nseg);
    m_throttles = std::move(throttles);
    m_nseg = static_cast<unsigned>(m_throttles.size()) / 3u;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
}
void sims_flanagan_hf::set_throttles(std::vector<double>::const_iterator it1, std::vector<double>::const_iterator it2)
{
    if (((std::distance(it1, it2) % 3) != 0) || std::distance(it1, it2) <= 0) {
        throw std::logic_error(
            "The throttles of a high-fidelity sims-flanagan leg leg are being set with invalid iterators.");
    }
    m_throttles.resize(static_cast<size_t>(std::distance(it1, it2)));
    std::copy(it1, it2, m_throttles.begin());
    m_nseg = static_cast<unsigned>(m_throttles.size()) / 3u;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
}
void sims_flanagan_hf::set_rvf(std::array<std::array<double, 3>, 2> rv)
{
    std::copy(rv[0].begin(), rv[0].end(), m_rvmf.begin());
    std::copy(rv[1].begin(), rv[1].end(), std::next(m_rvmf.begin(), 3));
}
void sims_flanagan_hf::set_mf(double mass)
{
    m_rvmf[6] = mass;
}
void sims_flanagan_hf::set_max_thrust(double max_thrust)
{
    _check_max_thrust(max_thrust);
    m_max_thrust = max_thrust;
}
void sims_flanagan_hf::set_isp(double isp)
{
    _check_isp(isp);
    m_isp = isp;
}
void sims_flanagan_hf::set_mu(double mu)
{
    _check_mu(mu);
    m_mu = mu;
}
void sims_flanagan_hf::set_cut(double cut)
{
    _check_cut(cut);
    m_cut = cut;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
}
void sims_flanagan_hf::set_tol(double tol)
{
    _check_tol(tol);
    m_tol = tol;
}
void sims_flanagan_hf::set_rvms(std::array<double, 7> rvms)
{
    m_rvms = rvms;
}
void sims_flanagan_hf::set_rvmf(std::array<double, 7> rvmf)
{
    m_rvmf = rvmf;
}

void sims_flanagan_hf::set(const std::array<std::array<double, 3>, 2> &rvs, double ms,
                           const std::vector<double> &throttles,
                           // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                           const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust,
                           double isp, double mu, double cut, double tol)
{
    _sanity_checks(throttles, tof, max_thrust, isp, mu, cut, tol, m_nseg, m_nseg_fwd, m_nseg_bck);
    // Set initial state
    set_rvs(rvs);
    set_ms(ms);
    // Set final state
    set_rvf(rvf);
    set_mf(mf);
    m_throttles = throttles;
    m_tof = tof;
    m_max_thrust = max_thrust;
    m_isp = isp;
    m_mu = mu;
    m_cut = cut;
    m_tol = tol;
    m_nseg = static_cast<unsigned>(m_throttles.size()) / 3u;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}

void sims_flanagan_hf::set(const std::array<double, 7> &rvms, const std::vector<double> &throttles,
                           const std::array<double, 7> &rvmf, double tof, double max_thrust, double isp, double mu,
                           double cut, double tol)
{
    _sanity_checks(throttles, tof, max_thrust, isp, mu, cut, tol, m_nseg, m_nseg_fwd, m_nseg_bck);
    set_rvms(rvms);
    m_throttles = throttles;
    set_rvmf(rvmf);
    m_tof = tof;
    m_max_thrust = max_thrust;
    m_isp = isp;
    m_mu = mu;
    m_cut = cut;
    m_tol = tol;
    m_nseg = static_cast<unsigned>(m_throttles.size()) / 3u;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}

void sims_flanagan_hf::set(const std::array<double, 7> &rvms, const std::vector<double> &throttles,
                           const std::array<double, 7> &rvmf, double time_of_flight)
{
    set_rvms(rvms);
    m_throttles = throttles;
    set_rvmf(rvmf);
    m_tof = time_of_flight;
    m_nseg = static_cast<unsigned>(m_throttles.size()) / 3u;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
    _sanity_checks(throttles, m_tof, m_max_thrust, m_isp, m_mu, m_cut, m_tol, m_nseg, m_nseg_fwd, m_nseg_bck);

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}

// Getters
double sims_flanagan_hf::get_tof() const
{
    return m_tof;
}
const std::array<std::array<double, 3>, 2> sims_flanagan_hf::get_rvs() const
{
    std::array<std::array<double, 3>, 2> rvs{};
    std::copy(m_rvms.begin(), std::next(m_rvms.begin(), 3), rvs[0].begin());
    std::copy(std::next(m_rvms.begin(), 3), std::next(m_rvms.begin(), 6), rvs[1].begin());
    return rvs;
}
double sims_flanagan_hf::get_ms() const
{
    return m_rvms[6];
}
const std::vector<double> &sims_flanagan_hf::get_throttles() const
{
    return m_throttles;
}
const std::array<std::array<double, 3>, 2> sims_flanagan_hf::get_rvf() const
{
    std::array<std::array<double, 3>, 2> rvf{{{0., 0., 0.}, {0., 0., 0.}}};
    std::copy(m_rvmf.begin(), std::next(m_rvmf.begin(), 3), rvf[0].begin());
    std::copy(std::next(m_rvmf.begin(), 3), std::next(m_rvmf.begin(), 6), rvf[1].begin());
    return rvf;
}
double sims_flanagan_hf::get_mf() const
{
    return m_rvmf[6];
}
double sims_flanagan_hf::get_max_thrust() const
{
    return m_max_thrust;
}
double sims_flanagan_hf::get_isp() const
{
    return m_isp;
}
double sims_flanagan_hf::get_mu() const
{
    return m_mu;
}
double sims_flanagan_hf::get_cut() const
{
    return m_cut;
}
double sims_flanagan_hf::get_tol() const
{
    return m_tol;
}
unsigned sims_flanagan_hf::get_nseg() const
{
    return m_nseg;
}
unsigned sims_flanagan_hf::get_nseg_fwd() const
{
    return m_nseg_fwd;
}
unsigned sims_flanagan_hf::get_nseg_bck() const
{
    return m_nseg_bck;
}
heyoka::taylor_adaptive<double> sims_flanagan_hf::get_tas() const
{
    return m_tas;
}
heyoka::taylor_adaptive<double> sims_flanagan_hf::get_tas_var() const
{
    return m_tas_var;
}
std::array<double, 7> sims_flanagan_hf::get_rvms() const
{
    return m_rvms;
}
std::array<double, 7> sims_flanagan_hf::get_rvmf() const
{
    return m_rvmf;
}

// The core routines
std::array<double, 7> sims_flanagan_hf::compute_mismatch_constraints()
{
    // General settings
    const double prop_seg_duration = (m_tof / m_nseg);

    // Forward pass
    // Initial state
    // Set the Taylor Integration initial conditions
    m_tas.set_time(0.);
    std::copy(m_rvms.begin(), m_rvms.end(), m_tas.get_state_data());

    // Loop through segments in forward pass of Sims-Flanagan transcription
    for (unsigned int i = 0u; i < m_nseg_fwd; ++i) {
        // Assign current thrusts to Taylor adaptive integrator
        if (static_cast<size_t>((i + 1) * 3) <= m_thrusts.size()) {
            std::copy(std::next(m_thrusts.begin(), static_cast<long>(i * 3)),
                      std::next(m_thrusts.begin(), static_cast<long>(3 * (i + 1))),
                      std::next(m_tas.get_pars_data(), 2));
        } else {
            throw std::runtime_error("The retrieved thrust index is larger than the size of the m_thrusts vector.");
        }
        // ... and integrate
        auto [status, min_h, max_h, nsteps, _1, _2] = m_tas.propagate_until((i + 1) * prop_seg_duration);
        if (status != heyoka::taylor_outcome::time_limit) {
            throw std::domain_error("stark_problem: failure to reach the final time requested during a propagation.");
        }
    }

    // Set fwd final state
    std::vector<double> rvm_fwd_final = m_tas.get_state();

    // Backward pass
    // Final state
    // Set the Taylor Integration final conditions
    m_tas.set_time(m_tof);
    std::copy(m_rvmf.begin(), m_rvmf.end(), m_tas.get_state_data());

    // Loop through segments in backward pass of Sims-Flanagan transcription
    for (unsigned int i = 0u; i < m_nseg_bck; ++i) {
        // Assign current_thrusts to Taylor adaptive integrator
        if (static_cast<size_t>((m_nseg - i) * 3) <= m_thrusts.size()) {
            // Copy thrust into Taylor-adaptive integrator
            std::copy(std::next(m_thrusts.begin(), static_cast<long>((m_nseg - (i + 1)) * 3)),
                      std::next(m_thrusts.begin(), static_cast<long>((m_nseg - i) * 3)),
                      std::next(m_tas.get_pars_data(), 2));
        } else {
            throw std::runtime_error("The retrieved thrust index is larger than the size of the m_thrusts vector.");
        }
        // ... and integrate
        auto [status, min_h, max_h, nsteps, _1, _2] = m_tas.propagate_until(m_tof - (i + 1) * prop_seg_duration);
        if (status != heyoka::taylor_outcome::time_limit) {
            throw std::domain_error("stark_problem: failure to reach the final time requested during a propagation.");
        }
    }

    return {rvm_fwd_final[0] - m_tas.get_state()[0], rvm_fwd_final[1] - m_tas.get_state()[1],
            rvm_fwd_final[2] - m_tas.get_state()[2], rvm_fwd_final[3] - m_tas.get_state()[3],
            rvm_fwd_final[4] - m_tas.get_state()[4], rvm_fwd_final[5] - m_tas.get_state()[5],
            rvm_fwd_final[6] - m_tas.get_state()[6]};
}

std::vector<double> sims_flanagan_hf::compute_throttle_constraints() const
{
    std::vector<double> retval(m_nseg);
    for (decltype(m_throttles.size()) i = 0u; i < m_nseg; ++i) {
        retval[i] = m_throttles[3 * i] * m_throttles[3 * i] + m_throttles[3 * i + 1] * m_throttles[3 * i + 1]
                    + m_throttles[3 * i + 2] * m_throttles[3 * i + 2] - 1.;
    }
    return retval;
}

std::vector<double> sims_flanagan_hf::compute_constraints()
{
    std::vector<double> retval(7 + m_nseg, 0.);
    // Fitness
    // Equality Constraints
    auto eq_con = compute_mismatch_constraints();
    retval[0] = eq_con[0];
    retval[1] = eq_con[1];
    retval[2] = eq_con[2];
    retval[3] = eq_con[3];
    retval[4] = eq_con[4];
    retval[5] = eq_con[5];
    retval[6] = eq_con[6];
    //  Inequality Constraints
    auto ineq_con = compute_throttle_constraints();
    std::copy(ineq_con.begin(), ineq_con.end(), retval.begin() + 7);
    return retval;
}

std::vector<double> sims_flanagan_hf::set_and_compute_constraints(std::vector<double> chromosome)
{
    std::array<double, 7> rvms;
    std::copy(chromosome.begin(), chromosome.begin() + 7, rvms.begin());
    std::vector<double> throttles(m_nseg * 3);
    std::copy(chromosome.begin() + 7, chromosome.begin() + 7 + m_nseg * 3, throttles.begin());
    std::array<double, 7> rvmf;
    std::copy(chromosome.begin() + 7 + m_nseg * 3, chromosome.begin() + 7 + m_nseg * 3 + 7, rvmf.begin());
    double time_of_flight = chromosome[29];
    // Set relevant quantities before evaluating constraints
    set(rvms, throttles, rvmf, time_of_flight);
    // Evaluate and return constraints
    return compute_constraints();
}

// Return specific two-body 'stark' dynamics state derivative
std::array<double, 7> sims_flanagan_hf::get_state_derivative(std::array<double, 7> state,
                                                             std::array<double, 3> throttles)
{

    std::array<double, 3> thrusts;
    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    std::transform(throttles.begin(), throttles.end(), thrusts.begin(), throttle_to_thrust);

    std::array<double, 7> dstatedt;
    // The square of the radius
    std::array<double, 3> state_squared = {std::pow(state[0], 2.), std::pow(state[1], 2.), std::pow(state[2], 2.)};
    const auto r2 = std::accumulate(state_squared.begin(), state_squared.end(), 0.0);
    double veff = get_isp() * kep3::G0;

    // The throttle magnitude
    std::array<double, 3> thrusts_squared
        = {std::pow(thrusts[0], 2.), std::pow(thrusts[1], 2.), std::pow(thrusts[2], 2.)};
    const auto u_norm = std::sqrt(std::accumulate(thrusts_squared.begin(), thrusts_squared.end(), 0.0));

    // The Equations of Motion
    dstatedt[0] = state[3];
    dstatedt[1] = state[4];
    dstatedt[2] = state[5];
    dstatedt[3] = -get_mu() * pow(r2, -3. / 2) * state[0] + thrusts[0] / state[6];
    dstatedt[4] = -get_mu() * pow(r2, -3. / 2) * state[1] + thrusts[1] / state[6];
    dstatedt[5] = -get_mu() * pow(r2, -3. / 2) * state[2] + thrusts[2] / state[6];
    dstatedt[6] = (u_norm != 0) ? -u_norm / veff : 0; // Conditional for if thrust is zero or not

    return dstatedt;
}

std::tuple<std::array<std::array<double, 7u>, 5u>, std::array<std::array<double, 49u>, 5u>,
           std::array<std::array<double, 21u>, 5u>>
sims_flanagan_hf::compute_mc_grad()
{
    // Initialise
    std::array<std::array<double, 7u>, 5u> xf_per_seg = {{{0}}};
    std::array<std::array<double, 49u>, 5u> dxdx_per_seg = {{{0}}};
    std::array<std::array<double, 21u>, 5u> dxdu_per_seg = {{{0}}};

    // General settings
    const double prop_seg_duration = (m_tof / m_nseg);

    // Forward loop
    // Set the Taylor Integration initial conditions
    m_tas_var.set_time(0.);
    std::copy(m_rvms.begin(), m_rvms.end(), m_tas_var.get_state_data());

    for (unsigned int i = 0u; i < m_nseg_fwd; ++i) {

        // Initialise var conditions
        std::copy(m_vars.begin(), m_vars.end(), m_tas_var.get_state_data() + 7);
        // Assign current thrusts to Taylor adaptive integrator
        if (static_cast<size_t>((i + 1) * 3) <= m_thrusts.size()) {
            std::copy(std::next(m_thrusts.begin(), static_cast<long>(i * 3)),
                      std::next(m_thrusts.begin(), static_cast<long>(3 * (i + 1))),
                      std::next(m_tas_var.get_pars_data(), 2));
        } else {
            throw std::runtime_error("The retrieved thrust index is larger than the size of the m_thrusts vector.");
        }
        // ... and integrate
        auto [status, min_h, max_h, nsteps, _1, _2] = m_tas_var.propagate_until((i + 1) * prop_seg_duration);
        if (status != heyoka::taylor_outcome::time_limit) {
            throw std::domain_error("stark_problem: failure to reach the final time requested during a propagation.");
        }
        // Save the variational state variables to respective arrays
        std::copy(m_tas_var.get_state().begin(), m_tas_var.get_state().begin() + 7, xf_per_seg[i].begin());
        for (auto j = 0; j < 7; ++j) {
            std::copy(std::next(m_tas_var.get_state().begin(), 7 + 10l * j),
                      std::next(m_tas_var.get_state().begin(), 7 + 10l * j + 7),
                      std::next(dxdx_per_seg[i].begin(), 7 * j));
            std::copy(m_tas_var.get_state().begin() + 14 + 10l * j, m_tas_var.get_state().begin() + 14 + 10l * j + 3,
                      dxdu_per_seg[i].begin() + 3l * j);
        }
    }

    // Backward loop
    // Set the Taylor Integration initial conditions
    m_tas_var.set_time(m_tof);
    std::copy(m_rvmf.begin(), m_rvmf.end(), m_tas_var.get_state_data());

    for (unsigned int i = 0u; i < m_nseg_bck; ++i) {

        // Initialise var conditions
        std::copy(m_vars.begin(), m_vars.end(), m_tas_var.get_state_data() + 7);
        // Assign current thrusts to Taylor adaptive integrator
        if (static_cast<size_t>((m_nseg - i) * 3) <= m_thrusts.size()) {
            // Copy thrust into Taylor-adaptive integrator
            std::copy(std::next(m_thrusts.begin(), static_cast<long>((m_nseg - (i + 1)) * 3)),
                      std::next(m_thrusts.begin(), static_cast<long>((m_nseg - i) * 3)),
                      std::next(m_tas_var.get_pars_data(), 2));
        } else {
            throw std::runtime_error("The retrieved thrust index is larger than the size of the m_thrusts vector.");
        }
        // ... and integrate
        auto [status, min_h, max_h, nsteps, _1, _2] = m_tas_var.propagate_until(m_tof - (i + 1) * prop_seg_duration);
        if (status != heyoka::taylor_outcome::time_limit) {
            throw std::domain_error("stark_problem: failure to reach the final time requested during a propagation.");
        }
        // Save the variational state variables to respective arrays
        std::copy(m_tas_var.get_state().begin(), m_tas_var.get_state().begin() + 7,
                  xf_per_seg[m_nseg - (i + 1)].begin());
        for (auto j = 0; j < 7; ++j) {
            std::copy(m_tas_var.get_state().begin() + 7 + 10l * j, m_tas_var.get_state().begin() + 7 + 10l * j + 7,
                      dxdx_per_seg[m_nseg - (i + 1)].begin() + 7 * j);
            std::copy(m_tas_var.get_state().begin() + 14 + 10l * j, m_tas_var.get_state().begin() + 14 + 10l * j + 3,
                      dxdu_per_seg[m_nseg - (i + 1)].begin() + 3l * j);
        }
    }

    return std::make_tuple(xf_per_seg, dxdx_per_seg, dxdu_per_seg);
}

std::vector<double> sims_flanagan_hf::compute_tc_grad() const
{
    std::vector<double> retval(static_cast<size_t>(m_nseg) * m_nseg * 3u, 0);
    for (decltype(m_throttles.size()) i = 0u; i < m_nseg; ++i) {
        retval[i * m_nseg * 3 + 3 * i] = 2 * m_throttles[3 * i];
        retval[i * m_nseg * 3 + 3 * i + 1] = 2 * m_throttles[3 * i + 1];
        retval[i * m_nseg * 3 + 3 * i + 2] = 2 * m_throttles[3 * i + 2];
    }
    return retval;
}

std::ostream &operator<<(std::ostream &s, const sims_flanagan_hf &sf)
{
    s << fmt::format("Number of segments: {}\n", sf.get_nseg());
    s << fmt::format("Number of fwd segments: {}\n", sf.get_nseg_fwd());
    s << fmt::format("Number of bck segments: {}\n", sf.get_nseg_bck());
    s << fmt::format("Maximum thrust: {}\n", sf.get_max_thrust());
    s << fmt::format("Central body gravitational parameter: {}\n", sf.get_mu());
    s << fmt::format("Specific impulse: {}\n\n", sf.get_isp());
    s << fmt::format("Time of flight: {}\n", sf.get_tof());
    s << fmt::format("Initial mass: {}\n", sf.get_ms());
    s << fmt::format("Final mass: {}\n", sf.get_mf());
    s << fmt::format("State at departure: {}\n", sf.get_rvs());
    s << fmt::format("State at arrival: {}\n", sf.get_rvf());
    s << fmt::format("Throttles values: {}\n\n", sf.get_throttles());
    // s << fmt::format("Mismatch constraints: {}\n", sf.compute_mismatch_constraints());
    s << fmt::format("Throttle constraints: {}\n\n", sf.compute_throttle_constraints());
    return s;
}

} // namespace kep3::leg