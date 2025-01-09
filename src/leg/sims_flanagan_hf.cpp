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
#include <sys/types.h>
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

// Constructors

sims_flanagan_hf::sims_flanagan_hf()
{
    // We perform some sanity checks on the user provided inputs
    kep3::leg::_sanity_checks(m_throttles, m_tof, m_max_thrust, m_isp, m_mu, m_cut, m_tol, m_nseg, m_nseg_fwd,
                              m_nseg_bck);

    // Initialize m_tas and m_tas_var
    const heyoka::taylor_adaptive<double> ta_cache = kep3::ta::get_ta_stark(m_tol);
    m_tas = ta_cache;
    const heyoka::taylor_adaptive<double> ta_var_cache = kep3::ta::get_ta_stark_var(m_tol);
    m_tas_var = ta_var_cache;

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
                                   const std::vector<double> &throttles,
                                   // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                                   const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof,
                                   double max_thrust, double isp, double mu, double cut, double tol)
    : m_throttles(throttles), m_tof(tof), m_max_thrust(max_thrust), m_isp(isp), m_mu(mu), m_cut(cut),
      m_tol(tol), m_nseg(static_cast<unsigned>(m_throttles.size()) / 3u),
      m_nseg_fwd(static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut)), m_nseg_bck(m_nseg - m_nseg_fwd)
{
    // We perform some sanity checks on the user provided inputs
    kep3::leg::_sanity_checks(m_throttles, m_tof, m_max_thrust, m_isp, m_mu, m_cut, m_tol, m_nseg, m_nseg_fwd,
                              m_nseg_bck);

    // Initialize m_tas and m_tas_var
    const heyoka::taylor_adaptive<double> ta_cache = kep3::ta::get_ta_stark(m_tol);
    m_tas = ta_cache;
    const heyoka::taylor_adaptive<double> ta_var_cache = kep3::ta::get_ta_stark_var(m_tol);
    m_tas_var = ta_var_cache;

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

sims_flanagan_hf::sims_flanagan_hf(const std::array<double, 7> &rvms, const std::vector<double> &throttles,
                                   // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                                   const std::array<double, 7> &rvmf, double tof, double max_thrust, double isp,
                                   double mu, double cut, double tol)
    : m_rvms(rvms), m_throttles(throttles), m_rvmf(rvmf), m_tof(tof), m_max_thrust(max_thrust), m_isp(isp),
      m_mu(mu), m_cut(cut), m_tol(tol), m_nseg(static_cast<unsigned>(m_throttles.size()) / 3u),
      m_nseg_fwd(static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut)), m_nseg_bck(m_nseg - m_nseg_fwd)
{
    // We perform some sanity checks on the user provided inputs
    kep3::leg::_sanity_checks(m_throttles, m_tof, m_max_thrust, m_isp, m_mu, m_cut, m_tol, m_nseg, m_nseg_fwd,
                              m_nseg_bck);

    // Initialize m_tas and m_tas_var
    const heyoka::taylor_adaptive<double> ta_cache = kep3::ta::get_ta_stark(m_tol);
    m_tas = ta_cache;
    const heyoka::taylor_adaptive<double> ta_var_cache = kep3::ta::get_ta_stark_var(m_tol);
    m_tas_var = ta_var_cache;

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
    kep3::leg::_check_tof(tof);
    m_tof = tof;
}
void sims_flanagan_hf::set_rvs(const std::array<std::array<double, 3>, 2> &rv)
{
    std::copy(rv[0].begin(), rv[0].end(), m_rvms.begin());
    std::copy(rv[1].begin(), rv[1].end(), std::next(m_rvms.begin(), 3));
}
void sims_flanagan_hf::set_ms(double mass)
{
    m_rvms[6] = mass;
}
void sims_flanagan_hf::set_throttles(const std::vector<double> &throttles)
{
    kep3::leg::_check_throttles(throttles);
    m_throttles = throttles;
    m_nseg = static_cast<unsigned>(m_throttles.size()) / 3u;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}
void sims_flanagan_hf::set_throttles(const std::vector<double>::const_iterator &it1,
                                     const std::vector<double>::const_iterator &it2)
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

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}
void sims_flanagan_hf::set_rvf(const std::array<std::array<double, 3>, 2> &rv)
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
    kep3::leg::_check_max_thrust(max_thrust);
    m_max_thrust = max_thrust;
}
void sims_flanagan_hf::set_isp(double isp)
{
    kep3::leg::_check_isp(isp);
    m_isp = isp;
}
void sims_flanagan_hf::set_mu(double mu)
{
    kep3::leg::_check_mu(mu);
    m_mu = mu;
}
void sims_flanagan_hf::set_cut(double cut)
{
    kep3::leg::_check_cut(cut);
    m_cut = cut;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
}
void sims_flanagan_hf::set_tol(double tol)
{
    kep3::leg::_check_tol(tol);
    m_tol = tol;
}
void sims_flanagan_hf::set_rvms(const std::array<double, 7> &rvms)
{
    m_rvms = rvms;
}
void sims_flanagan_hf::set_rvmf(const std::array<double, 7> &rvmf)
{
    m_rvmf = rvmf;
}
// void sims_flanagan_hf::set_tas(const heyoka::taylor_adaptive<double> &tas)
// {
//     m_tas = tas;
// }
// void sims_flanagan_hf::set_tas_var(const heyoka::taylor_adaptive<double> &tas_var)
// {
//     m_tas_var = tas_var;
// }

void sims_flanagan_hf::set(const std::array<std::array<double, 3>, 2> &rvs, double ms,
                           const std::vector<double> &throttles,
                           // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                           const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust,
                           double isp, double mu, double cut, double tol)
{
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
    kep3::leg::_sanity_checks(throttles, tof, max_thrust, isp, mu, cut, tol, m_nseg, m_nseg_fwd, m_nseg_bck);

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}

void sims_flanagan_hf::set(const std::array<double, 7> &rvms, const std::vector<double> &throttles,
                           const std::array<double, 7> &rvmf, double tof, double max_thrust, double isp, double mu,
                           double cut, double tol)
{
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
    kep3::leg::_sanity_checks(throttles, tof, max_thrust, isp, mu, cut, tol, m_nseg, m_nseg_fwd, m_nseg_bck);

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
std::array<std::array<double, 3>, 2> sims_flanagan_hf::get_rvs() const
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
std::array<std::array<double, 3>, 2> sims_flanagan_hf::get_rvf() const
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
// LCOV_EXCL_START
const heyoka::taylor_adaptive<double> &sims_flanagan_hf::get_tas() const
{
    return m_tas;
}
const heyoka::taylor_adaptive<double> &sims_flanagan_hf::get_tas_var() const
{
    return m_tas_var;
}
// LCOV_EXCL_END
const std::array<double, 7> &sims_flanagan_hf::get_rvms() const
{
    return m_rvms;
}
const std::array<double, 7> &sims_flanagan_hf::get_rvmf() const
{
    return m_rvmf;
}

// The core routines
std::array<double, 7> sims_flanagan_hf::compute_mismatch_constraints() const
{
    // General settings
    const double prop_seg_duration = (m_tof / m_nseg);

    // Forward pass
    // Initial state
    // Set the Taylor Integration initial conditions
    m_tas.set_time(0.);
    std::copy(m_rvms.begin(), m_rvms.end(), m_tas.get_state_data());

    // Loop through segments in forward pass of Sims-Flanagan transcription
    for (auto i = 0u; i < m_nseg_fwd; ++i) {
        // Assign current thrusts to Taylor adaptive integrator
        std::copy(std::next(m_thrusts.begin(), static_cast<long>(i * 3)),
                  std::next(m_thrusts.begin(), static_cast<long>(3 * (i + 1))), std::next(m_tas.get_pars_data(), 2));
        // ... and integrate
        auto [status, min_h, max_h, nsteps, _1, _2] = m_tas.propagate_until((i + 1) * prop_seg_duration);
        if (status != heyoka::taylor_outcome::time_limit) {
            throw std::domain_error(
                "stark_problem: failure to reach the final time requested during a propagation."); // LCOV_EXCL_LINE
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
    for (auto i = 0u; i < m_nseg_bck; ++i) {
        // Assign current_thrusts to Taylor adaptive integrator
        std::copy(std::next(m_thrusts.begin(), static_cast<long>((m_nseg - (i + 1)) * 3)),
                  std::next(m_thrusts.begin(), static_cast<long>((m_nseg - i) * 3)),
                  std::next(m_tas.get_pars_data(), 2));
        // ... and integrate
        auto [status, min_h, max_h, nsteps, _1, _2] = m_tas.propagate_until(m_tof - (i + 1) * prop_seg_duration);
        if (status != heyoka::taylor_outcome::time_limit) {
            throw std::domain_error(
                "stark_problem: failure to reach the final time requested during a propagation."); // LCOV_EXCL_LINE
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

std::vector<double> sims_flanagan_hf::compute_constraints() const
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
    // Inequality Constraints
    auto ineq_con = compute_throttle_constraints();
    std::copy(ineq_con.begin(), ineq_con.end(), retval.begin() + 7);
    return retval;
}

std::vector<double> sims_flanagan_hf::set_and_compute_constraints(const std::vector<double> &chromosome)
{
    std::array<double, 7> rvms;
    std::copy(chromosome.begin(), chromosome.begin() + 7, rvms.begin());
    std::vector<double> throttles(m_nseg * 3);
    std::copy(chromosome.begin() + 7, chromosome.begin() + 7 + m_nseg * 3, throttles.begin());
    std::array<double, 7> rvmf;
    std::copy(chromosome.begin() + 7 + m_nseg * 3, chromosome.begin() + 7 + m_nseg * 3 + 7, rvmf.begin());
    double time_of_flight = chromosome[(7 + m_nseg * 3 + 7 + 1) - 1];
    // Set relevant quantities before evaluating constraints
    set(rvms, throttles, rvmf, time_of_flight);
    // Evaluate and return constraints
    return compute_constraints();
}

// Return specific two-body 'stark' dynamics state derivative
std::array<double, 7> sims_flanagan_hf::get_state_derivative(const std::array<double, 7> &state,
                                                             const std::array<double, 3> &throttles) const
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
    dstatedt[3] = -get_mu() * std::pow(r2, -3. / 2) * state[0] + thrusts[0] / state[6];
    dstatedt[4] = -get_mu() * std::pow(r2, -3. / 2) * state[1] + thrusts[1] / state[6];
    dstatedt[5] = -get_mu() * std::pow(r2, -3. / 2) * state[2] + thrusts[2] / state[6];
    dstatedt[6] = (u_norm != 0) ? -u_norm / veff : 0; // Conditional for if thrust is zero or not

    return dstatedt;
}

std::tuple<std::vector<std::array<double, 49u>>, std::vector<std::array<double, 21u>>,
           std::vector<std::array<double, 7u>>>
sims_flanagan_hf::compute_all_gradients() const
{
    // Initialise
    std::vector<std::array<double, 7u>> xf_per_seg(m_nseg, {0});
    std::vector<std::array<double, 49u>> dxdx_per_seg(m_nseg, {0});
    std::vector<std::array<double, 21u>> dxdu_per_seg(m_nseg, {0});
    // For ToF gradient
    std::vector<std::array<double, 7u>> x0_per_seg(m_nseg, {0});
    std::vector<std::array<double, 7u>> dxdtof_per_seg(m_nseg, {0});

    // General settings
    const double prop_seg_duration = (m_tof / m_nseg);

    // Forward loop
    // Set the Taylor Integration initial conditions
    m_tas_var.set_time(0.);
    std::copy(m_rvms.begin(), m_rvms.end(), m_tas_var.get_state_data());

    for (auto i = 0u; i < m_nseg_fwd; ++i) {

        // Initialise var conditions
        std::copy(m_vars.begin(), m_vars.end(), m_tas_var.get_state_data() + 7);
        // Assign current thrusts to Taylor adaptive integrator
        std::copy(std::next(m_thrusts.begin(), static_cast<long>(i * 3)),
                  std::next(m_thrusts.begin(), static_cast<long>(3 * (i + 1))),
                  std::next(m_tas_var.get_pars_data(), 2));
        // ... and integrate
        auto [status, min_h, max_h, nsteps, _1, _2] = m_tas_var.propagate_until((i + 1) * prop_seg_duration);
        if (status != heyoka::taylor_outcome::time_limit) {
            throw std::domain_error(
                "stark_problem: failure to reach the final time requested during a propagation."); // LCOV_EXCL_LINE
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

    for (auto i = 0u; i < m_nseg_bck; ++i) {

        // Initialise var conditions
        std::copy(m_vars.begin(), m_vars.end(), m_tas_var.get_state_data() + 7);
        // Assign current thrusts to Taylor adaptive integrator
        std::copy(std::next(m_thrusts.begin(), static_cast<long>((m_nseg - (i + 1)) * 3)),
                  std::next(m_thrusts.begin(), static_cast<long>((m_nseg - i) * 3)),
                  std::next(m_tas_var.get_pars_data(), 2));
        // ... and integrate
        auto [status, min_h, max_h, nsteps, _1, _2] = m_tas_var.propagate_until(m_tof - (i + 1) * prop_seg_duration);
        if (status != heyoka::taylor_outcome::time_limit) {
            throw std::domain_error(
                "stark_problem: failure to reach the final time requested during a propagation."); // LCOV_EXCL_LINE
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

    // Get ToF gradients
    // Initialize initial state matrix
    if (m_nseg_fwd > 0) {
        x0_per_seg[0] = m_rvms;
    }
    for (decltype(m_nseg_fwd) i = 1; i < m_nseg_fwd; ++i) {
        x0_per_seg[i] = xf_per_seg[i - 1];
    }
    if (m_nseg_bck > 0) {
        x0_per_seg[m_nseg - 1] = m_rvmf;
    }
    for (decltype(m_nseg_bck) i = 1; i < m_nseg_bck; ++i) {
        x0_per_seg[(m_nseg - 1) - i] = xf_per_seg[(m_nseg - 1) - (i - 1)];
    }

    for (decltype(dxdtof_per_seg.size()) i = 0; i < dxdtof_per_seg.size(); ++i) {
        std::array<double, 3> current_throttles = {m_throttles[i * 3], m_throttles[i * 3 + 1], m_throttles[i * 3 + 2]};
        dxdtof_per_seg[i] = get_state_derivative(x0_per_seg[i], current_throttles);
    }

    return std::make_tuple(dxdx_per_seg, dxdu_per_seg, dxdtof_per_seg);
}

std::tuple<std::array<double, 49>, std::array<double, 49>, std::vector<double>>
sims_flanagan_hf::get_relevant_gradients(const std::vector<std::array<double, 49u>> &dxdx_per_seg,
                                         const std::vector<std::array<double, 21u>> &dxdu_per_seg,
                                         const std::vector<std::array<double, 7u>> &dxdtof_per_seg) const
{

    auto xt_dxdx_per_seg = xt::adapt(reinterpret_cast<const double *>(dxdx_per_seg.data()), {m_nseg, 49u});
    // Mn_o will contain [Mnf-1, Mnf-1@Mnf-2, Mnf-2@Mnf-3, Mnf-1@M0, Mnf, Mnf@Mnf+1, Mnf@Mnf+2, Mnf@Mn]
    std::vector<xt::xarray<double>> Mn_o(m_nseg, xt::zeros<double>({7u, 7u}));
    // Fwd leg
    xt::xarray<double> final_M;
    xt::xarray<double> current_M;
    if (m_nseg_fwd > 0) {
        Mn_o[0] = xt::reshape_view(xt::view(xt_dxdx_per_seg, m_nseg_fwd - 1, xt::all()), {7, 7});
        for (decltype(m_nseg_fwd) i = 0; i < m_nseg_fwd - 1; ++i) {
            current_M = xt::reshape_view(xt::view(xt_dxdx_per_seg, m_nseg_fwd - 1 - (i + 1), xt::all()), {7, 7});
            if (i == 0) {
                final_M = xt::reshape_view(xt::view(xt_dxdx_per_seg, m_nseg_fwd - 1, xt::all()), {7, 7});
            } else {
                final_M = Mn_o[i];
            }
            Mn_o[i + 1] = xt::linalg::dot(final_M, current_M);
        }
    }
    // Bck leg
    if (m_nseg_bck > 0) {
        Mn_o[m_nseg_fwd] = xt::reshape_view(xt::view(xt_dxdx_per_seg, m_nseg_fwd, xt::all()), {7, 7});
        for (decltype(m_nseg_fwd) i(0); i < m_nseg_bck - 1; ++i) {
            current_M = xt::reshape_view(xt::view(xt_dxdx_per_seg, m_nseg_fwd + (i + 1), xt::all()), {7, 7});
            if (i == 0) {
                final_M = xt::reshape_view(xt::view(xt_dxdx_per_seg, m_nseg_fwd, xt::all()), {7, 7});
            } else {
                final_M = Mn_o[m_nseg_fwd + i];
            }
            Mn_o[m_nseg_fwd + i + 1] = xt::linalg::dot(final_M, current_M);
        }
    }

    // Initial and final displacements
    std::array<double, 49> grad_rvm = {0};
    auto xgrad_rvm = xt::adapt(grad_rvm, {7u, 7u});
    if (m_nseg_fwd > 0) {
        xt::view(xgrad_rvm, xt::all(), xt::all()) = xt::view(Mn_o[m_nseg_fwd - 1], xt::all(), xt::all());
    } else {
        xt::view(xgrad_rvm, xt::all(), xt::all()) = xt::eye(7);
    }

    std::array<double, 49> grad_rvm_bck = {0};
    auto xgrad_rvm_bck = xt::adapt(grad_rvm_bck, {7u, 7u});
    if (m_nseg_bck > 0) {
        xt::view(xgrad_rvm_bck, xt::all(), xt::all())
            = xt::view(Mn_o[m_nseg - 1], xt::all(), xt::all()) * -1; // Multiple by -1 because mass correlation is -1.
    } else {
        xt::view(xgrad_rvm_bck, xt::all(), xt::all()) = xt::eye(7) * -1;
    }

    // Throttle derivatives
    xt::xarray<double> xt_dxdu_per_seg
        = xt::adapt(reinterpret_cast<const double *>(dxdu_per_seg.data()), {m_nseg, 21u});
    std::vector<double> grad_final_throttle(static_cast<size_t>(7) * (m_nseg * 3u), 0.);
    auto xgrad_final_throttle = xt::adapt(grad_final_throttle, {7u, static_cast<unsigned>(m_nseg) * 3u});
    xt::xarray<double> corresponding_M;
    xt::xarray<double> current_U;
    for (decltype(m_nseg_fwd) i(0); i < m_nseg; ++i) {
        current_U = xt::reshape_view(xt::view(xt_dxdu_per_seg, i, xt::all()), {7, 3});
        if (i == m_nseg_fwd - 1) {
            corresponding_M = xt::eye(7);
        } else if (i == m_nseg_fwd) {
            corresponding_M = xt::eye(7) * -1; // Multiple by -1 because mass correlation is -1.
        } else if (i <= m_nseg_fwd - 2 && m_nseg_fwd >= 2) {
            corresponding_M = Mn_o[m_nseg_fwd - 2 - i];
        } else if (i > m_nseg_fwd) {
            corresponding_M = Mn_o[i - 1] * -1; // Multiple by -1 because mass correlation is -1.
        } else {
            throw std::runtime_error(
                "During calculation of the throttle derivatives, the index doesn't correspond to "
                "any leg and therefore cannot find the corresponding gradients."); // LCOV_EXCL_LINE
        }
        xt::view(xgrad_final_throttle, xt::all(), xt::range(3 * i, 3 * (i + 1)))
            = xt::linalg::dot(corresponding_M, current_U);
    }

    // ToF derivatives
    xt::xarray<double> xt_dxdtof_per_seg
        = xt::adapt(reinterpret_cast<const double *>(dxdtof_per_seg.data()), {m_nseg, 7u});
    std::vector<double> grad_final_tof(static_cast<size_t>(7), 0.);
    auto xgrad_final_tof = xt::adapt(grad_final_tof, {7u, 1u});
    for (decltype(m_nseg_fwd) i(0); i < m_nseg; ++i) {
        xt::xarray<double> current_F = xt::reshape_view(xt::view(xt_dxdtof_per_seg, i, xt::all()), {7, 1});
        if ((i <= m_nseg_fwd - 1) && m_nseg_fwd > 0) {
            corresponding_M = Mn_o
                [m_nseg_fwd - 1
                 - i]; // +1 w.r.t. throttle derivatives because dx/dtof is defined at begin of leg rather than end
        } else if ((static_cast<int>(i) > static_cast<int>(m_nseg_fwd) - 1) && m_nseg_bck > 0) {
            corresponding_M = Mn_o[i]; // Idem
        } else {
            throw std::runtime_error(
                "During calculation of the tof derivatives, the index doesn't correspond to "
                "any leg and therefore cannot find the corresponding gradients."); // LCOV_EXCL_LINE
        }
        xgrad_final_tof += xt::linalg::dot(corresponding_M, current_F);
    }
    xgrad_final_tof /= m_nseg;

    // Combine throttle and tof matrices
    std::vector<double> grad_final(static_cast<size_t>(7) * (m_nseg * 3u + 1u), 0.);
    auto xgrad_final = xt::adapt(grad_final, {7u, static_cast<unsigned>(m_nseg) * 3u + 1u});
    xt::view(xgrad_final, xt::all(), xt::range(0, m_nseg * 3)) = xt::view(xgrad_final_throttle, xt::all(), xt::all());
    xt::view(xgrad_final, xt::all(), m_nseg * 3) = xt::view(xgrad_final_tof, xt::all(), 0);

    return {std::move(grad_rvm), std::move(grad_rvm_bck), std::move(grad_final)};
}

std::tuple<std::array<double, 49>, std::array<double, 49>, std::vector<double>>
sims_flanagan_hf::compute_mc_grad() const
{
    // Initialise
    std::vector<std::array<double, 49u>> dxdx_per_seg;
    std::vector<std::array<double, 21u>> dxdu_per_seg;
    std::vector<std::array<double, 7u>> dxdtof_per_seg;
    std::tie(dxdx_per_seg, dxdu_per_seg, dxdtof_per_seg) = compute_all_gradients();

    std::array<double, 49> grad_rvm = {0};
    std::array<double, 49> grad_rvm_bck = {0};
    std::vector<double> grad_final(static_cast<size_t>(7) * (m_nseg * 3u + 1u), 0.);
    std::tie(grad_rvm, grad_rvm_bck, grad_final) = get_relevant_gradients(dxdx_per_seg, dxdu_per_seg, dxdtof_per_seg);

    return {grad_rvm, grad_rvm_bck, std::move(grad_final)};
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

std::vector<std::vector<double>> sims_flanagan_hf::get_state_history(unsigned grid_points_per_segment) const
{
    // Get time grid
    const double prop_seg_duration = (m_tof / m_nseg);
    std::vector<double> leg_time_grid;
    // Initial time
    double timestep = 0.0;
    leg_time_grid.push_back(timestep);

    for (decltype(m_nseg) i = 0; i < grid_points_per_segment * m_nseg - 2; ++i) {
        timestep += (prop_seg_duration / (grid_points_per_segment - 1));
        leg_time_grid.push_back(timestep);
    }
    // leg_time_grid.push_back(m_tof);
    std::vector<double> current_leg_time_grid(grid_points_per_segment);

    // Forward pass
    // Initial state
    // Set the Taylor Integration initial conditions
    m_tas.set_time(0.);
    std::copy(m_rvms.begin(), m_rvms.end(), m_tas.get_state_data());
    std::vector<std::vector<double>> output_per_seg(m_nseg);

    // Loop through segments in forward pass of Sims-Flanagan transcription
    for (decltype(m_nseg_fwd) i = 0u; i < m_nseg_fwd; ++i) {
        // Assign current thrusts to Taylor adaptive integrator
        std::copy(std::next(m_thrusts.begin(), static_cast<long>(i * 3)),
                  std::next(m_thrusts.begin(), static_cast<long>(3 * (i + 1))), std::next(m_tas.get_pars_data(), 2));

        // Current leg time grid
        std::copy(std::next(leg_time_grid.begin(), i * (grid_points_per_segment - 1)),
                  std::next(leg_time_grid.begin(), (i + 1) * (grid_points_per_segment - 1) + 1),
                  current_leg_time_grid.begin());
        m_tas.set_time(current_leg_time_grid.at(0));
        // ... and integrate
        auto [status, min_h, max_h, nsteps, _1, output_states] = m_tas.propagate_grid(current_leg_time_grid);
        if (status != heyoka::taylor_outcome::time_limit) {
            throw std::domain_error(
                "stark_problem: failure to reach the final time requested during a propagation."); // LCOV_EXCL_LINE
        }
        output_per_seg[i] = output_states;
    }

    // Backward pass
    // Final state
    // Set the Taylor Integration final conditions
    m_tas.set_time(m_tof);
    std::copy(m_rvmf.begin(), m_rvmf.end(), m_tas.get_state_data());
    std::vector<double> back_time_grid(grid_points_per_segment);

    // Loop through segments in backward pass of Sims-Flanagan transcription
    for (decltype(m_nseg) i = 0u; i < m_nseg_bck; ++i) {
        // Assign current_thrusts to Taylor adaptive integrator
        std::copy(std::next(m_thrusts.begin(), static_cast<long>((m_nseg - (i + 1)) * 3)),
                  std::next(m_thrusts.begin(), static_cast<long>((m_nseg - i) * 3)),
                  std::next(m_tas.get_pars_data(), 2));

        // Current leg time grid
        std::reverse_copy(leg_time_grid.begin() + (m_nseg - (i + 1)) * (grid_points_per_segment - 1),
                          leg_time_grid.begin() + (m_nseg - i) * (grid_points_per_segment - 1) + 1,
                          back_time_grid.begin());
        m_tas.set_time(back_time_grid.at(0));

        // ... and integrate
        auto [status, min_h, max_h, nsteps, _1, output_states] = m_tas.propagate_grid(back_time_grid);
        if (status != heyoka::taylor_outcome::time_limit) {
            throw std::domain_error(
                "stark_problem: failure to reach the final time requested during a propagation."); // LCOV_EXCL_LINE
        }
        output_per_seg[m_nseg - 1 - i] = output_states;
    }

    return output_per_seg;
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