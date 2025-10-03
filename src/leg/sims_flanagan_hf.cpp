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
#include <heyoka/kw.hpp>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <sys/types.h>
#include <vector>

#include <boost/range/algorithm.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/io/xio.hpp>
#include <xtensor/views/xview.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/epoch.hpp>
#include <kep3/leg/sf_checks.hpp>
#include <kep3/leg/sims_flanagan_hf.hpp>
#include <kep3/linalg.hpp>
#include <kep3/ta/zero_hold_kep.hpp>

#include <heyoka/taylor.hpp>

// Anonymous namespace for functions
namespace
{
// Check that the dynamics is a zero hold one
void _check_zero_hold(const std::pair<const heyoka::taylor_adaptive<double> &, const heyoka::taylor_adaptive<double> &> &tas) {
    const heyoka::taylor_adaptive<double> &ta = tas.first;
    const heyoka::taylor_adaptive<double> &ta_var = tas.second;
    // 1) the ta must have dimension 7
    if (!ta.get_state().size() == 7) {
         throw std::domain_error(
            "The zero-hold Taylor adaptive integrator in the high fidelity sf leg, does not have a dimension of 7.");
    }
    // 2) the ta must have at least 5 parameters mu, veff, T1,T2,T3
    if (ta.get_pars().size() < 5) {
         throw std::domain_error(
            "The zero-hold Taylor adaptive integrator in the high fidelity sf leg, must have at least 5 parameters");
    }
    // 3) the ta is not a variational integrator
    if (ta.is_variational()) {
        throw std::domain_error(
            "The zero-hold Taylor adaptive integrator in the high fidelity sf leg cannot be variational");
    }
    // 4) the ta_var must have dimension 7 + 7 * 7 + 3 * 7
    if (!(ta_var.get_state().size() == 77)) {
         throw std::domain_error(
            "The zero-hold variational Taylor adaptive integrator in the high fidelity sf leg, does not have a dimension of 28.");
    }
    // 2) the ta_var must have the same parameters as ts
    if (ta.get_pars().size() != ta_var.get_pars().size()) {
         throw std::domain_error(
            "The zero-hold Taylor adaptive integrator and its varitional counterpart seem to have different number of parameters?");
    }
    // 3) the ta_var is a variational integrator
    if (!ta_var.is_variational()) {
        throw std::domain_error(
            "The zero-hold variartional Taylor adaptive integrator in the high fidelity sf must be variational");
    }
    return;
}
// Utilty
std::array<double, 7> make_rvm(const std::array<std::array<double, 3>, 2> &posvel, double m)
{
    return {posvel[0][0], posvel[0][1], posvel[0][2], posvel[1][0], posvel[1][1], posvel[1][2], m};
}

// Returning the dynamics excluding mass
heyoka::cfunc<double> dynamic_cfunc_factory(const heyoka::taylor_adaptive<double> &ta)
{
    auto dyn = ta.get_sys();
    // Collect variables and RHS expressions
    std::vector<heyoka::expression> variable;
    std::vector<heyoka::expression> rhs;
    variable.reserve(dyn.size());
    rhs.reserve(dyn.size());
    for (auto &row : dyn) {
        variable.push_back(row.first);
        rhs.push_back(row.second);
    }
    // Build compiled function, include parameters in var list
    return heyoka::cfunc<double>(
        {rhs[0], rhs[1], rhs[2], rhs[3], rhs[4], rhs[5], rhs[6]},                                     // outputs
        {variable[0], variable[1], variable[2], variable[3], variable[4], variable[5], variable[6]}); // inputs (state)
}

} // namespace
namespace kep3::leg
{

// Default constructor
sims_flanagan_hf::sims_flanagan_hf()
{
    // We perform some sanity checks on the user provided inputs
    kep3::leg::_sanity_checks(m_throttles, m_tof, m_max_thrust, m_veff, m_mu, m_cut, m_tol, m_nseg, m_nseg_fwd,
                              m_nseg_bck);

    // Initialize m_ta and m_ta_var
    const heyoka::taylor_adaptive<double> ta_cache = kep3::ta::get_ta_zero_hold_kep(m_tol);
    m_ta = ta_cache;
    const heyoka::taylor_adaptive<double> ta_var_cache = kep3::ta::get_ta_zero_hold_kep_var(m_tol);
    m_ta_var = ta_var_cache;
    m_cf_dyn = dynamic_cfunc_factory(m_ta);

    // We set mu and veff for the non variational
    *m_ta.get_pars_data() = m_mu;
    *(m_ta.get_pars_data() + 1) = m_veff;

    // ... and variational version of the integrator
    *(m_ta_var.get_pars_data()) = m_mu;
    *(m_ta_var.get_pars_data() + 1) = m_veff;
    // We copy the initial conditions for the variational equations
    std::copy(m_ta_var.get_state().begin() + 7, m_ta_var.get_state().end(), m_vars.begin());

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}

// Main constructor
sims_flanagan_hf::sims_flanagan_hf(
    const std::array<double, 7> &rvms, const std::vector<double> &throttles,
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    const std::array<double, 7> &rvmf, double tof, double max_thrust, double veff, double mu, double cut, double tol,
    std::optional<std::pair<const heyoka::taylor_adaptive<double> &, const heyoka::taylor_adaptive<double> &>> tas)
    : m_rvms(rvms), m_throttles(throttles), m_rvmf(rvmf), m_tof(tof), m_max_thrust(max_thrust), m_veff(veff), m_mu(mu),
      m_cut(cut), m_tol(tol), m_nseg(static_cast<unsigned>(m_throttles.size()) / 3u),
      m_nseg_fwd(static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut)), m_nseg_bck(m_nseg - m_nseg_fwd)
{
    // We perform some sanity checks on the user provided inputs
    kep3::leg::_sanity_checks(m_throttles, m_tof, m_max_thrust, m_veff, m_mu, m_cut, m_tol, m_nseg, m_nseg_fwd,
                              m_nseg_bck);

    // Initialize m_ta and m_ta_var
    if (tas) {
        _check_zero_hold(tas.value());
        m_ta = tas.value().first;
        m_ta_var = tas.value().second;
    } else {
        const heyoka::taylor_adaptive<double> ta_cache = kep3::ta::get_ta_zero_hold_kep(m_tol);
        m_ta = ta_cache;
        const heyoka::taylor_adaptive<double> ta_var_cache = kep3::ta::get_ta_zero_hold_kep_var(m_tol);
        m_ta_var = ta_var_cache;
    }
    // Build the cfunc for the dynamics
    m_cf_dyn = dynamic_cfunc_factory(m_ta);

    // We set mu and veff for the non variational
    *m_ta.get_pars_data() = m_mu;
    *(m_ta.get_pars_data() + 1) = m_veff;

    // ... and variational version of the integrator
    *(m_ta_var.get_pars_data()) = m_mu;
    *(m_ta_var.get_pars_data() + 1) = m_veff;
    // We copy the initial conditions for the variational equations
    std::copy(m_ta_var.get_state().begin() + 7, m_ta_var.get_state().end(), m_vars.begin());

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}

// Convenience constructor from posvel, m
sims_flanagan_hf::sims_flanagan_hf(
    const std::array<std::array<double, 3>, 2> &rvs, double ms, const std::vector<double> &throttles,
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust, double veff, double mu,
    double cut, double tol,
    std::optional<std::pair<const heyoka::taylor_adaptive<double> &, const heyoka::taylor_adaptive<double> &>> tas)
    : sims_flanagan_hf(make_rvm(rvs, ms), throttles, make_rvm(rvf, mf), tof, max_thrust, veff, mu, cut, tol, tas)
{
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
    set_throttles(throttles.begin(), throttles.end());
}

void sims_flanagan_hf::set_throttles(const std::vector<double>::const_iterator &it1,
                                     const std::vector<double>::const_iterator &it2)
{
    if (((std::distance(it1, it2) % 3) != 0) || std::distance(it1, it2) <= 0) {
        throw std::logic_error(
            "The throttles of a high-fidelity sims-flanagan leg are being set with invalid iterators.");
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
void sims_flanagan_hf::set_veff(double veff)
{
    kep3::leg::_check_veff(veff);
    m_veff = veff;
    *(m_ta.get_pars_data() + 1l) = veff;
    *(m_ta_var.get_pars_data() + 1l) = veff;
}
void sims_flanagan_hf::set_mu(double mu)
{
    kep3::leg::_check_mu(mu);
    m_mu = mu;
    *m_ta.get_pars_data() = mu;
    *m_ta_var.get_pars_data() = mu;
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

void sims_flanagan_hf::set(const std::array<std::array<double, 3>, 2> &rvs, double ms,
                           const std::vector<double> &throttles,
                           // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                           const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust,
                           double veff, double mu, double cut, double tol)
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
    m_veff = veff;
    m_mu = mu;
    m_cut = cut;
    m_tol = tol;
    m_nseg = static_cast<unsigned>(m_throttles.size()) / 3u;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
    kep3::leg::_sanity_checks(throttles, tof, max_thrust, veff, mu, cut, tol, m_nseg, m_nseg_fwd, m_nseg_bck);

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}

void sims_flanagan_hf::set(const std::array<double, 7> &rvms, const std::vector<double> &throttles,
                           const std::array<double, 7> &rvmf, double tof, double max_thrust, double veff, double mu,
                           double cut, double tol)
{
    set_rvms(rvms);
    m_throttles = throttles;
    set_rvmf(rvmf);
    m_tof = tof;
    m_max_thrust = max_thrust;
    m_veff = veff;
    m_mu = mu;
    m_cut = cut;
    m_tol = tol;
    m_nseg = static_cast<unsigned>(m_throttles.size()) / 3u;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
    kep3::leg::_sanity_checks(throttles, tof, max_thrust, veff, mu, cut, tol, m_nseg, m_nseg_fwd, m_nseg_bck);

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
    _sanity_checks(throttles, m_tof, m_max_thrust, m_veff, m_mu, m_cut, m_tol, m_nseg, m_nseg_fwd, m_nseg_bck);

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
double sims_flanagan_hf::get_veff() const
{
    return m_veff;
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
    return m_ta;
}
const heyoka::taylor_adaptive<double> &sims_flanagan_hf::get_tas_var() const
{
    return m_ta_var;
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
    const double mass_thresh = 1e-12 * (*(m_rvmf.begin() + 6l));

    // Forward pass
    // Initial state
    // Set the Taylor Integration initial conditions
    m_ta.set_time(0.);
    std::copy(m_rvms.begin(), m_rvms.end(), m_ta.get_state_data());

    // Loop through segments in forward pass of Sims-Flanagan transcription
    for (auto i = 0u; i < m_nseg_fwd; ++i) {
        // Assign current thrusts to Taylor adaptive integrator
        std::copy(m_thrusts.begin() + i * 3l, m_thrusts.begin() + 3 * (i + 1l), m_ta.get_pars_data() + 2l);

        if (!std::isfinite(prop_seg_duration)) {
            // fmt::print("Non-finitite propagation duration requested in forward pass\n");
            break;
        } else {
            // ... and integrate
            double norm_thrusts = std::sqrt(std::inner_product(
                m_thrusts.begin() + i * 3l, m_thrusts.begin() + 3 * (i + 1l), m_thrusts.begin() + i * 3l, 0.0));
            double mass_est = m_ta.get_state()[6] - norm_thrusts * prop_seg_duration / (m_veff);
            if (mass_est < mass_thresh) {
            } else {
                auto [status, min_h, max_h, nsteps, _1, _2] = m_ta.propagate_until((i + 1) * prop_seg_duration);
            }
        }
    }

    // Reset veff
    *(m_ta.get_pars_data() + 1l) = get_veff();

    // Set fwd final state
    std::vector<double> rvm_fwd_final = m_ta.get_state();

    // Backward pass
    // Final state
    // Set the Taylor Integration final conditions
    m_ta.set_time(m_tof);
    std::copy(m_rvmf.begin(), m_rvmf.end(), m_ta.get_state_data());

    // Loop through segments in backward pass of Sims-Flanagan transcription
    for (auto i = 0u; i < m_nseg_bck; ++i) {
        // Assign current_thrusts to Taylor adaptive integrator
        std::copy(m_thrusts.begin() + (m_nseg - (i + 1)) * 3l, m_thrusts.begin() + 3l * (m_nseg - i),
                  m_ta.get_pars_data() + 2l);
        if (!std::isfinite(prop_seg_duration)) {
            // fmt::print("Non-finitite propagation duration requested in backward pass\n");
            break;
        } else {
            // ... and integrate
            auto [status, min_h, max_h, nsteps, _1, _2] = m_ta.propagate_until(m_tof - (i + 1) * prop_seg_duration);
            if (status != heyoka::taylor_outcome::time_limit) { // LCOV_EXCL_START
                // fmt::print("mismatch bck: {}\n", status);
                break;
            } // LCOV_EXCL_STOP
        }
    }

    // Set fwd final state
    std::vector<double> rvm_bck_final = m_ta.get_state();

    if (!std::all_of(rvm_fwd_final.begin(), rvm_fwd_final.end(), [](double x) { return std::isfinite(x); })) {
        fmt::print("rvm_fwd_final {}\n", rvm_fwd_final);
    }

    if (!std::all_of(rvm_bck_final.begin(), rvm_bck_final.end(), [](double x) { return std::isfinite(x); })) {
        fmt::print("rvm_bck_final {}\n", rvm_bck_final);
    }

    return {rvm_fwd_final[0] - m_ta.get_state()[0], rvm_fwd_final[1] - m_ta.get_state()[1],
            rvm_fwd_final[2] - m_ta.get_state()[2], rvm_fwd_final[3] - m_ta.get_state()[3],
            rvm_fwd_final[4] - m_ta.get_state()[4], rvm_fwd_final[5] - m_ta.get_state()[5],
            rvm_fwd_final[6] - m_ta.get_state()[6]};
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

// Return specific two-body 'zero_hold_kep' dynamics state derivative
std::array<double, 7> sims_flanagan_hf::get_state_derivative(const std::array<double, 7> &state,
                                                             const std::array<double, 3> &throttles) const
{
    std::array<double, 3> thrusts{};
    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    std::transform(throttles.begin(), throttles.end(), thrusts.begin(), throttle_to_thrust);
    auto outputs = std::array<double, 7>{};
    m_cf_dyn(outputs, state, heyoka::kw::pars = {m_mu, m_veff, thrusts[0], thrusts[1], thrusts[2]});
    return outputs;
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
    const double mass_thresh = 1e-4 * (*(m_rvmf.begin() + 6l));

    // Forward loop
    // Set the Taylor Integration initial conditions
    m_ta_var.set_time(0.);
    std::copy(m_rvms.begin(), m_rvms.end(), m_ta_var.get_state_data());

    for (auto i = 0u; i < m_nseg_fwd; ++i) {

        // Initialise var conditions
        std::copy(m_vars.begin(), m_vars.end(), m_ta_var.get_state_data() + 7);
        // Assign current thrusts to Taylor adaptive integrator
        std::copy(m_thrusts.begin() + i * 3l, m_thrusts.begin() + 3 * (i + 1l), m_ta_var.get_pars_data() + 2l);

        // LCOV_EXCL_START
        if (!std::isfinite(prop_seg_duration)) {
            fmt::print("Non-finite propagation duration requested in forwards step\n");
            break;
        } else { // LCOV_EXCL_STOP
            // ... and integrate
            double norm_thrusts = std::sqrt(std::inner_product(
                m_thrusts.begin() + i * 3l, m_thrusts.begin() + 3 * (i + 1l), m_thrusts.begin() + i * 3l, 0.0));
            double final_mass = m_ta_var.get_state()[6] - norm_thrusts * prop_seg_duration / (m_veff);
            // Perform the integration only if the final mass is above a certain threshold
            if (final_mass > mass_thresh) {
                auto [status, min_h, max_h, nsteps, _1, _2] = m_ta_var.propagate_until((i + 1) * prop_seg_duration);
            }
        }

        // Save the variational state variables to respective arrays
        std::copy(m_ta_var.get_state().begin(), m_ta_var.get_state().begin() + 7, xf_per_seg[i].begin());
        for (auto j = 0; j < 7; ++j) {
            std::copy(std::next(m_ta_var.get_state().begin(), 7 + 10l * j),
                      std::next(m_ta_var.get_state().begin(), 7 + 10l * j + 7),
                      std::next(dxdx_per_seg[i].begin(), 7l * j));
            std::copy(m_ta_var.get_state().begin() + 14l + 10l * j, m_ta_var.get_state().begin() + 14l + 10l * j + 3l,
                      dxdu_per_seg[i].begin() + 3l * j);
        }
    }

    // Backward loop
    // Set the Taylor Integration initial conditions
    m_ta_var.set_time(m_tof);
    std::copy(m_rvmf.begin(), m_rvmf.end(), m_ta_var.get_state_data());

    for (auto i = 0u; i < m_nseg_bck; ++i) {

        // Initialise var conditions
        std::copy(m_vars.begin(), m_vars.end(), m_ta_var.get_state_data() + 7);
        // Assign current thrusts to Taylor adaptive integrator
        std::copy(m_thrusts.begin() + (m_nseg - (i + 1)) * 3l, m_thrusts.begin() + 3l * (m_nseg - i),
                  m_ta_var.get_pars_data() + 2l);
        // LCOV_EXCL_START
        if (!std::isfinite(prop_seg_duration)) {
            fmt::print("Non-finite propagation duration requested in forwards step\n");
            break;
        } else { // LCOV_EXCL_STOP
            auto [status, min_h, max_h, nsteps, _1, _2] = m_ta_var.propagate_until(m_tof - (i + 1) * prop_seg_duration);
        }
        // Save the variational state variables to respective arrays
        std::copy(m_ta_var.get_state().begin(), m_ta_var.get_state().begin() + 7, xf_per_seg[m_nseg - (i + 1)].begin());
        for (auto j = 0; j < 7; ++j) {
            std::copy(m_ta_var.get_state().begin() + 7 + 10l * j, m_ta_var.get_state().begin() + 7 + 10l * j + 7l,
                      dxdx_per_seg[m_nseg - (i + 1)].begin() + 7l * j);
            std::copy(m_ta_var.get_state().begin() + 14 + 10l * j, m_ta_var.get_state().begin() + 14 + 10l * j + 3l,
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

    // We must multiple the derivatives wrt thrust into throttles
    for (auto &item : dxdu_per_seg) {
        std::transform(item.begin(), item.end(), item.begin(), [this](double x) { return x * m_max_thrust; });
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
    std::vector<double> grad_final_throttle(static_cast<size_t>(7) * (m_nseg * 3l), 0.);
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

    return {grad_rvm, grad_rvm_bck, std::move(grad_final)};
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
    std::vector<double> grad_utof(static_cast<size_t>(7) * (m_nseg * 3u + 1u), 0.);
    std::tie(grad_rvm, grad_rvm_bck, grad_utof) = get_relevant_gradients(dxdx_per_seg, dxdu_per_seg, dxdtof_per_seg);

    return {grad_rvm, grad_rvm_bck, std::move(grad_utof)};
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
    m_ta.set_time(0.);
    std::copy(m_rvms.begin(), m_rvms.end(), m_ta.get_state_data());
    std::vector<std::vector<double>> output_per_seg(m_nseg);

    // Loop through segments in forward pass of Sims-Flanagan transcription
    for (decltype(m_nseg_fwd) i = 0u; i < m_nseg_fwd; ++i) {
        // Assign current thrusts to Taylor adaptive integrator
        std::copy(m_thrusts.begin() + i * 3l, m_thrusts.begin() + 3 * (i + 1l), m_ta.get_pars_data() + 2l);

        // Current leg time grid
        std::copy(leg_time_grid.begin() + i * (grid_points_per_segment - 1l),
                  leg_time_grid.begin() + (i + 1l) * (grid_points_per_segment - 1l) + 1l,
                  current_leg_time_grid.begin());
        m_ta.set_time(current_leg_time_grid.at(0));
        // ... and integrate
        auto [status, min_h, max_h, nsteps, _1, output_states] = m_ta.propagate_grid(current_leg_time_grid);
        if (status != heyoka::taylor_outcome::time_limit) { // LCOV_EXCL_START
            fmt::print("thrust: [{}, {}, {}]\n", *(m_ta.get_pars_data() + 2l), *(m_ta.get_pars_data() + 3l),
                       *(m_ta.get_pars_data() + 4l));
            fmt::print(": {}\n", status);
            fmt::print("state: {}\n", m_ta.get_state());
            fmt::print("reached time: {}\n", m_ta.get_time());
            fmt::print("requested time: {}\n", (i + 1) * prop_seg_duration);
            throw std::domain_error(
                "zero_hold_kep_problem: failure to reach the final time requested during a propagation.");
        } // LCOV_EXCL_STOP
        output_per_seg[i] = output_states;
    }

    // Backward pass
    // Final state
    // Set the Taylor Integration final conditions
    m_ta.set_time(m_tof);
    std::copy(m_rvmf.begin(), m_rvmf.end(), m_ta.get_state_data());
    std::vector<double> back_time_grid(grid_points_per_segment);

    // Loop through segments in backward pass of Sims-Flanagan transcription
    for (decltype(m_nseg) i = 0u; i < m_nseg_bck; ++i) {
        // Assign current_thrusts to Taylor adaptive integrator
        std::copy(m_thrusts.begin() + (m_nseg - (i + 1)) * 3l, m_thrusts.begin() + 3l * (m_nseg - i),
                  m_ta.get_pars_data() + 2l);
        // Current leg time grid
        std::reverse_copy(leg_time_grid.begin() + (m_nseg - (i + 1l)) * (grid_points_per_segment - 1l),
                          leg_time_grid.begin() + (m_nseg - i) * (grid_points_per_segment - 1l) + 1l,
                          back_time_grid.begin());
        m_ta.set_time(back_time_grid.at(0));

        // ... and integrate
        auto [status, min_h, max_h, nsteps, _1, output_states] = m_ta.propagate_grid(back_time_grid);
        if (status != heyoka::taylor_outcome::time_limit) { // LCOV_EXCL_START
            fmt::print("thrust: [{}, {}, {}]\n", *(m_ta.get_pars_data() + 2l), *(m_ta.get_pars_data() + 3l),
                       *(m_ta.get_pars_data() + 4l));
            fmt::print(": {}\n", status);
            fmt::print("state: {}\n", m_ta.get_state());
            fmt::print("reached time: {}\n", m_ta.get_time());
            fmt::print("requested time: {}\n", (i + 1) * prop_seg_duration);
            throw std::domain_error(
                "zero_hold_kep_problem: failure to reach the final time requested during a propagation.");
        } // LCOV_EXCL_STOP
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
    s << fmt::format("Specific impulse: {}\n\n", sf.get_veff());
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