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
#include <kep3/leg/sims_flanagan_hf_alpha.hpp>
#include <kep3/linalg.hpp>
#include <kep3/ta/stark.hpp>

#include <heyoka/taylor.hpp>

namespace kep3::leg
{

// Constructors

sims_flanagan_hf_alpha::sims_flanagan_hf_alpha()
{
    // We perform some sanity checks on the user provided inputs
    kep3::leg::_sanity_checks_alpha(m_throttles, m_talphas, m_tof, m_max_thrust, m_isp, m_mu, m_cut, m_tol, m_nseg, m_nseg_fwd,
                              m_nseg_bck);

    // Initialize m_tas and m_tas_var
    const heyoka::taylor_adaptive<double> ta_cache = kep3::ta::get_ta_stark(m_tol);
    m_tas = ta_cache;
    const heyoka::taylor_adaptive<double> ta_var_cache = kep3::ta::get_ta_stark_var(m_tol);
    m_tas_var = ta_var_cache;

    // We set mu and veff for the non-variational
    *m_tas.get_pars_data() = m_mu;
    *(m_tas.get_pars_data() + 1) = m_isp * kep3::G0;

    // ... and variational version of the integrator
    *(m_tas_var.get_pars_data()) = m_mu;
    *(m_tas_var.get_pars_data() + 1) = m_isp * kep3::G0;

    // We copy the initial conditions for the variational equations
    std::copy(m_tas_var.get_state().begin() + 7, m_tas_var.get_state().end(), m_vars.begin());

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is the same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}

sims_flanagan_hf_alpha::sims_flanagan_hf_alpha(const std::array<std::array<double, 3>, 2> &rvs, double ms,
                                    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                                   const std::vector<double> &throttles, const std::vector<double> &talphas,
                                    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                                   const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof,
                                   
                                   double max_thrust, double isp, double mu, double cut, double tol)
    : m_throttles(throttles), m_talphas(talphas), m_tof(tof), m_max_thrust(max_thrust), m_isp(isp), m_mu(mu), m_cut(cut), m_tol(tol),
      m_nseg(static_cast<unsigned>(m_throttles.size()) / 3u),
      m_nseg_fwd(static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut)), m_nseg_bck(m_nseg - m_nseg_fwd)
{
    // We perform some sanity checks on the user provided inputs
    kep3::leg::_sanity_checks_alpha(m_throttles, m_talphas, m_tof, m_max_thrust, m_isp, m_mu, m_cut, m_tol, m_nseg, m_nseg_fwd,
                              m_nseg_bck);

    // Initialize m_tas and m_tas_var
    const heyoka::taylor_adaptive<double> ta_cache = kep3::ta::get_ta_stark(m_tol);
    m_tas = ta_cache;
    const heyoka::taylor_adaptive<double> ta_var_cache = kep3::ta::get_ta_stark_var(m_tol);
    m_tas_var = ta_var_cache;

    // We set mu and veff for the non-variational
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

sims_flanagan_hf_alpha::sims_flanagan_hf_alpha(const std::array<double, 7> &rvms, 
                                   // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                                   const std::vector<double> &throttles, const std::vector<double> &talphas, 
                                   const std::array<double, 7> &rvmf, 
                                   // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                                   double tof, double max_thrust, double isp, double mu, double cut, double tol)
    : m_rvms(rvms), m_throttles(throttles), m_talphas(talphas), m_rvmf(rvmf), m_tof(tof), m_max_thrust(max_thrust), 
      m_isp(isp), m_mu(mu), m_cut(cut), m_tol(tol), m_nseg(static_cast<unsigned>(m_throttles.size()) / 3u),
      m_nseg_fwd(static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut)), m_nseg_bck(m_nseg - m_nseg_fwd)
{
    // We perform some sanity checks on the user provided inputs
    kep3::leg::_sanity_checks_alpha(m_throttles, m_talphas, m_tof, m_max_thrust, m_isp, m_mu, m_cut, m_tol, m_nseg, m_nseg_fwd,
                              m_nseg_bck);

    // Initialize m_tas and m_tas_var
    const heyoka::taylor_adaptive<double> ta_cache = kep3::ta::get_ta_stark(m_tol);
    m_tas = ta_cache;
    const heyoka::taylor_adaptive<double> ta_var_cache = kep3::ta::get_ta_stark_var(m_tol);
    m_tas_var = ta_var_cache;

    // We set mu and veff for the non-variational
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
void sims_flanagan_hf_alpha::set_tof(double tof)
{
    kep3::leg::_check_tof(tof);
    m_tof = tof;
}
void sims_flanagan_hf_alpha::set_rvs(const std::array<std::array<double, 3>, 2> &rv)
{
    std::copy(rv[0].begin(), rv[0].end(), m_rvms.begin());
    std::copy(rv[1].begin(), rv[1].end(), std::next(m_rvms.begin(), 3));
}
void sims_flanagan_hf_alpha::set_ms(double mass)
{
    m_rvms[6] = mass;
}
void sims_flanagan_hf_alpha::set_throttles(const std::vector<double> &throttles)
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

    *(m_tas.get_pars_data()+2l) = m_thrusts[0];
    *(m_tas.get_pars_data()+3l) = m_thrusts[1];
    *(m_tas.get_pars_data()+4l) = m_thrusts[2];

    *(m_tas_var.get_pars_data()+2l) = m_thrusts[0];
    *(m_tas_var.get_pars_data()+3l) = m_thrusts[1];
    *(m_tas_var.get_pars_data()+4l) = m_thrusts[2];

}
void sims_flanagan_hf_alpha::set_throttles(const std::vector<double>::const_iterator &it1,
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

    *(m_tas.get_pars_data()+2l) = m_thrusts[0];
    *(m_tas.get_pars_data()+3l) = m_thrusts[1];
    *(m_tas.get_pars_data()+4l) = m_thrusts[2];

    *(m_tas_var.get_pars_data()+2l) = m_thrusts[0];
    *(m_tas_var.get_pars_data()+3l) = m_thrusts[1];
    *(m_tas_var.get_pars_data()+4l) = m_thrusts[2];

}
void sims_flanagan_hf_alpha::set_rvf(const std::array<std::array<double, 3>, 2> &rv)
{
    std::copy(rv[0].begin(), rv[0].end(), m_rvmf.begin());
    std::copy(rv[1].begin(), rv[1].end(), std::next(m_rvmf.begin(), 3));
}
void sims_flanagan_hf_alpha::set_mf(double mass)
{
    m_rvmf[6] = mass;
}
void sims_flanagan_hf_alpha::set_max_thrust(double max_thrust)
{
    kep3::leg::_check_max_thrust(max_thrust);
    m_max_thrust = max_thrust;
}
void sims_flanagan_hf_alpha::set_isp(double isp)
{
    kep3::leg::_check_isp(isp);
    m_isp = isp;
    *(m_tas.get_pars_data()+1l) = isp * kep3::G0;
    *(m_tas_var.get_pars_data()+1l) = isp * kep3::G0;
}
void sims_flanagan_hf_alpha::set_mu(double mu)
{
    kep3::leg::_check_mu(mu);
    m_mu = mu;
    *m_tas.get_pars_data() = mu;
    *m_tas_var.get_pars_data() = mu;
}
void sims_flanagan_hf_alpha::set_cut(double cut)
{
    kep3::leg::_check_cut(cut);
    m_cut = cut;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
}
void sims_flanagan_hf_alpha::set_tol(double tol)
{
    kep3::leg::_check_tol(tol);
    m_tol = tol;
}
void sims_flanagan_hf_alpha::set_rvms(const std::array<double, 7> &rvms)
{
    m_rvms = rvms;
}
void sims_flanagan_hf_alpha::set_rvmf(const std::array<double, 7> &rvmf)
{
    m_rvmf = rvmf;
}
// Setter for m_talphas
void sims_flanagan_hf_alpha::set_talphas(const std::vector<double> &talphas)
{
    m_talphas = talphas;
}
// void sims_flanagan_hf_alpha::set_talphas(const std::vector<double>::const_iterator &it1,
//     const std::vector<double>::const_iterator &it2)
// {
//     if ( std::distance(it1, it2) <= 0) {
//         throw std::logic_error("The talphas of a sims_flanagan_hf_alpha leg are being set with invalid iterators.");
//     }
//     m_talphas.resize(static_cast<size_t>(std::distance(it1, it2)));
//     std::copy(it1, it2, m_talphas.begin());
// }

// void sims_flanagan_hf_alpha::set_tas(const heyoka::taylor_adaptive<double> &tas)
// {
//     m_tas = tas;
// }
// void sims_flanagan_hf_alpha::set_tas_var(const heyoka::taylor_adaptive<double> &tas_var)
// {
//     m_tas_var = tas_var;
// }
// 

void sims_flanagan_hf_alpha::set(const std::array<std::array<double, 3>, 2> &rvs, double ms,
    const std::vector<double> &throttles, const std::vector<double> &talphas,
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
    m_talphas = talphas;  // Set the talphas vector
    m_tof = tof;
    m_max_thrust = max_thrust;
    m_isp = isp;
    m_mu = mu;
    m_cut = cut;
    m_tol = tol;
    m_nseg = static_cast<unsigned>(m_throttles.size()) / 3u;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
    kep3::leg::_sanity_checks_alpha(throttles, talphas, tof, max_thrust, isp, mu, cut, tol, m_nseg, m_nseg_fwd, m_nseg_bck);

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}

void sims_flanagan_hf_alpha::set(const std::array<double, 7> &rvms, const std::vector<double> &throttles,
    const std::vector<double> &talphas, const std::array<double, 7> &rvmf, double tof, double max_thrust,
    double isp, double mu, double cut, double tol)
{
    set_rvms(rvms);
    m_throttles = throttles;
    m_talphas = talphas;  // Set the talphas vector
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
    kep3::leg::_sanity_checks_alpha(throttles, talphas, tof, max_thrust, isp, mu, cut, tol, m_nseg, m_nseg_fwd, m_nseg_bck);

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}

void sims_flanagan_hf_alpha::set(const std::array<double, 7> &rvms, const std::vector<double> &throttles,
    const std::vector<double> &talphas, const std::array<double, 7> &rvmf, double time_of_flight)
{
    set_rvms(rvms);
    m_throttles = throttles;
    m_talphas = talphas;  // Set the talphas vector
    set_rvmf(rvmf);
    m_tof = time_of_flight;
    m_nseg = static_cast<unsigned>(m_throttles.size()) / 3u;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
    _sanity_checks_alpha(throttles, talphas, m_tof, m_max_thrust, m_isp, m_mu, m_cut, m_tol, m_nseg, m_nseg_fwd, m_nseg_bck);

    // Convert throttles to current_thrusts.
    auto throttle_to_thrust = [this](double throttle) { return throttle * get_max_thrust(); };
    m_thrusts.resize(m_throttles.size()); // Ensure that std::vector m_thrusts is same size as m_throttles
    std::transform(m_throttles.begin(), m_throttles.end(), m_thrusts.begin(), throttle_to_thrust);
}


// Getters
double sims_flanagan_hf_alpha::get_tof() const
{
    return m_tof;
}
std::array<std::array<double, 3>, 2> sims_flanagan_hf_alpha::get_rvs() const
{
    std::array<std::array<double, 3>, 2> rvs{};
    std::copy(m_rvms.begin(), std::next(m_rvms.begin(), 3), rvs[0].begin());
    std::copy(std::next(m_rvms.begin(), 3), std::next(m_rvms.begin(), 6), rvs[1].begin());
    return rvs;
}
double sims_flanagan_hf_alpha::get_ms() const
{
    return m_rvms[6];
}
const std::vector<double> &sims_flanagan_hf_alpha::get_throttles() const
{
    return m_throttles;
}
std::array<std::array<double, 3>, 2> sims_flanagan_hf_alpha::get_rvf() const
{
    std::array<std::array<double, 3>, 2> rvf{{{0., 0., 0.}, {0., 0., 0.}}};
    std::copy(m_rvmf.begin(), std::next(m_rvmf.begin(), 3), rvf[0].begin());
    std::copy(std::next(m_rvmf.begin(), 3), std::next(m_rvmf.begin(), 6), rvf[1].begin());
    return rvf;
}
double sims_flanagan_hf_alpha::get_mf() const
{
    return m_rvmf[6];
}
double sims_flanagan_hf_alpha::get_max_thrust() const
{
    return m_max_thrust;
}
double sims_flanagan_hf_alpha::get_isp() const
{
    return m_isp;
}
double sims_flanagan_hf_alpha::get_mu() const
{
    return m_mu;
}
double sims_flanagan_hf_alpha::get_cut() const
{
    return m_cut;
}
double sims_flanagan_hf_alpha::get_tol() const
{
    return m_tol;
}
unsigned sims_flanagan_hf_alpha::get_nseg() const
{
    return m_nseg;
}
unsigned sims_flanagan_hf_alpha::get_nseg_fwd() const
{
    return m_nseg_fwd;
}
unsigned sims_flanagan_hf_alpha::get_nseg_bck() const
{
    return m_nseg_bck;
}
// Getter for m_talphas
const std::vector<double> &sims_flanagan_hf_alpha::get_talphas() const
{
    return m_talphas;
}
// LCOV_EXCL_START
const heyoka::taylor_adaptive<double> &sims_flanagan_hf_alpha::get_tas() const
{
    return m_tas;
}
const heyoka::taylor_adaptive<double> &sims_flanagan_hf_alpha::get_tas_var() const
{
    return m_tas_var;
}
// LCOV_EXCL_END
const std::array<double, 7> &sims_flanagan_hf_alpha::get_rvms() const
{
    return m_rvms;
}
const std::array<double, 7> &sims_flanagan_hf_alpha::get_rvmf() const
{
    return m_rvmf;
}

// The core routines
std::array<double, 7> sims_flanagan_hf_alpha::compute_mismatch_constraints() const
{
    // General settings
    const double mass_thresh = 1e-12 * (*(m_rvmf.begin()+6l));


    // Forward pass
    // Initial state
    // Set the Taylor Integration initial conditions
    double prop_seg_temp = 0.0;

    m_tas.set_time(prop_seg_temp);
    std::copy(m_rvms.begin(), m_rvms.end(), m_tas.get_state_data());

    
    // Loop through segments in forward pass of Sims-Flanagan transcription
    for (auto i = 0u; i < m_nseg_fwd; ++i) {
        // Assign current thrusts to Taylor adaptive integrator
        std::copy(m_thrusts.begin() + i * 3l, m_thrusts.begin() + 3 * (i + 1l), m_tas.get_pars_data() + 2l);

        // ... and integrate
        prop_seg_temp += m_talphas[i];

        // ... and integrate
        double norm_thrusts = std::sqrt(std::inner_product(m_thrusts.begin() + i * 3l, m_thrusts.begin() + 3 * (i + 1l), m_thrusts.begin() + i * 3l, 0.0));
        double mass_est = m_tas.get_state()[6] - norm_thrusts * m_talphas[i] / (m_isp * kep3::G0);
        double isp_est = norm_thrusts * m_talphas[i] / (-kep3::G0 * (m_tas.get_state()[6] - mass_est ));
        if (mass_est < mass_thresh) { // Set Isp to zero
            // fmt::print("Warning Mismatch: sf hf leg will run out of mass, setting Isp to inf. Mass estimate {} m0 {} T {} tof {}\n", mass_est, m_tas.get_state()[6], norm_thrusts, prop_seg_duration);
            *(m_tas.get_pars_data()+1l) = isp_est * kep3::G0;
        } else {
            auto [status, min_h, max_h, nsteps, _1, _2] = m_tas.propagate_until(prop_seg_temp);
            if (status != heyoka::taylor_outcome::time_limit) { // LCOV_EXCL_START
                fmt::print("mismatch fwd: {}\n", status);
                break;
            } // LCOV_EXCL_STOP
        }
    }

    // Reset ISP
    *(m_tas.get_pars_data()+1l) = get_isp() * kep3::G0;

    // Set fwd final state
    std::vector<double> rvm_fwd_final = m_tas.get_state();

    // Backward pass
    // Final state
    // Set the Taylor Integration final conditions
    prop_seg_temp = m_tof;

    m_tas.set_time(prop_seg_temp);
    std::copy(m_rvmf.begin(), m_rvmf.end(), m_tas.get_state_data());

    // Loop through segments in backward pass of Sims-Flanagan transcription
    for (auto i = 0u; i < m_nseg_bck; ++i) {
        // Assign current_thrusts to Taylor adaptive integrator
        std::copy(m_thrusts.begin() + (m_nseg - (i + 1)) * 3l, m_thrusts.begin() + 3l * (m_nseg - i),
                  m_tas.get_pars_data() + 2l);
        // ... and integrate
        prop_seg_temp -= m_talphas[m_talphas.size() - (i+1)];
        auto [status, min_h, max_h, nsteps, _1, _2] = m_tas.propagate_until(prop_seg_temp);
        if (status != heyoka::taylor_outcome::time_limit) { // LCOV_EXCL_START
            // fmt::print("mismatch bck: {}\n", status);
            break;
        } // LCOV_EXCL_STOP
    }

    return {rvm_fwd_final[0] - m_tas.get_state()[0], rvm_fwd_final[1] - m_tas.get_state()[1],
            rvm_fwd_final[2] - m_tas.get_state()[2], rvm_fwd_final[3] - m_tas.get_state()[3],
            rvm_fwd_final[4] - m_tas.get_state()[4], rvm_fwd_final[5] - m_tas.get_state()[5],
            rvm_fwd_final[6] - m_tas.get_state()[6]};
}

std::vector<double> sims_flanagan_hf_alpha::compute_throttle_constraints() const
{
    std::vector<double> retval(m_nseg);
    for (decltype(m_throttles.size()) i = 0u; i < m_nseg; ++i) {
        retval[i] = m_throttles[3 * i] * m_throttles[3 * i] + m_throttles[3 * i + 1] * m_throttles[3 * i + 1]
                    + m_throttles[3 * i + 2] * m_throttles[3 * i + 2] - 1.;
    }
    return retval;
}

std::vector<double> sims_flanagan_hf_alpha::compute_constraints() const
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
    // talphas
    return retval;
}

// std::vector<double> sims_flanagan_hf_alpha::set_and_compute_constraints(const std::vector<double> &chromosome)
// {
//     std::array<double, 7> rvms{};
//     std::copy(chromosome.begin(), chromosome.begin() + 7, rvms.begin());

//     std::vector<double> talphas(m_nseg);
//     std::copy(chromosome.begin() + 7, chromosome.begin() + 7 + m_nseg, talphas.begin());


//     std::vector<double> throttles(m_nseg * 3l);
//     std::copy(chromosome.begin() + 7 + m_nseg, chromosome.begin() + 7 + m_nseg + m_nseg * 3l, throttles.begin());
//     std::array<double, 7> rvmf{};
//     std::copy(chromosome.begin() + 7 + m_nseg + m_nseg * 3l, chromosome.begin() + 7 + m_nseg + m_nseg * 3l + 7, rvmf.begin());
//     double time_of_flight = chromosome[(7 + m_nseg + m_nseg * 3 + 7 + 1) - 1];
//     // Set relevant quantities before evaluating constraints
//     set(rvms, throttles, talphas, rvmf, time_of_flight);
//     // Evaluate and return constraints
//     return compute_constraints();
// }

std::vector<std::vector<double>> sims_flanagan_hf_alpha::get_state_history(unsigned grid_points_per_segment) const
{

    // Get time grid
    const double prop_seg_duration = (m_tof / m_nseg);
    std::vector<double> leg_time_grid;
    // Initial time
    double timestep = 0.0;
    leg_time_grid.push_back(timestep);
    for (decltype(m_nseg) j = 0; j < m_nseg; ++j) {
        for (decltype(m_nseg) i = 0; i < grid_points_per_segment-1; ++i) {
            timestep += m_talphas[j] / (grid_points_per_segment - 1);
            leg_time_grid.push_back(timestep);
        }
    }

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
        std::copy(m_thrusts.begin() + i * 3l, m_thrusts.begin() + 3 * (i + 1l), m_tas.get_pars_data() + 2l);

        // Current leg time grid
        std::copy(leg_time_grid.begin() + i * (grid_points_per_segment - 1l),
                  leg_time_grid.begin() + (i + 1l) * (grid_points_per_segment - 1l) + 1l,
                  current_leg_time_grid.begin());

        m_tas.set_time(current_leg_time_grid.at(0));
        // ... and integrate
        auto [status, min_h, max_h, nsteps, _1, output_states] = m_tas.propagate_grid(current_leg_time_grid);

        if (status != heyoka::taylor_outcome::time_limit) { // LCOV_EXCL_START
            fmt::print("thrust: [{}, {}, {}]\n", *(m_tas.get_pars_data()+2l), *(m_tas.get_pars_data() + 3l),
                       *(m_tas.get_pars_data() + 4l));
            fmt::print(": {}\n", status);
            fmt::print("state: {}\n", m_tas.get_state());
            fmt::print("reached time: {}\n", m_tas.get_time());
            fmt::print("requested time: {}\n", (i + 1) * prop_seg_duration);
            throw std::domain_error(
                "stark_problem: failure to reach the final time requested during a propagation."); 
        } // LCOV_EXCL_STOP
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
        std::copy(m_thrusts.begin() + (m_nseg - (i + 1)) * 3l, m_thrusts.begin() + 3l * (m_nseg - i),
                  m_tas.get_pars_data() + 2l);
        // Current leg time grid
        std::reverse_copy(leg_time_grid.begin() + (m_nseg - (i + 1l)) * (grid_points_per_segment - 1l),
                          leg_time_grid.begin() + (m_nseg - i) * (grid_points_per_segment - 1l) + 1l,
                          back_time_grid.begin());

        m_tas.set_time(back_time_grid.at(0));

        // ... and integrate
        auto [status, min_h, max_h, nsteps, _1, output_states] = m_tas.propagate_grid(back_time_grid);

        
        if (status != heyoka::taylor_outcome::time_limit) { // LCOV_EXCL_START
            fmt::print("thrust: [{}, {}, {}]\n", *(m_tas.get_pars_data()+2l), *(m_tas.get_pars_data() + 3l),
                       *(m_tas.get_pars_data() + 4l));
            fmt::print(": {}\n", status);
            fmt::print("state: {}\n", m_tas.get_state());
            fmt::print("reached time: {}\n", m_tas.get_time());
            fmt::print("requested time: {}\n", (i + 1) * prop_seg_duration);
            throw std::domain_error(
                "stark_problem: failure to reach the final time requested during a propagation."); 
        } // LCOV_EXCL_STOP
        output_per_seg[m_nseg - 1 - i] = output_states;
    }

    return output_per_seg;
}

std::ostream &operator<<(std::ostream &s, const sims_flanagan_hf_alpha &sf)
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
    s << fmt::format("Mismatch constraints: {}\n", sf.compute_mismatch_constraints());
    s << fmt::format("Throttle constraints: {}\n\n", sf.compute_throttle_constraints());
    return s;
}

} // namespace kep3::leg