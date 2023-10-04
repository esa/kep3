// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include <limits>

#include <chrono>
#include <stdexcept>
#include <string>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/core_astro/eq2par2eq.hpp>
#include <kep3/core_astro/ic2eq2ic.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/epoch.hpp>
#include <kep3/planet.hpp>
#include <kep3/planets/keplerian.hpp>

namespace kep3::udpla
{

keplerian::keplerian(const epoch &ref_epoch, const std::array<std::array<double, 3>, 2> &pos_vel,
                     double mu_central_body, std::string name, std::array<double, 3> added_params)
    : m_ref_epoch(ref_epoch), m_name(std::move(name)), m_mu_central_body(mu_central_body), m_mu_self(added_params[0]),
      m_radius(added_params[1]), m_safe_radius(added_params[2]), m_period(), m_ellipse(), m_pos_vel_0(pos_vel)
{
    double R = std::sqrt(pos_vel[0][0] * pos_vel[0][0] + pos_vel[0][1] * pos_vel[0][1] + pos_vel[0][2] * pos_vel[0][2]);
    double en = (pos_vel[1][0] * pos_vel[1][0] + pos_vel[1][1] * pos_vel[1][1] + pos_vel[1][2] * pos_vel[1][2]) / 2.
                - mu_central_body / R;
    (en > 0) ? m_ellipse = false : m_ellipse = true;
    double a = -m_mu_central_body / 2. / en;
    if (m_ellipse) {
        m_period = kep3::pi * 2. * std::sqrt(a * a * a / m_mu_central_body);
    } else {
        m_period = std::numeric_limits<double>::quiet_NaN();
    }
}

keplerian::keplerian(const epoch &ref_epoch, const std::array<double, 6> &par_in, double mu_central_body,
                     std::string name, std::array<double, 3> added_params, kep3::elements_type el_type)
    : m_ref_epoch(ref_epoch), m_name(std::move(name)), m_mu_central_body(mu_central_body), m_mu_self(added_params[0]),
      m_radius(added_params[1]), m_safe_radius(added_params[2]), m_period(), m_ellipse(), m_pos_vel_0()
{
    // orbital parameters a,e,i,W,w,f will be stored here
    std::array<double, 6> par(par_in);
    // we convert according to the chosen input
    switch (el_type) {
        case kep3::elements_type::KEP_F:
            break;
        case kep3::elements_type::KEP_M:
            if (par[0] < 0) {
                throw std::logic_error("Mean anomaly is only available for ellipses.");
            }
            par[5] = kep3::m2f(par[5], par[1]);
            break;
        case kep3::elements_type::MEQ:
            par = kep3::eq2par(par, false);
            break;
        case kep3::elements_type::MEQ_R:
            par = kep3::eq2par(par, true);
            break;
        default:
            throw std::logic_error("You should not go here!");
    }

    if (par[0] * (1 - par[1]) <= 0) {
        throw std::domain_error("A Keplerian planet constructor was called with with non compatible "
                                "a,e:"
                                "The following must hold: a<0 -> e>1 [a>0 -> e<1].");
    }

    m_pos_vel_0 = kep3::par2ic(par, mu_central_body);

    (par[0] < 0) ? m_ellipse = false : m_ellipse = true;
    if (m_ellipse) {
        m_period = kep3::pi * 2. * std::sqrt(par[0] * par[0] * par[0] / m_mu_central_body);
    } else {
        m_period = std::numeric_limits<double>::quiet_NaN();
    }
}

std::array<std::array<double, 3>, 2> keplerian::eph(const kep3::epoch &ep) const
{
    // 1 - We compute the dt
    auto dt = epoch::as_sec(ep - m_ref_epoch);
    // 2 - We propagate (make a copy as we do not want to change m_pos_vel_0)
    auto retval(m_pos_vel_0);
    kep3::propagate_lagrangian(retval, dt, m_mu_central_body);
    return retval;
}

std::string keplerian::get_name() const
{
    return m_name;
}

double keplerian::get_mu_central_body() const
{
    return m_mu_central_body;
}

double keplerian::get_mu_self() const
{
    return m_mu_self;
}

double keplerian::get_radius() const
{
    return m_radius;
}

double keplerian::get_safe_radius() const
{
    return m_safe_radius;
}

double keplerian::period(const kep3::epoch &) const
{
    return m_period;
}

kep3::epoch keplerian::get_ref_epoch() const
{
    return m_ref_epoch;
}

std::array<double, 6> keplerian::elements(kep3::elements_type el_type) const
{
    std::array<double, 6> retval{};
    switch (el_type) {
        case kep3::elements_type::KEP_F:
            retval = kep3::ic2par(m_pos_vel_0, m_mu_central_body);
            break;
        case kep3::elements_type::KEP_M:
            if (!m_ellipse) {
                throw std::logic_error("Mean anomaly is only available for ellipses.");
            }
            retval = kep3::ic2par(m_pos_vel_0, m_mu_central_body);
            retval[5] = kep3::f2m(retval[5], retval[1]);
            break;
        case kep3::elements_type::MEQ:
            retval = kep3::ic2eq(m_pos_vel_0, m_mu_central_body);
            break;
        case kep3::elements_type::MEQ_R:
            retval = kep3::ic2eq(m_pos_vel_0, m_mu_central_body, true);
            break;
        default:
            throw std::logic_error("You should not go here!");
    }
    return retval;
}

std::string keplerian::get_extra_info() const
{
    auto par = elements();
    std::string retval
        = fmt::format("Keplerian planet elements: \n") + fmt::format("Semi major axis (AU): {}\n", par[0] / kep3::AU)
          + fmt::format("Eccentricity: {}\n", par[1]) + fmt::format("Inclination (deg.): {}\n", par[2] * kep3::RAD2DEG)
          + fmt::format("Big Omega (deg.): {}\n", par[3] * kep3::RAD2DEG)
          + fmt::format("Small omega (deg.): {}\n", par[4] * kep3::RAD2DEG)
          + fmt::format("True anomly (deg.): {}\n", par[5] * kep3::RAD2DEG);
    if (m_ellipse) {
        retval += fmt::format("Mean anomly (deg.): {}\n", kep3::f2m(par[5], par[1]) * kep3::RAD2DEG);
    }
    retval += fmt::format("Elements reference epoch (MJD2000): {}\n", m_ref_epoch)
              + fmt::format("Elements reference epoch (date): {}\n", m_ref_epoch)
              + fmt::format("r at ref. = {}\n", m_pos_vel_0[0]) + fmt::format("v at ref. = {}\n", m_pos_vel_0[1]);
    return retval;
}

std::ostream &operator<<(std::ostream &os, const kep3::udpla::keplerian &udpla)
{
    os << udpla.get_extra_info() << std::endl;
    return os;
}

} // namespace kep3::udpla

// NOLINTNEXTLINE
TANUKI_S11N_WRAP_EXPORT_IMPLEMENT(kep3::udpla::keplerian, kep3::detail::planet_iface)
