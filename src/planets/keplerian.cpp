// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <string>

#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>

#include <kep3/core_astro/kepler_equations.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/epoch.hpp>
#include <kep3/planets/keplerian.hpp>
#include <utility>

namespace kep3::udpla {

constexpr double pi{boost::math::constants::pi<double>()};
const double DAY2SEC = 24.*60.*60.;

keplerian::keplerian(const epoch &ref_epoch,
                     const std::array<std::array<double, 3>, 2> &pos_vel,
                     double mu_central_body, std::string name,
                     std::array<double, 3> added_params)
    : m_ref_epoch(ref_epoch), m_name(std::move(name)),
      m_mu_central_body(mu_central_body), m_mu_self(added_params[0]),
      m_radius(added_params[1]), m_ellipse(), m_safe_radius(added_params[2]),
      m_pos_vel_0(pos_vel) {
  double R =
      std::sqrt(pos_vel[0][0] * pos_vel[0][0] + pos_vel[0][1] * pos_vel[0][1] +
                pos_vel[0][2] * pos_vel[0][2]);
  double en = (pos_vel[1][0] * pos_vel[1][0] + pos_vel[1][1] * pos_vel[1][1] +
              pos_vel[1][2] * pos_vel[1][2]) /
                 2. -
             mu_central_body / R;
  (en>0) ? m_ellipse=false : m_ellipse=true;
}

std::array<std::array<double, 3>, 2>
keplerian::eph(const kep3::epoch &ep) const {
  // 1 - We compute the dt
  double dt = (ep.mjd2000() - m_ref_epoch.mjd2000()) * DAY2SEC;
  // 2 - We propagate (make a copy as we do not want to change m_pos_vel_0)
  auto retval(m_pos_vel_0);
  kep3::propagate_lagrangian(retval, dt, m_mu_central_body);
  return retval;
}

const std::array<double, 6> keplerian::default_elements = {1.0, 0.0, 0.0,
                                                           0.0, 0.0, 0.0};

} // namespace kep3::udpla
