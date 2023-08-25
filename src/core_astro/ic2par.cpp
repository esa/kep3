// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>

#include <boost/math/constants/constants.hpp>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xtensor.hpp>

#include <kep3/core_astro/ic2par.hpp>

using boost::math::constants::pi;
using xt::linalg::cross;
using xt::linalg::dot;

namespace kep3 {

// r,v,mu -> keplerian osculating elements [a,e,i,W,w,E]. The last
// is the eccentric anomaly or the Gudermannian according to e. The
// semi-major axis a is positive for ellipses, negative for hyperbolae.
// The anomalies W, w and E are in [0, 2pi]. Inclination is in [0, pi].

std::array<double, 6> ic2par(const std::array<double, 3> &rin,
                             const std::array<double, 3> &vin, double mu) {
  // Return value
  std::array<double, 6> retval{};
  // 0 - We prepare a few xtensor constants.
  xt::xtensor_fixed<double, xt::xshape<3>> k{0.0, 0.0, 1.0};
  xt::xtensor_fixed<double, xt::xshape<3>> r0 = xt::adapt(rin);
  xt::xtensor_fixed<double, xt::xshape<3>> v0 = xt::adapt(vin);

  // 1 - We compute the orbital angular momentum vector
  auto h = cross(r0, v0); // h = r0 x v0

  // 2 - We compute the orbital parameter
  auto p = dot(h, h) / mu; // p = h^2 / mu

  // 3 - We compute the vector of the node line
  // This operation is singular when inclination is zero, in which case the
  // Keplerian orbital parameters are not well defined
  auto n = cross(k, h);
  n = n / xt::linalg::norm(n); // n = (k x h) / |k x h|

  // 4 - We compute the eccentricity vector
  auto R0 = xt::linalg::norm(r0);
  auto evett = cross(v0, h) / mu - r0 / R0; // e = (v x h)/mu - r0/R0;

  // The eccentricity is calculated and stored as the second orbital element
  retval[1] = xt::linalg::norm(evett);

  // The semi-major axis (positive for ellipses, negative for hyperbolas) is
  // calculated and stored as the first orbital element a = p / (1 - e^2)
  retval[0] = p(0) / (1 - retval[1] * retval[1]);

  // Inclination is calculated and stored as the third orbital element
  // i = acos(hy/h)
  retval[2] = std::acos(h(2) / xt::linalg::norm(h));

  // Argument of pericentrum is calculated and stored as the fifth orbital
  // elemen.t w = acos(n.e)\|n||e|
  auto temp = dot(n, evett);
  retval[4] = std::acos(temp(0) / retval[1]);
  if (evett(2) < 0) {
    retval[4] = 2 * pi<double>() - retval[4];
  }
  // Argument of longitude is calculated and stored as the fourth orbital
  // element
  retval[3] = std::acos(n(0));
  if (n(1) < 0) {
    retval[3] = 2 * boost::math::constants::pi<double>() - retval[3];
  }

  // 4 - We compute ni: the true anomaly (in 0, 2*PI)
  temp = dot(evett, r0);
  auto ni = std::acos(temp(0) / retval[1] / R0);

  temp = dot(r0, v0);
  if (temp(0) < 0.0) {
    ni = 2 * boost::math::constants::pi<double>() - ni;
  }

  // Eccentric anomaly or the gudermannian is calculated and stored as the
  // sixth orbital element
  if (retval[1] < 1.0) {
    retval[5] = 2.0 * atan(sqrt((1 - retval[1]) / (1 + retval[1])) *
                           tan(ni / 2.0)); // algebraic Kepler's equation
  } else {
    retval[5] =
        2.0 * atan(sqrt((retval[1] - 1) / (retval[1] + 1)) *
                   tan(ni / 2.0)); // algebraic equivalent of Kepler's
                                   // equation in terms of the Gudermannian
  }
  return retval;
}
} // namespace kep3