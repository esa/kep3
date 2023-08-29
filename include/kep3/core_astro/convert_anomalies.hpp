// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_CONVERT_ANOMALIES_H
#define kep3_CONVERT_ANOMALIES_H

#include <cmath>

#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>

#include <kep3/core_astro/kepler_equations.hpp>

namespace kep3 {

// mean to eccentric (only ellipses) e<1. Preserves the sign and integer number
// of revolutions.
inline double m2e(double M, double ecc) {
  // We use as intial guess the Mean Anomaly
  // (tests indicated that any higher order expansion does not really improve)
  double IG = M;
  const int digits = std::numeric_limits<double>::digits;
  std::uintmax_t max_iter = 100u;
  double sol = boost::math::tools::halley_iterate(
      [M, ecc](double E) {
        return std::make_tuple(kepE(E, M, ecc), d_kepE(E, ecc),
                               dd_kepE(E, ecc));
      },
      IG, IG - boost::math::constants::pi<double>(),
      IG + boost::math::constants::pi<double>(), digits, max_iter);
  if (max_iter == 100u) {
    throw std::domain_error(
        "Maximum number of iterations exceeded when solving Kepler's "
        "equation for the eccentric anomaly in m2e.");
  }
  return sol;
}
// eccentric to mean (only ellipses) e<1
inline double e2m(double E, double ecc) { return (E - ecc * std::sin(E)); }

// eccentric to true (only ellipses) e<1 (returns in range [-pi,pi])
inline double e2f(double E, double ecc) {
  return 2 * std::atan(std::sqrt((1 + ecc) / (1 - ecc)) * std::tan(E / 2));
}
// true to eccentric (only ellipses) e<1 (returns in range [-pi,pi])
inline double f2e(double f, double ecc) {
  return 2 * std::atan(std::sqrt((1 - ecc) / (1 + ecc)) * std::tan(f / 2));
}

// mean to true (only ellipses) e<1 (returns in range [-pi,pi])
inline double m2f(double M, double ecc) {
  return e2f(m2e(M, ecc), ecc);
}
// true to mean (only ellipses) e<1 (returns in range [-pi,pi])
inline double f2m(double f, double ecc) {
  return e2m(f2e(f, ecc), ecc);
}

// gudermannian to true (only hyperbolas) e>1 (returns in range [-pi,pi])
inline double zeta2f(double f, double ecc) {
  return 2 * std::atan(std::sqrt((1 + ecc) / (ecc - 1)) * std::tan(f / 2));
}
// true to gudermannian (only hyperbolas) e>1 (returns in range [-pi,pi])
inline double f2zeta(double zeta, double ecc) {
  return 2 * std::atan(std::sqrt((ecc - 1) / (1 + ecc)) * std::tan(zeta / 2));
}
} // namespace kep3
#endif // kep3_TOOLBOX_M2E_H
