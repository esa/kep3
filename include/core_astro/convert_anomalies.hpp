/*****************************************************************************
 *   Copyright (C) 2023 The pykep development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://gitter.im/esa/pykep                                             *
 *   https://github.com/esa/pykep                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#ifndef kep3_CONVERT_ANOMALIES_H
#define kep3_CONVERT_ANOMALIES_H

#include <cmath>

#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>

#include <core_astro/kepler_equations.hpp>

namespace kep3 {

// mean to eccentric (will return a result in [0, 2pi]) (only ellipses) e<1
inline double m2e(double M, double ecc) {
  // We sent M in the [0, 2pi] range
  M = std::fmod(M, 2. * boost::math::constants::pi<double>());
  if (M < 0) {
    M += 2. * boost::math::constants::pi<double>();
  }
  // We use as intial guess the expansion on e of Kepler's equation
  double IG = M; //+ ecc * std::sin(M);
  const int digits = std::numeric_limits<double>::digits;
  double sol = boost::math::tools::halley_iterate(
      [M, ecc](double E) {
        return std::make_tuple(kepE(E, M, ecc), d_kepE(E, ecc), dd_kepE(E, ecc));
      },
      IG, 0., 2. * boost::math::constants::pi<double>(), digits);
  return sol;
}
// eccentric to mean (only ellipses) e<1
inline double e2m(double E, double e) { return (E - e * std::sin(E)); }
// eccentric to true (only ellipses) e<1
inline double e2f(double E, double e) {
  return 2 * std::atan(std::sqrt((1 + e) / (1 - e)) * std::tan(E / 2));
}
// true to eccentric (only ellipses) e<1
inline double f2e(double f, double e) {
  return 2 * std::atan(std::sqrt((1 - e) / (1 + e)) * std::tan(f / 2));
}
// gudermannian to true (only hyperbolas) e>1
inline double zeta2f(double E, double e) {
  return 2 * std::atan(std::sqrt((1 + e) / (e - 1)) * std::tan(E / 2));
}
// true to gudermannian (only hyperbolas) e>1
inline double f2zeta(double zeta, double e) {
  return 2 * std::atan(std::sqrt((e - 1) / (1 + e)) * std::tan(zeta / 2));
}
} // namespace kep3
#endif // kep3_TOOLBOX_M2E_H
