// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <cmath>

#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>

#include <kep3/core_astro/kepler_equations.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/core_astro/special_functions.hpp>

namespace kep3 {

constexpr double pi{boost::math::constants::pi<double>()};

/// Lagrangian propagation
/**
 * This function propagates an initial Cartesian state for a time t assuming a
 * central body and a keplerian motion. Lagrange coefficients are used as basic
 * numerical technique. All units systems can be used, as long
 * as the input parameters are all expressed in the same system.
 */

void propagate_lagrangian(std::array<std::array<double, 3>, 2> &pos_vel_0,
                          const double dt, const double mu) {
  auto &[r0, v0] = pos_vel_0;
  double R = std::sqrt(r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2]);
  double V = std::sqrt(v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]);
  double energy = (V * V / 2 - mu / R);
  double a = -mu / 2.0 / energy; // will be negative for hyperbolae
  double sqrta = 0.;
  double F = 0., G = 0., Ft = 0., Gt = 0.;

  double sigma0 =
      (r0[0] * v0[0] + r0[1] * v0[1] + r0[2] * v0[2]) / std::sqrt(mu);

  if (a > 0) { // Solve Kepler's equation, elliptical case
    sqrta = std::sqrt(a);
    double DM = std::sqrt(mu / pow(a, 3)) * dt;

    // Solve Kepler Equation for ellipses in DE (eccentric anomaly difference)
    const int digits = std::numeric_limits<double>::digits;
    double DE = boost::math::tools::halley_iterate(
        [DM, sigma0, sqrta, a, R](double DE) {
          return std::make_tuple(kepDE(DE, DM, sigma0, sqrta, a, R),
                                 d_kepDE(DE, sigma0, sqrta, a, R),
                                 dd_kepDE(DE, sigma0, sqrta, a, R));
        },
        DM, DM - pi, DM + pi, digits);

    double r = a + (R - a) * cos(DE) + sigma0 * sqrta * sin(DE);

    // Lagrange coefficients
    F = 1 - a / R * (1 - cos(DE));
    G = a * sigma0 / sqrt(mu) * (1 - cos(DE)) + R * sqrt(a / mu) * sin(DE);
    Ft = -sqrt(mu * a) / (r * R) * sin(DE);
    Gt = 1 - a / r * (1 - cos(DE));
  } else { // Solve Kepler's equation, hyperbolic case
    sqrta = sqrt(-a);
    double DN = sqrt(-mu / pow(a, 3)) * dt;
    double IG = 0.;
    dt > 0 ? IG = 1 : IG = -1; // TODO(darioizzo): find a better initial guess.
                               // I tried with 0 and D (both have numercial
                               // problems and result in exceptions)

    // Solve Kepler Equation for ellipses in DH (hyperbolic anomaly difference)
    const int digits = std::numeric_limits<double>::digits;
    double DH = boost::math::tools::halley_iterate(
        [DN, sigma0, sqrta, a, R](double DH) {
          return std::make_tuple(kepDH(DH, DN, sigma0, sqrta, a, R),
                                 d_kepDH(DH, sigma0, sqrta, a, R),
                                 dd_kepDH(DH, sigma0, sqrta, a, R));
        },
        IG, IG - pi, IG + pi, digits);

    double r = a + (R - a) * cosh(DH) + sigma0 * sqrta * sinh(DH);

    // Lagrange coefficients
    F = 1 - a / R * (1 - cosh(DH));
    G = a * sigma0 / sqrt(mu) * (1 - cosh(DH)) + R * sqrt(-a / mu) * sinh(DH);
    Ft = -sqrt(-mu * a) / (r * R) * sinh(DH);
    Gt = 1 - a / r * (1 - cosh(DH));
  }

  double temp[3] = {r0[0], r0[1], r0[2]};
  for (int i = 0; i < 3; i++) {
    r0[i] = F * r0[i] + G * v0[i];
    v0[i] = Ft * temp[i] + Gt * v0[i];
  }
}

/// Universial Variales version
/**
 * This function has the same prototype as kep3::propagate_lgrangian, but
 * internally makes use of universal variables formulation for the Lagrange
 * Coefficients.
 */
void propagate_lagrangian_u(std::array<std::array<double, 3>, 2> &pos_vel_0,
                            const double dt, const double mu) {
  // If time is negative we need to invert time and velocities. Unlike the other
  // formulation of the propagate lagrangian we cannot rely on negative times to
  // automatically mean back-propagation
  double dt_copy = dt;
  auto &[r0, v0] = pos_vel_0;

  if (dt < 0) {
    dt_copy = -dt;
    v0[0] = -v0[0];
    v0[1] = -v0[1];
    v0[2] = -v0[2];
  }

  double F = 0., G = 0., Ft = 0., Gt = 0.;
  double R0 = sqrt(r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2]);
  double V0 = sqrt(v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]);
  // the reciprocal of the semi-major axis
  double alpha = 2 / R0 - V0 * V0 / mu;
  // initial radial velocity
  double VR0 = (r0[0] * v0[0] + r0[1] * v0[1] + r0[2] * v0[2]) / R0;

  // solve kepler's equation in universal variables
  double IG = 1;
  alpha > 0 ? IG = sqrt(mu) * dt_copy * std::abs(alpha)
            : IG = 1; // initial guess for the universal anomaly. For
                      // hyperbolas it is 1.... can be better?

  // Solve Kepler Equation in DS (univrsal anomaly difference)
  const int digits = std::numeric_limits<double>::digits;
  double DS = boost::math::tools::halley_iterate(
      [dt_copy, R0, VR0, alpha, mu](double DS) {
        return std::make_tuple(kepDS(DS, dt_copy, R0, VR0, alpha, mu),
                               d_kepDS(DS, R0, VR0, alpha, mu),
                               dd_kepDS(DS, R0, VR0, alpha, mu));
      },
      IG, IG - pi, IG + pi, digits);

  // evaluate the lagrangian coefficients F and G
  double S = stumpff_s(alpha * DS * DS);
  double C = stumpff_c(alpha * DS * DS);
  //
  double z = alpha * DS * DS;
  F = 1 - DS * DS / R0 * C;
  G = dt_copy - 1 / sqrt(mu) * DS * DS * DS * S;

  double r0_copy[3] = {r0[0], r0[1], r0[2]};
  // compute the final position
  r0[0] = F * r0[0] + G * v0[0];
  r0[1] = F * r0[1] + G * v0[1];
  r0[2] = F * r0[2] + G * v0[2];
  double RF = sqrt(r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2]);

  // compute the lagrangian coefficients Ft, Gt
  Ft = sqrt(mu) / RF / R0 * (z * S - 1) * DS;
  Gt = 1 - DS * DS / RF * C;

  // compute the final velocity
  v0[0] = Ft * r0_copy[0] + Gt * v0[0];
  v0[1] = Ft * r0_copy[1] + Gt * v0[1];
  v0[2] = Ft * r0_copy[2] + Gt * v0[2];

  if (dt < 0) {
    v0[0] = -v0[0];
    v0[1] = -v0[1];
    v0[2] = -v0[2];
  }
}

} // namespace kep3