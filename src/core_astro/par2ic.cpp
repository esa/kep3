// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include <exception>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>

#include <boost/math/constants/constants.hpp>

#include <kep3/core_astro/par2ic.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

namespace kep3 {

// keplerian osculating elements [a,e,i,W,w,f] -> r,v.
// The last osculating elements is the true anomaly.
// The semi-major axis a needs to be positive
// for ellipses, negative for hyperbolae.
// The anomalies W, w and E must be in [0, 2pi] and inclination in [0, pi].

std::array<std::array<double, 3>, 2> par2ic(const std::array<double, 6> &par,
                                            double mu) {
  // Return values
  std::array<double, 3> pos{};
  std::array<double, 3> vel{};
  auto pos_xt = xt::adapt(pos, {3u, 1u});
  auto vel_xt = xt::adapt(vel, {3u, 1u});

  // Rename some variables for readibility
  double sma = par[0];
  double ecc = par[1];
  double inc = par[2];
  double omg = par[3];
  double omp = par[4];
  double f = par[5];

  if (sma * (1 - ecc) < 0) {
    throw std::domain_error("par2ic was called with ecc and sma not compatible "
                            "with the convention a<0 -> e>1 [a>0 -> e<1].");
  }
  double cosf = std::cos(f);
  if (ecc > 1 && cosf < -1 / ecc) {
    throw std::domain_error("par2ic was called for an hyperbola but the true "
                            "anomaly is beyond asymptotes (cosf<-1/e)");
  }

  // 1 - We start by evaluating position and velocity in the perifocal reference
  // system
  double p = sma * (1.0 - ecc * ecc);
  double r = p / (1.0 + ecc * std::cos(f));
  double h = std::sqrt(p * mu);
  double sinf = std::sin(f);
  double x_per = r * cosf;
  double y_per = r * sinf;
  double xdot_per = -mu / h * sinf;
  double ydot_per = mu / h * (ecc + cosf);

  // 2 - We then built the rotation matrix from perifocal reference frame to
  // inertial
  double cosomg = std::cos(omg);
  double cosomp = std::cos(omp);
  double sinomg = std::sin(omg);
  double sinomp = std::sin(omp);
  double cosi = std::cos(inc);
  double sini = std::sin(inc);

  xt::xtensor_fixed<double, xt::xshape<3, 3>> R = {
      {cosomg * cosomp - sinomg * sinomp * cosi,
       -cosomg * sinomp - sinomg * cosomp * cosi, sinomg * sini},
      {sinomg * cosomp + cosomg * sinomp * cosi,
       -sinomg * sinomp + cosomg * cosomp * cosi, -cosomg * sini},
      {sinomp * sini, cosomp * sini, cosi}};

  // 3 - We end by transforming according to this rotation matrix
  xt::xtensor_fixed<double, xt::xshape<3, 1>> pos_per{{x_per}, {y_per}, {0.0}};
  xt::xtensor_fixed<double, xt::xshape<3, 1>> vel_per{
      {xdot_per}, {ydot_per}, {0.0}};

  // The following lines, since use xtensors adapted to pos and vel, will change
  // pos and vel.
  pos_xt = xt::linalg::dot(R, pos_per);
  vel_xt = xt::linalg::dot(R, vel_per);

  return {pos, vel};
}

} // namespace kep3