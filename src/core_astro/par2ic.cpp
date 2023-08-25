// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>

#include <boost/math/constants/constants.hpp>

#include <kep3/core_astro/par2ic.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

namespace kep3 {

constexpr double pi4{boost::math::constants::quarter_pi<double>()};

// keplerian osculating elements [a,e,i,W,w,E] -> r,v.
// The last osculating elements needs to be the eccentric anomaly or
// the Gudermannian according to e. The semi-major axis a needs to be positive
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
  double acc = par[0];
  double ecc = par[1];
  double inc = par[2];
  double omg = par[3];
  double omp = par[4];
  double EA = par[5];

  // TODO(darioizzo): Check a<0 if e>1

  // 1 - We start by evaluating position and velocity in the perifocal reference
  // system
  double xper = 0., yper = 0., xdotper = 0., ydotper = 0.;
  if (ecc < 1.0) // EA is the eccentric anomaly
  {
    double b = acc * std::sqrt(1 - ecc * ecc);
    double n = std::sqrt(mu / (acc * acc * acc));
    xper = acc * (std::cos(EA) - ecc);
    yper = b * std::sin(EA);
    xdotper = -(acc * n * std::sin(EA)) / (1 - ecc * std::cos(EA));
    ydotper = (b * n * std::cos(EA)) / (1 - ecc * std::cos(EA));
  } else // EA is the Gudermannian
  {
    double b = -acc * std::sqrt(ecc * ecc - 1);
    double n = std::sqrt(-mu / (acc * acc * acc));

    double dNdZeta = ecc * (1 + std::tan(EA) * std::tan(EA)) -
                     (0.5 + 0.5 * std::pow(std::tan(0.5 * EA + pi4), 2)) /
                         std::tan(0.5 * EA + pi4);

    xper = acc / std::cos(EA) - acc * ecc;
    yper = b * std::tan(EA);

    xdotper = acc * tan(EA) / cos(EA) * n / dNdZeta;
    ydotper = b / pow(cos(EA), 2) * n / dNdZeta;
  }

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
  xt::xtensor_fixed<double, xt::xshape<3, 1>> temp1{{xper}, {yper}, {0.0}};
  xt::xtensor_fixed<double, xt::xshape<3, 1>> temp2{
      {xdotper}, {ydotper}, {0.0}};

  pos_xt = R * temp1;
  vel_xt = R * temp2;
  fmt::print("\npar2ic: {}, {}", pos_xt, vel_xt);
  fmt::print("\npar2ic: {}, {}", pos, vel);

  return {pos, vel};
}

} // namespace kep3