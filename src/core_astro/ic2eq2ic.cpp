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

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xadapt.hpp>

#include <kep3/core_astro/ic2eq2ic.hpp>

using xt::linalg::cross;
using xt::linalg::dot;

namespace kep3 {

constexpr double half_pi{boost::math::constants::half_pi<double>()};

std::array<double, 6> ic2eq(const std::array<std::array<double, 3>, 2> &pos_vel,
                            double mu, bool retrogade) {
  {
    // Switch between the element types.
    int I = 0;
    if (retrogade) {
      I = -1;
    } else {
      I = 1;
    }

    // 0 - We prepare a few xtensor constants.
    auto r0 = xt::adapt(pos_vel[0]);
    auto v0 = xt::adapt(pos_vel[1]);
    // The equinoctial reference frame
    xt::xtensor_fixed<double, xt::xshape<3>> fv = {0.0, 0.0, 0.0};
    xt::xtensor_fixed<double, xt::xshape<3>> gv = {0.0, 0.0, 0.0};

    // angular momentum
    auto ang = cross(r0, v0);

    // 0 - We compute the semi-major axis
    double R0 = norm(r0)(0);
    double V0 = norm(v0)(0);
    auto a = std::abs(1. / (2. / R0 - V0 * V0 / mu));

    // 1 - We compute the equinoctial frame
    auto w = cross(r0, v0);
    w = w / xt::linalg::norm(w);

    double k = w(0) / (1 + I * w(2));
    double h = -w(1) / (1 + I * w(2));
    double den = k * k + h * h + 1;
    fv(0) = (1. - k * k + h * h) / den;
    fv(1) = (2. * k * h) / den;
    fv(2) = (-2. * I * k) / den;

    gv(0) = (2. * I * k * h) / den;
    gv(1) = (1. + k * k - h * h) * I / den;
    gv(2) = (2. * h) / den;

    // 2 - We compute evett: the eccentricity vector
    auto evett = cross(v0, ang) / mu - r0 / R0; // e = (v x h)/mu - r0/R0;

    double g = dot(evett, gv)(0);
    double f = dot(evett, fv)(0);
    double ecc = norm(evett)(0);

    // 3 - We compute the true longitude L
    // This solution is certainly not the most elegant, but it works and will
    // never be singular.

    double det1 = (gv(1) * fv(0) - fv(1) * gv(0)); // xy
    double det2 = (gv(2) * fv(0) - fv(2) * gv(0)); // xz
    double det3 = (gv(2) * fv(1) - fv(2) * gv(1)); // yz
    double max = std::max({std::abs(det1), std::abs(det2), std::abs(det3)});

    double X = 0., Y = 0.;
    if (std::abs(det1) == max) {
      X = (gv(1) * r0(0) - gv(0) * r0(1)) / det1;
      Y = (-fv(1) * r0(0) + fv(0) * r0(1)) / det1;
    } else if (std::abs(det2) == max) {
      X = (gv(2) * r0(0) - gv(0) * r0(2)) / det2;
      Y = (-fv(2) * r0(0) + fv(0) * r0(2)) / det2;
    } else {
      X = (gv(2) * r0(1) - gv(1) * r0(2)) / det3;
      Y = (-fv(2) * r0(1) + fv(1) * r0(2)) / det3;
    }

    double L = std::atan2(Y / R0, X / R0);

    // 5 - We assign the results
    return {a * (1. - ecc * ecc),f, g, h, k, L};
  }
}

} // namespace kep3