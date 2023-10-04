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

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xadapt.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>

namespace kep3
{

using xt::linalg::cross;
using xt::linalg::dot;

// r,v,mu -> keplerian osculating elements [a,e,i,W,w,f]. The last
// is the true anomaly. The semi-major axis a is positive for ellipses, negative
// for hyperbolae. The anomalies W, w, f are in [0, 2pi]. Inclination is in [0,
// pi].

std::array<double, 6> ic2par(const std::array<std::array<double, 3>, 2> &pos_vel, double mu)
{
    // Return value
    std::array<double, 6> retval{};
    // 0 - We prepare a few xtensor constants.
    xt::xtensor_fixed<double, xt::xshape<3>> k{0.0, 0.0, 1.0};
    auto r0 = xt::adapt(pos_vel[0]);
    auto v0 = xt::adapt(pos_vel[1]);

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
    // i = acos(hy/h) in [0, pi]
    retval[2] = std::acos(h(2) / xt::linalg::norm(h));

    // Argument of pericentrum is calculated and stored as the fifth orbital
    // element. w = acos(n.e)\|n||e| in [0, 2pi]
    auto temp = dot(n, evett);
    retval[4] = std::acos(temp(0) / retval[1]);
    if (evett(2) < 0) {
        retval[4] = 2 * pi - retval[4];
    }
    // Argument of longitude is calculated and stored as the fourth orbital
    // element in [0, 2pi]
    retval[3] = std::acos(n(0));
    if (n(1) < 0) {
        retval[3] = 2 * pi - retval[3];
    }

    // 4 - We compute ni: the true anomaly in [0, 2pi]
    temp = dot(evett, r0);
    auto f = std::acos(temp(0) / retval[1] / R0);

    temp = dot(r0, v0);
    if (temp(0) < 0.0) {
        f = 2 * pi - f;
    }
    retval[5] = f;
    return retval;
}

// keplerian osculating elements [a,e,i,W,w,f] -> r,v.
// The last osculating elements is the true anomaly.
// The semi-major axis a needs to be positive
// for ellipses, negative for hyperbolae.
// The anomalies W, w and E must be in [0, 2pi] and inclination in [0, pi].

std::array<std::array<double, 3>, 2> par2ic(const std::array<double, 6> &par, double mu)
{
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

    xt::xtensor_fixed<double, xt::xshape<3, 3>> R
        = {{cosomg * cosomp - sinomg * sinomp * cosi, -cosomg * sinomp - sinomg * cosomp * cosi, sinomg * sini},
           {sinomg * cosomp + cosomg * sinomp * cosi, -sinomg * sinomp + cosomg * cosomp * cosi, -cosomg * sini},
           {sinomp * sini, cosomp * sini, cosi}};

    // 3 - We end by transforming according to this rotation matrix
    xt::xtensor_fixed<double, xt::xshape<3, 1>> pos_per{{x_per}, {y_per}, {0.0}};
    xt::xtensor_fixed<double, xt::xshape<3, 1>> vel_per{{xdot_per}, {ydot_per}, {0.0}};

    // The following lines, since use xtensors adapted to pos and vel, will change
    // pos and vel.
    pos_xt = xt::linalg::dot(R, pos_per);
    vel_xt = xt::linalg::dot(R, vel_per);

    return {pos, vel};
}
} // namespace kep3