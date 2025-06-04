// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <array>
#include <cmath>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xadapt.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/linalg.hpp>

namespace kep3
{

using xt::linalg::dot;

inline double safe_acos(double x)
{
    return std::acos(std::clamp(x, -1.0, 1.0));
}

// Type alias for clarity
using vec3 = std::array<double, 3>;

inline vec3 operator+(const vec3 &a, const vec3 &b)
{
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

inline vec3 operator-(const vec3 &a, const vec3 &b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

inline vec3 operator*(const vec3 &a, double s)
{
    return {a[0] * s, a[1] * s, a[2] * s};
}

inline vec3 operator*(double s, const vec3 &a)
{
    return a * s;
}

inline vec3 operator/(const vec3 &a, double s)
{
    return {a[0] / s, a[1] / s, a[2] / s};
}

inline double dot(const vec3 &a, const vec3 &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline double norm(const vec3 &a)
{
    return std::sqrt(dot(a, a));
}

inline vec3 cross(const vec3 &a, const vec3 &b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// r,v,mu -> keplerian osculating elements [a,e,i,W,w,f]. The last
// is the true anomaly. The semi-major axis a is positive for ellipses, negative
// for hyperbolae. The anomalies W, w, f are in [0, 2pi]. Inclination is in [0,
// pi].

std::array<double, 6> ic2par(const std::array<vec3, 2> &pos_vel, double mu)
{
    std::array<double, 6> retval{};
    const vec3 &r = pos_vel[0];
    const vec3 &v = pos_vel[1];

    vec3 h = cross(r, v);
    double h_norm = norm(h);
    double p = dot(h, h) / mu;

    // Node vector n = k × h
    vec3 k = {0.0, 0.0, 1.0};
    vec3 n = cross(k, h);
    double n_norm = norm(n);
    n = n / n_norm;

    // Eccentricity vector
    vec3 e_vec = cross(v, h) / mu - r / norm(r);
    double e = norm(e_vec);

    retval[0] = p / (1.0 - e * e);        // semi-major axis
    retval[1] = e;                        // eccentricity
    retval[2] = safe_acos(h[2] / h_norm); // inclination

    // Argument of pericenter ω
    double cos_w = dot(n, e_vec) / e;
    retval[4] = safe_acos(cos_w);
    if (e_vec[2] < 0.0)
        retval[4] = 2.0 * pi - retval[4];

    // Longitude of ascending node Ω
    retval[3] = safe_acos(n[0]);
    if (n[1] < 0.0)
        retval[3] = 2.0 * pi - retval[3];

    // True anomaly f
    double cosf = dot(e_vec, r) / (e * norm(r));
    double sinf = dot(r, v) * h_norm / (e * norm(r) * mu);
    double f = std::atan2(sinf, cosf);
    if (f < 0.0)
        f += 2.0 * pi;
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