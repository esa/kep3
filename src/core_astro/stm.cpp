// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "kep3/core_astro/kepler_equations.hpp"
#include <array>
#include <cmath>

#include <boost/math/tools/roots.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xadapt.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/core_astro/stm.hpp>

using xt::linalg::dot;
using xt::linalg::inv;
using mat31 = xt::xtensor_fixed<double, xt::xshape<3, 1>>;
using mat13 = xt::xtensor_fixed<double, xt::xshape<1, 3>>;
using mat33 = xt::xtensor_fixed<double, xt::xshape<3, 3>>;
using mat36 = xt::xtensor_fixed<double, xt::xshape<3, 6>>;
using mat66 = xt::xtensor_fixed<double, xt::xshape<6, 6>>;
using mat32 = xt::xtensor_fixed<double, xt::xshape<3, 2>>;
using mat62 = xt::xtensor_fixed<double, xt::xshape<6, 2>>;
using mat61 = xt::xtensor_fixed<double, xt::xshape<6, 1>>;
using mat16 = xt::xtensor_fixed<double, xt::xshape<1, 6>>;
using mat63 = xt::xtensor_fixed<double, xt::xshape<6, 3>>;

namespace kep3
{

template <std::size_t m, std::size_t k, std::size_t n>
xt::xtensor_fixed<double, xt::xshape<m, n>> _dot(const xt::xtensor_fixed<double, xt::xshape<m, k>> &A,
                                                 const xt::xtensor_fixed<double, xt::xshape<k, n>> &B)
{
    xt::xtensor_fixed<double, xt::xshape<m, n>> C{};
    for (decltype(m) i = 0u; i < m; ++i) {
        for (decltype(n) j = 0u; j < n; ++j) {
            C(i, j) = 0;
            for (decltype(k) l = 0u; l < k; ++l) {
                {
                    C(i, j) += A(i, l) * B(l, j);
                }
            }
        }
    }
    return C;
}

mat33 _skew(const mat31 &v)
{
    return {{0., -v(2, 0), v(1, 0)}, {v(2, 0), 0., -v(0, 0)}, {-v(1, 0), v(0, 0), 0.}};
}

mat31 _cross(const mat31 &v1, const mat31 &v2)
{
    return {{v1(1, 0) * v2(2, 0) - v2(1, 0) * v1(2, 0), v2(0, 0) * v1(2, 0) - v1(0, 0) * v2(2, 0),
             v1(0, 0) * v2(1, 0) - v2(0, 0) * v1(1, 0)}};
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
std::array<double, 36> stm_l(const std::array<std::array<double, 3>, 2> &pos_vel0,                        // NOLINT
                             const std::array<std::array<double, 3>, 2> &pos_velf, double tof, double mu, // NOLINT
                             double R0,                                                                   // NOLINT
                             double V02, double energy, double sigma0, double a, double s0, double c0,    // NOLINT
                             double DE,                                                                   // NOLINT
                             double F, double G,                                                          // NOLINT
                             double Ft, double Gt)                                                        // NOLINT
{
    // Create xtensor fixed arrays from input (we avoid adapt as its slower in this case since all is fixed size)
    // We use row vectors (not column) as its then more conventional for gradients as the differential goes to the
    // end
    mat13 r0 = {{pos_vel0[0][0]}, {pos_vel0[0][1]}, {pos_vel0[0][2]}};
    mat13 v0 = {{pos_vel0[1][0]}, {pos_vel0[1][1]}, {pos_vel0[1][2]}};
    mat13 rf = {{pos_velf[0][0]}, {pos_velf[0][1]}, {pos_velf[0][2]}};
    mat13 vf = {{pos_velf[1][0]}, {pos_velf[1][1]}, {pos_velf[1][2]}};

    // We seed the gradients with the initial dr0/dx0 and dv0/dx0
    mat36 dr0 = xt::zeros<double>({3, 6});
    mat36 dv0 = xt::zeros<double>({3, 6});
    dr0(0, 0) = 1;
    dr0(1, 1) = 1;
    dr0(2, 2) = 1;
    dv0(0, 3) = 1;
    dv0(1, 4) = 1;
    dv0(2, 5) = 1;

    // 1 - We start computing the differentials of basic quantities common at all eccentricities
    double Rf = std::sqrt(rf(0, 0) * rf(0, 0) + rf(0, 1) * rf(0, 1) + rf(0, 2) * rf(0, 2));
    double sqrtmu = std::sqrt(mu);
    mat16 dV02 = 2. * _dot(v0, dv0);
    mat16 dR0 = 1. / R0 * _dot(r0, dr0);
    mat16 denergy = 0.5 * dV02 + mu / R0 / R0 * dR0;
    mat16 dsigma0 = ((_dot(r0, dv0) + _dot(v0, dr0))) / sqrtmu;
    mat16 da = mu / 2. / energy / energy * denergy; // a = -mu / 2 / energy
    mat16 dF, dFt, dG, dGt;

    if (a > 0) { // ellipses
        double sqrta = std::sqrt(a);
        double sinDE = std::sin(DE);
        double cosDE = std::cos(DE);

        mat16 ds0 = dsigma0 / sqrta - 0.5 * sigma0 / sqrta / sqrta / sqrta * da; // s0 = sigma0 / sqrta
        mat16 dc0 = -1. / a * dR0 + R0 / a / a * da;                             // c0 = (1- R/a)
        mat16 dDM = -1.5 * sqrtmu * tof / std::pow(sqrta, 5) * da;               // M = sqrt(mu/a**3) tof
        mat16 dDE = (dDM - (1 - cosDE) * ds0 + sinDE * dc0) / (1 + s0 * sinDE - c0 * cosDE);
        mat16 dRf = (1 - cosDE + 0.5 / sqrta * sigma0 * sinDE) * da + cosDE * dR0
                    + (sigma0 * sqrta * cosDE - (R0 - a) * sinDE) * dDE
                    + sqrta * sinDE * dsigma0; // r = a + (r0 - a) * cosDE + sigma0 * sqrta * sinDE

        // 2 - We may now compute the differentials of the Lagrange coefficients
        dF = -(1 - cosDE) / R0 * da + a / R0 / R0 * (1 - cosDE) * dR0 - a / R0 * sinDE * dDE;
        dG = (1 - F) * (R0 * dsigma0 + sigma0 * dR0) - (sigma0 * R0) * dF + (sqrta * R0 * cosDE) * dDE
             + (sqrta * sinDE) * dR0 + (0.5 * R0 * sinDE / sqrta) * da; // sqrtmu G = sigma0 r0 (1-F) + r0 sqrta sinDE
        dG = dG / sqrtmu;
        dFt = (-sqrta / R0 / Rf * cosDE) * dDE - (0.5 / sqrta / R0 / Rf * sinDE) * da
              + (sqrta / Rf / R0 / R0 * sinDE) * dR0 + (sqrta / Rf / Rf / R0 * sinDE) * dRf;
        dFt = dFt * sqrtmu;
        dGt = -(1 - cosDE) / Rf * da + a / Rf / Rf * (1 - cosDE) * dRf - a / Rf * sinDE * dDE;
    } else { // hyperbolas (sqrta is sqrt(-a))
        double sqrta = std::sqrt(-a);
        double sinhDH = std::sinh(DE);
        double coshDH = std::cosh(DE);

        mat16 ds0 = dsigma0 / sqrta + 0.5 * sigma0 / sqrta / sqrta / sqrta * da; // s0 = sigma0 / sqrta
        mat16 dc0 = -1. / a * dR0 + R0 / a / a * da;                             // c0 = (1- R/a)
        mat16 dDN = 1.5 * sqrtmu * tof / std::pow(sqrta, 5) * da;                // N = sqrt(-mu/a**3) tof
        mat16 dDH = (dDN - (coshDH - 1) * ds0 - sinhDH * dc0) / (s0 * sinhDH + c0 * coshDH - 1);
        mat16 dRf = (coshDH - 1 - 0.5 / sqrta * sigma0 * sinhDH) * da + coshDH * dR0
                    + (sigma0 * sqrta * coshDH + (R0 + a) * sinhDH) * dDH
                    + sqrta * sinhDH * dsigma0; // r = -a + (r0 + a) * coshDH + sigma0 * sqrta * sinhDH

        // 2 - We may now compute the differentials of the Lagrange coefficients
        dF = -(1 - coshDH) / R0 * da + a / R0 / R0 * (1 - coshDH) * dR0 + a / R0 * sinhDH * dDH;
        dG = (1 - F) * (R0 * dsigma0 + sigma0 * dR0) - (sigma0 * R0) * dF + (sqrta * R0 * coshDH) * dDH
             + (sqrta * sinhDH) * dR0
             - (0.5 * R0 * sinhDH / sqrta) * da; // sqrtmu G = sigma0 r0 (1-F) + r0 sqrta sinhDH
        dG = dG / sqrtmu;
        dFt = (-sqrta / R0 / Rf * coshDH) * dDH + (0.5 / sqrta / R0 / Rf * sinhDH) * da
              + (sqrta / Rf / R0 / R0 * sinhDH) * dR0 + (sqrta / Rf / Rf / R0 * sinhDH) * dRf;
        dFt = dFt * sqrtmu;
        dGt = (1 - coshDH) / Rf * da + a / Rf / Rf * (1 - coshDH) * dRf + a / Rf * sinhDH * dDH;
    }
    // 3 - And finally assemble the state transition matrix
    mat36 Mr = F * dr0 + _dot(mat31(xt::transpose(r0)), dF) + G * dv0 + _dot(mat31(xt::transpose(v0)), dG);
    mat36 Mv = Ft * dr0 + _dot(mat31(xt::transpose(r0)), dFt) + Gt * dv0 + _dot(mat31(xt::transpose(v0)), dGt);
    mat66 M{};
    xt::view(M, xt::range(0, 3), xt::all()) = Mr;
    xt::view(M, xt::range(3, 6), xt::all()) = Mv;
    // ... and flatten it
    std::array<double, 36> retval{};
    std::copy(M.begin(), M.end(), retval.begin());
    return retval;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
mat66 _compute_Y(const mat31 &r0, const mat31 &v0, const mat31 &r, const mat31 &v, double tof, double mu)
{
    mat31 h = _cross(r, v);
    double r0_mod = std::sqrt(r0(0, 0) * r0(0, 0) + r0(1, 0) * r0(1, 0) + r0(2, 0) * r0(2, 0));
    double r_mod = std::sqrt(r(0, 0) * r(0, 0) + r(1, 0) * r(1, 0) + r(2, 0) * r(2, 0));
    double r3 = r_mod * r_mod * r_mod;
    mat32 B{};
    xt::view(B, xt::all(), 0) = xt::view(r0 / std::sqrt(mu * r0_mod), xt::all(), 0);
    xt::view(B, xt::all(), 1) = xt::view(v0 * r0_mod / mu, xt::all(), 0);
    mat63 fc{};
    xt::view(fc, xt::range(0, 3), xt::all()) = _skew(r);
    xt::view(fc, xt::range(3, 6), xt::all()) = _skew(v);
    auto sct = -_dot(xt::eval(_dot(_skew(r), _skew(v)) + _skew(h)), B);
    auto scb = _dot(xt::eval(mu / r3 * _dot(_skew(r), _skew(r)) - _dot(_skew(v), _skew(v))), B);
    mat62 sc{};
    xt::view(sc, xt::range(0, 3), xt::all()) = sct;
    xt::view(sc, xt::range(3, 6), xt::all()) = scb;
    auto tct = (-r + 1.5 * v * tof);
    auto tcb = (v / 2. - 1.5 * mu / r3 * r * tof);
    mat61 tc{};
    xt::view(tc, xt::range(0, 3), xt::all()) = tct;
    xt::view(tc, xt::range(3, 6), xt::all()) = tcb;
    mat66 Y{};
    xt::view(Y, xt::all(), xt::range(0, 3)) = fc;
    xt::view(Y, xt::all(), xt::range(3, 5)) = sc;
    xt::view(Y, xt::all(), xt::range(5, 6)) = tc;
    return Y;
}

// From:
// Reynolds, Reid G. "Direct Solution of the Keplerian State Transition Matrix." Journal of Guidance, Control, and
// Dynamics 45, no. 6 (2022): 1162-1165.
std::array<double, 36> stm(const std::array<std::array<double, 3>, 2> &pos_vel0,
                           const std::array<std::array<double, 3>, 2> &pos_velf, double tof, double mu)
{
    mat31 r0 = {{pos_vel0[0][0], pos_vel0[0][1], pos_vel0[0][2]}};
    mat31 v0 = {{pos_vel0[1][0], pos_vel0[1][1], pos_vel0[1][2]}};
    mat31 rf = {{pos_velf[0][0], pos_velf[0][1], pos_velf[0][2]}};
    mat31 vf = {{pos_velf[1][0], pos_velf[1][1], pos_velf[1][2]}};
    // And the output (6,6)
    // std::array<double, 36> retval{};
    mat66 retval{};
    // auto retval_xt = xt::adapt(retval, shape66);
    //  Compute the STM using Reynolds' Cartesian Representation
    auto Y = _compute_Y(r0, v0, rf, vf, tof, mu);
    auto Y0 = _compute_Y(r0, v0, r0, v0, 0., mu);
    retval = dot(Y, inv(Y0));
    // return retval;
    std::array<double, 36> ret{};
    std::copy(retval.begin(), retval.end(), ret.begin());
    return ret;
}

std::pair<std::array<std::array<double, 3>, 2>, std::array<double, 36>>
propagate_stm(const std::array<std::array<double, 3>, 2> &pos_vel0, double tof, double mu)
{
    auto pos_velf = pos_vel0;
    kep3::propagate_lagrangian(pos_velf, tof, mu);
    auto retval_stm = stm(pos_vel0, pos_velf, tof, mu);
    return {pos_velf, retval_stm};
}

std::pair<std::array<std::array<double, 3>, 2>, std::array<double, 36>>
propagate_stm2(const std::array<std::array<double, 3>, 2> &pos_vel0, double tof, double mu)
{
    auto pos_velf = pos_vel0;
    kep3::propagate_lagrangian(pos_velf, tof, mu);
    const auto &[r0, v0] = pos_vel0;
    double R = std::sqrt(r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2]);
    double V2 = v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2];
    double energy = (V2 / 2 - mu / R);
    double a = -mu / 2.0 / energy; // will be negative for hyperbolae
    double sqrta = 0.;
    double F = 0., G = 0., Ft = 0., Gt = 0.;
    double sigma0 = (r0[0] * v0[0] + r0[1] * v0[1] + r0[2] * v0[2]) / std::sqrt(mu);

    sqrta = std::sqrt(a);
    double DM = std::sqrt(mu / std::pow(a, 3)) * tof;
    double sinDM = std::sin(DM), cosDM = std::cos(DM);
    // Here we use the atan2 to recover the mean anomaly difference in the
    // [0,2pi] range. This makes sure that for high value of M no catastrophic
    // cancellation occurs, as would be the case using std::fmod(DM, 2pi)
    double DM_cropped = std::atan2(sinDM, cosDM);
    if (DM_cropped < 0) {
        DM_cropped += 2 * kep3::pi;
    }
    double s0 = sigma0 / sqrta;
    double c0 = (1 - R / a);
    // This initial guess was developed applying Lagrange expansion theorem to
    // the Kepler's equation in DE. We stopped at 3rd order.
    double IG = DM_cropped + c0 * sinDM - s0 * (1 - cosDM) + (c0 * cosDM - s0 * sinDM) * (c0 * sinDM + s0 * cosDM - s0)
                + 0.5 * (c0 * sinDM + s0 * cosDM - s0)
                      * (2 * std::pow(c0 * cosDM - s0 * sinDM, 2)
                         - (c0 * sinDM + s0 * cosDM - s0) * (c0 * sinDM + s0 * cosDM));

    // Solve Kepler Equation for ellipses in DE (eccentric anomaly difference)
    const int digits = std::numeric_limits<double>::digits;
    std::uintmax_t max_iter = 100u;
    // NOTE: Halley iterates may result into instabilities (specially with a
    // poor IG)

    double DE = boost::math::tools::newton_raphson_iterate(
        [DM_cropped, sigma0, sqrta, a, R](double DE) {
            return std::make_tuple(kepDE(DE, DM_cropped, sigma0, sqrta, a, R), d_kepDE(DE, sigma0, sqrta, a, R));
        },
        IG, IG - pi, IG + pi, digits, max_iter);
    if (max_iter == 100u) {
        throw std::domain_error(fmt::format("Maximum number of iterations exceeded when solving Kepler's "
                                            "equation for the eccentric anomaly in propagate_lagrangian.\n"
                                            "DM={}\nsigma0={}\nsqrta={}\na={}\nR={}\nDE={}",
                                            DM, sigma0, sqrta, a, R, DE));
    }
    double r = a + (R - a) * std::cos(DE) + sigma0 * sqrta * std::sin(DE);

    // Lagrange coefficients
    F = 1 - a / R * (1 - std::cos(DE));
    G = a * sigma0 / std::sqrt(mu) * (1 - std::cos(DE)) + R * std::sqrt(a / mu) * std::sin(DE);
    Ft = -std::sqrt(mu * a) / (r * R) * std::sin(DE);
    Gt = 1 - a / r * (1 - std::cos(DE));

    auto retval_stm = stm_l(pos_vel0, pos_velf, tof, mu, R, V2, energy, sigma0, a, s0, c0, DE, F, G, Ft, Gt);
    return {pos_velf, retval_stm};
}
} // namespace kep3