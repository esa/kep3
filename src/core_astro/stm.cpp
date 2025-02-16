// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <array>
#include <cmath>

#include <boost/math/tools/roots.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xadapt.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/kepler_equations.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/core_astro/stm.hpp>
#include <kep3/linalg.hpp>

namespace kep3
{

using xt::linalg::inv;
using kep3::linalg::mat36;
using kep3::linalg::mat16;
using kep3::linalg::mat31;
using kep3::linalg::mat66;
using kep3::linalg::mat63;
using kep3::linalg::mat32;
using kep3::linalg::mat62;
using kep3::linalg::mat61;
using kep3::linalg::mat13;
using kep3::linalg::_dot;
using kep3::linalg::_cross;
using kep3::linalg::_skew;

// Here we take the lagrangian coefficient expressions for rf and vf as function of r0 and v0, and manually,
// differentiate it to obtain the state transition matrix.
std::array<double, 36> stm_lagrangian(const std::array<std::array<double, 3>, 2> &pos_vel0, double tof, // NOLINT
                                      double mu,                                                        // NOLINT
                                      double R0, double Rf, double energy,                              // NOLINT
                                      double sigma0,                                                    // NOLINT
                                      double a, double s0, double c0,                                   // NOLINT
                                      double DX, double F, double G, double Ft, double Gt)
{
    // Create xtensor fixed arrays from input (we avoid adapt as its slower in this case since all is fixed size)
    // We use row vectors (not column) as its then more conventional for gradients as the differential goes to the
    // end
    mat13 r0 = {{pos_vel0[0][0]}, {pos_vel0[0][1]}, {pos_vel0[0][2]}};
    mat13 v0 = {{pos_vel0[1][0]}, {pos_vel0[1][1]}, {pos_vel0[1][2]}};

    // We seed the gradients with the initial dr0/dx0 and dv0/dx0
    mat36 dr0 = xt::zeros<double>({3, 6});
    mat36 dv0 = xt::zeros<double>({3, 6});
    dr0(0, 0) = 1;
    dr0(1, 1) = 1;
    dr0(2, 2) = 1;
    dv0(0, 3) = 1;
    dv0(1, 4) = 1;
    dv0(2, 5) = 1;

    // 1 - We start computing the differentials of basic quantities. A differential for a scalar will be a 16 mat
    // (gradient).
    double sqrtmu = std::sqrt(mu);
    mat16 dV02 = 2. * _dot(v0, dv0);
    mat16 dR0 = 1. / R0 * _dot(r0, dr0);
    mat16 denergy = 0.5 * dV02 + mu / R0 / R0 * dR0;
    mat16 dsigma0 = ((_dot(r0, dv0) + _dot(v0, dr0))) / sqrtmu;
    mat16 da = mu / 2. / energy / energy * denergy; // a = -mu / 2 / energy
    mat16 dF, dFt, dG, dGt;

    if (a > 0) { // ellipses
        double sqrta = std::sqrt(a);
        double sinDE = std::sin(DX);
        double cosDE = std::cos(DX);

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
        double sinhDH = std::sinh(DX); // DX is here the hyperbolic anomaly.
        double coshDH = std::cosh(DX);

        mat16 ds0 = dsigma0 / sqrta + 0.5 * sigma0 / sqrta / sqrta / sqrta * da; // s0 = sigma0 / sqrta
        mat16 dc0 = -1. / a * dR0 + R0 / a / a * da;                             // c0 = (1- R/a)
        mat16 dDN = 1.5 * sqrtmu * tof / std::pow(sqrta, 5) * da;                // N = sqrt(-mu/a**3) tof
        mat16 dDH = (dDN - (coshDH - 1) * ds0 - sinhDH * dc0) / (s0 * sinhDH + c0 * coshDH - 1);
        mat16 dRf = (1 - coshDH - 0.5 / sqrta * sigma0 * sinhDH) * da + coshDH * dR0
                    + (sigma0 * sqrta * coshDH + (R0 - a) * sinhDH) * dDH
                    + sqrta * sinhDH * dsigma0; // r = a + (r0 - a) * coshDH + sigma0 * sqrta * sinhDH

        // 2 - We may now compute the differentials of the Lagrange coefficients
        dF = -(1 - coshDH) / R0 * da + a / R0 / R0 * (1 - coshDH) * dR0 + a / R0 * sinhDH * dDH;
        dG = (1 - F) * (R0 * dsigma0 + sigma0 * dR0) - (sigma0 * R0) * dF + (sqrta * R0 * coshDH) * dDH
             + (sqrta * sinhDH) * dR0
             - (0.5 * R0 * sinhDH / sqrta) * da; // sqrtmu G = sigma0 r0 (1-F) + r0 sqrta sinhDH
        dG = dG / sqrtmu;
        dFt = (-sqrta / R0 / Rf * coshDH) * dDH + (0.5 / sqrta / R0 / Rf * sinhDH) * da
              + (sqrta / Rf / R0 / R0 * sinhDH) * dR0 + (sqrta / Rf / Rf / R0 * sinhDH) * dRf;
        dFt = dFt * sqrtmu;
        dGt = -(1 - coshDH) / Rf * da + a / Rf / Rf * (1 - coshDH) * dRf + a / Rf * sinhDH * dDH;
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
std::array<double, 36> stm_reynolds(const std::array<std::array<double, 3>, 2> &pos_vel0,
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
    retval = _dot(Y, mat66(inv(Y0)));
    // return retval;
    std::array<double, 36> ret{};
    std::copy(retval.begin(), retval.end(), ret.begin());
    return ret;
}

std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>
propagate_stm_reynolds(const std::array<std::array<double, 3>, 2> &pos_vel0, double tof, double mu, bool stm)
{
    auto res = kep3::propagate_lagrangian(pos_vel0, tof, mu);
    if (stm) {
        // NOLINTNEXTLINE(readability-suspicious-call-argument)
        auto retval_stm = stm_reynolds(pos_vel0, res.first, tof, mu);
        return {res.first, retval_stm};
    } else {
        return {res.first, std::nullopt};
    }
}

} // namespace kep3