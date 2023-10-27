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
using mat33 = xt::xtensor_fixed<double, xt::xshape<3, 3>>;
using mat66 = xt::xtensor_fixed<double, xt::xshape<6, 6>>;
using mat32 = xt::xtensor_fixed<double, xt::xshape<3, 2>>;
using mat62 = xt::xtensor_fixed<double, xt::xshape<6, 2>>;
using mat61 = xt::xtensor_fixed<double, xt::xshape<6, 1>>;
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
                    C(i,j) += A(i,l) * B(l,j);
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
    // Shapes needed
    std::vector<std::size_t> shape66 = {6, 6};
    std::vector<std::size_t> shapev = {3, 1};
    // Lets create an xtensor interface for the inputs (3,1)
    auto r0 = xt::adapt(pos_vel0[0], shapev);
    auto v0 = xt::adapt(pos_vel0[1], shapev);
    auto rf = xt::adapt(pos_velf[0], shapev);
    auto vf = xt::adapt(pos_velf[1], shapev);
    // And the output (6,6)
    std::array<double, 36> retval{};
    auto retval_xt = xt::adapt(retval, shape66);
    // Compute the STM using Reynolds' Cartesian Representation
    auto Y = _compute_Y(r0, v0, rf, vf, tof, mu);
    auto Y0 = _compute_Y(r0, v0, r0, v0, 0., mu);
    retval_xt = _dot(Y, mat66(inv(Y0)));
    return retval;
}

std::pair<std::array<std::array<double, 3>, 2>, std::array<double, 36>>
propagate_stm(const std::array<std::array<double, 3>, 2> &pos_vel0, double tof, double mu)
{
    auto pos_velf = pos_vel0;
    kep3::propagate_lagrangian(pos_velf, tof, mu);
    auto retval_stm = stm(pos_vel0, pos_velf, tof, mu);
    return {pos_velf, retval_stm};
}

} // namespace kep3