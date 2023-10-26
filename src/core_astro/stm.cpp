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

namespace kep3
{

xt::xarray<double> _skew(const xt::xarray<double> &v)
{
    xt::xarray<double> retval{{0., -v(2, 0), v(1, 0)}, {v(2, 0), 0., -v(0, 0)}, {-v(1, 0), v(0, 0), 0.}};
    return retval;
}

xt::xarray<double> _cross(const xt::xarray<double> &v1, const xt::xarray<double> &v2)
{
    xt::xarray<double> retval{{v1(1, 0) * v2(2, 0) - v2(1, 0) * v1(2, 0), v2(0, 0) * v1(2, 0) - v1(0, 0) * v2(2, 0),
                               v1(0, 0) * v2(1, 0) - v2(0, 0) * v1(1, 0)}};
    return retval.reshape({-1, 1});
}

xt::xarray<double> _compute_Y(const xt::xarray<double> &r0, const xt::xarray<double> &v0, const xt::xarray<double> &r,
                              const xt::xarray<double> &v, double tof, double mu)
{
    auto h = _cross(r, v);
    auto r0_mod = xt::linalg::norm(r0);
    auto r3 = std::pow(xt::linalg::norm(r), 3);
    auto B = xt::concatenate(xt::xtuple(r0 / std::sqrt(mu * r0_mod), v0 * r0_mod / mu), 1);
    auto fc = xt::concatenate(xt::xtuple(_skew(r), _skew(v)));
    auto tmp = dot(_skew(r), _skew(v));
    auto sct = -dot((tmp + _skew(h)), B);
    auto scb = dot((mu / r3 * dot(_skew(r), _skew(r)) - dot(_skew(v), _skew(v))), B);
    auto sc = xt::concatenate(xt::xtuple(sct, scb));
    auto tct = (-r + 1.5 * v * tof);
    auto tcb = (v / 2. - 1.5 * mu / r3 * r * tof);
    auto tc = xt::concatenate(xt::xtuple(tct, tcb));
    auto Y = xt::concatenate(xt::xtuple(fc, sc, tc), 1);
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
    retval_xt = dot(Y, inv(Y0));
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