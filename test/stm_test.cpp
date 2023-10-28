// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/core_astro/stm.hpp>
#include <kep3/planets/jpl_lp.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using xt::linalg::dot;

TEST_CASE("stm")
{
    // Test that the STM is the identity if tof = 0
    {
        std::array<std::array<double, 3>, 2> pos_vel0 = {{{1., 0., 0.}, {0., 1., 0.}}};
        std::array<std::array<double, 3>, 2> pos_velf = {{{1, 0., 0.}, {0., 1., 0.}}};
        double tof = 0;
        double mu = 1.;
        auto computed = kep3::stm(pos_vel0, pos_velf, tof, mu);
        std::array<double, 36> real{1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
        REQUIRE(kep3_tests::L_infinity_norm(computed, real) < 1e-16);
    }
    // Test a random fake case (ground truth obtained from the python code equivalent)
    {
        std::array<std::array<double, 3>, 2> pos_vel0 = {{{1.1, 0.1, 0.1}, {0.1, 1.1, 0.1}}};
        std::array<std::array<double, 3>, 2> pos_velf = {{{-1, 0.2, 0.2}, {0., -1., 0.}}};
        double tof = 1.34;
        double mu = 2.34;
        auto computed = kep3::stm(pos_vel0, pos_velf, tof, mu);
        std::array<double, 36> real{-1.5691824598131583, 0.1537355920717801,  -0.2846205723117815, 0.3538015104605221,
                                    -1.3684573575369972, -0.0845546539230396, 2.8449067152872787,  2.7467972480896288,
                                    0.5659753302814089,  3.572868357411019,   1.9316889178225798,  0.1920464396027999,
                                    0.4105288906043572,  0.2266109701296614,  -0.8802383449388316, 0.2157344694942041,
                                    0.2908873836143274,  0.308885154425711,   -5.209787769406611,  -1.9865454950329102,
                                    -0.6830277720366267, -2.6952058580669065, -3.0050517627805746, -0.3916881350706234,
                                    2.417990336024268,   0.7183181583660175,  0.2613590411991904,  0.7249455411331721,
                                    2.490891546462969,   0.2679864239663451,  1.0473579218305524,  0.4380489600644716,
                                    0.2071172401579186,  0.4088107529497264,  0.7257376435683549,  -0.8221209669568265};
        REQUIRE(kep3_tests::L_infinity_norm(computed, real) < 1e-14);
    }
}

TEST_CASE("propagate_stm")
{
    // We test the identity stm02 = stm12stm01
    std::array<std::array<double, 3>, 2> pos_vel = {{{1.2, -0.3, 0.}, {0., 2.3, -1.2}}};
    double mu = 0.89732;
    double dt1 = 1.1;
    double dt2 = 2.3;
    auto res01 = kep3::propagate_stm(pos_vel, dt1, mu);
    auto res02 = kep3::propagate_stm(pos_vel, dt2, mu);
    // We propagate in-place pos_vel ...
    kep3::propagate_lagrangian(pos_vel, dt1, mu);
    // ... and compute the stm from 1 to 2.
    auto res12 = kep3::propagate_stm(pos_vel, dt2 - dt1, mu);
    auto M01 = xt::adapt(res01.second, {6, 6});
    auto M02 = xt::adapt(res02.second, {6, 6});
    auto M12 = xt::adapt(res12.second, {6, 6});
    REQUIRE(xt::linalg::norm(M02 - dot(M12, M01)) < 1e-13);
}

TEST_CASE("propagate_stm2")
{
    // We test the identity stm02 = stm12stm01
    std::array<std::array<double, 3>, 2> pos_vel = {{{1., 0., 0.}, {0., 1.0, 0.}}};
    double mu = 1;
    double dt = 2*kep3::pi;
    auto res2 = kep3::propagate_stm2(pos_vel, dt, mu);
    fmt::print("{}", res2);
    auto res = kep3::propagate_stm(pos_vel, dt, mu);
    fmt::print("\n\n{}", res);

}