// Copyright 2023 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the cascade library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include<iostream>

#include<core_astro/convert_anomalies.hpp>

#include "catch.hpp"

using namespace kep3;

TEST_CASE("m2e")

{
    using Catch::Detail::Approx;
    std::random_device rd;
    //
    // Engines 
    //
    std::mt19937 rng_engine(rd());
    //
    // Distribtuions
    //
    std::uniform_real_distribution<> ecc_difficult_d(0.9, 0.99);
    std::uniform_real_distribution<> ecc_easy_d(0., 0.9);
    std::uniform_real_distribution<> M_d(0., 1.);

    // Testing on N random calls
    unsigned N = 10000;
    for (auto i = 0u; i < 10000; ++i) {
        auto mean_anom = M_d(rng_engine);
        auto ecc = ecc_easy_d(rng_engine);
        std::cout << e2m(m2e(mean_anom, ecc), ecc) - mean_anom << std::endl;
        REQUIRE(e2m(m2e(mean_anom, ecc), ecc) == Approx(mean_anom).epsilon(0.).margin(1e-11));
    }
}