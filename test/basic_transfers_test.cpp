// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <kep3/core_astro/basic_transfers.hpp>

#include "catch.hpp"

TEST_CASE("hohmann")
{
    // Test cases ground tryths from (https://www.omnicalculator.com/physics/hohmann-transfer)
    // thus the low precision requested
    {
        auto h_res = kep3::hohmann(1., 2., 1.);
        REQUIRE(h_res.first == Approx(0.28445705).epsilon(1e-6));
        REQUIRE(h_res.second[0] == Approx(0.15470054).epsilon(1e-6));
        REQUIRE(h_res.second[1] == Approx(0.1297565).epsilon(1e-6));
    }
    {
        auto h_res = kep3::hohmann(1.1, 2.2, 1.3);
        REQUIRE(h_res.first == Approx(0.3092374).epsilon(1e-6));
        REQUIRE(h_res.second[0] == Approx(0.1681772).epsilon(1e-6));
        REQUIRE(h_res.second[1] == Approx(0.1410602).epsilon(1e-6));
    }
}

TEST_CASE("bielliptic") {
    // Test cases ground truths from (https://www.omnicalculator.com/physics/hohmann-transfer)
    // thus the low precision requested. Playing with the fact that a bielliptic transfer collpses 
    // to Hohmann if rb=r2
    {
        auto h_res = kep3::bielliptic(1., 2., 2., 1.);
        REQUIRE(h_res.first == Approx(0.28445705).epsilon(1e-6));
        REQUIRE(h_res.second[0] == Approx(0.15470054).epsilon(1e-6));
        REQUIRE(h_res.second[1] == Approx(0.1297565).epsilon(1e-6));
    }
    {
        auto h_res = kep3::bielliptic(1.1, 2.2, 2.2, 1.3);
        REQUIRE(h_res.first == Approx(0.3092374).epsilon(1e-6));
        REQUIRE(h_res.second[0] == Approx(0.1681772).epsilon(1e-6));
        REQUIRE(h_res.second[1] == Approx(0.1410602).epsilon(1e-6));
    }
}