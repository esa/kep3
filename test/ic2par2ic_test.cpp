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

#include <kep3/core_astro/ic2par.hpp>
#include <kep3/core_astro/par2ic.hpp>

#include "catch.hpp"

TEST_CASE("ic2par") {
    fmt::print("\nResult: {}", kep3::ic2par({1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0));
    fmt::print("\nResult: {}", kep3::par2ic({1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 1.0));
    fmt::print("\nResult: {}", kep3::ic2par({1.0, 0.0, 0.0}, {0.0, 0.0, 1.1}, 1.0));
    fmt::print("\nResult: {}", kep3::par2ic({1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 1.0));

}
