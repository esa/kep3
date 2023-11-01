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
#include <stdexcept>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/core_astro/ic2eq2ic.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/exceptions.hpp>
#include <kep3/planets/jpl_lp.hpp>
#include <kep3/planets/vsop2013.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using kep3::udpla::jpl_lp;
using kep3::udpla::vsop2013;

TEST_CASE("constructor")
{
    vsop2013 pl1("mercury", 1e-5);
    fmt::println("{}", pl1.eph(0));

    jpl_lp pl2("mercury");
    fmt::println("{}", pl2.eph(0));
}
