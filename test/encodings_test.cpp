// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <kep3/core_astro/encodings.hpp>

#include "catch.hpp"

TEST_CASE("encodings")
{
    // We only call these here: their correctness is tested in python
    kep3::alpha2direct({0.1,0.2,0.3},  21.2);
    kep3::direct2alpha({0.1,0.2,0.3});
    kep3::eta2direct({0.1,0.2,0.3},  21.2);
    kep3::direct2eta({0.1,0.2,0.3},  21.2);
}