// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>

#include <kep3/lambert_problem.hpp>

#include "catch.hpp"


TEST_CASE("construct") {
  kep3::lambert_problem lp{{1.,0.,0.}, {0.,1.,0.}, kep3::pi/2 + 2*kep3::pi,1.,false,100};
  std::cout << lp << std::endl;

}