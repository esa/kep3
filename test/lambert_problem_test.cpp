// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>

#include <fmt/core.h>
#include <fmt/ranges.h>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xadapt.hpp>

#include <kep3/lambert_problem.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

TEST_CASE("construct") {
  // Here we test construction of a few simple geometries
  kep3::lambert_problem lp{{1., 0., 0.}, {0., 1., 0.}, 3 * kep3::pi / 2,
                           1.,           true,         100};
}

TEST_CASE("delta_guidance") {
  // Here we test that in a number of randomly generated Lambert Problems
  // The boundary data must satisfy the Delta Guidance law

  // Preamble
  std::array<double, 3> r1{{0, 0, 0}}, r2{{0, 0, 0}};
  double tof = 0.;

  // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
  std::mt19937 rng_engine(12201203u);
  std::uniform_int_distribution<unsigned> cw_d(0, 1);
  std::uniform_real_distribution<double> r_d(-2, 2);
  std::uniform_real_distribution<double> tof_d(2., 40.);
  std::uniform_real_distribution<double> mu_d(0.9, 1.1);
  unsigned revs_max = 20u;

  unsigned trials = 10000u;

  for (auto i = 0u; i < trials; ++i) {
    // 1 - generate a random problem geometry
    r1[0] = r_d(rng_engine);
    r1[1] = r_d(rng_engine);
    r1[2] = r_d(rng_engine);
    r2[0] = r_d(rng_engine);
    r2[1] = r_d(rng_engine);
    r2[2] = r_d(rng_engine);
    tof = tof_d(rng_engine);
    bool cw = static_cast<bool>(cw_d(rng_engine));
    double mu = mu_d(rng_engine);

    // 2 - Solve the lambert problem
    kep3::lambert_problem lp(r1, r2, tof, mu, cw, revs_max);

    // 3 - Check the Delta guidance error
    for (const auto &v1 : lp.get_v1()) {
      double dg_err = kep3_tests::delta_guidance_error(r1, r2, v1, mu);
      if (!(dg_err < 1e-12)) {
        std::cout << lp << std::endl;
        fmt::print("\nr1= {}\nr2= {}\ntof= {}\nmu= {}\ncw= {}\nrevs_max= {}", r1, r2,tof, mu, cw, revs_max);
      }
      REQUIRE(dg_err < 1e-12);
    }
  }
}