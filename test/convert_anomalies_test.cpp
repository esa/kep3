// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include <random>

#include <kep3/core_astro/convert_anomalies.hpp>

#include "catch.hpp"

using kep3::e2m;
using kep3::m2e;
using kep3::e2f;
using kep3::f2e;
using kep3::zeta2f;
using kep3::f2zeta;


TEST_CASE("m2e") {
  using Catch::Detail::Approx;
  std::random_device rd;
  //
  // Engines
  //
  std::mt19937 rng_engine(rd());
  //
  // Distributions
  //
  std::uniform_real_distribution<double> ecc_difficult_d(0.9, 0.99);
  std::uniform_real_distribution<double> ecc_easy_d(0., 0.9);
  std::uniform_real_distribution<double> M_d(-100, 100.);

  // Testing on N random calls (easy)
  unsigned N = 10000;
  for (auto i = 0u; i < N; ++i) {
    double mean_anom = M_d(rng_engine);
    double ecc = ecc_easy_d(rng_engine);
    double res = e2m(m2e(mean_anom, ecc), ecc);
    // We require 1e-13 and not 1e-14 as the range for m is [-100,100]. For high
    // values of M some precision loss is to be expected as to maintain digits
    // of the higher number.
    REQUIRE(std::sin(res) ==
            Approx(std::sin(mean_anom)).epsilon(0.).margin(1e-13));
    REQUIRE(std::cos(res) ==
            Approx(std::cos(mean_anom)).epsilon(0.).margin(1e-13));
  }
  // Testing on N random calls (difficult)
  for (auto i = 0u; i < N; ++i) {
    double mean_anom = M_d(rng_engine);
    double ecc = ecc_difficult_d(rng_engine);
    double res = e2m(m2e(mean_anom, ecc), ecc);
    REQUIRE(std::sin(res) ==
            Approx(std::sin(mean_anom)).epsilon(0.).margin(1e-13));
    REQUIRE(std::cos(res) ==
            Approx(std::cos(mean_anom)).epsilon(0.).margin(1e-13));
  }
}

TEST_CASE("f2e") {
  using Catch::Detail::Approx;
  std::random_device rd;
  //
  // Engines
  //
  std::mt19937 rng_engine(rd());
  //
  // Distributions
  //
  std::uniform_real_distribution<double> ecc_difficult_d(0.9, 0.99);
  std::uniform_real_distribution<double> ecc_easy_d(0., 0.9);
  std::uniform_real_distribution<double> f_d(-100, 100.);

  // Testing on N random calls (easy)
  unsigned N = 10000;
  for (auto i = 0u; i < N; ++i) {
    double true_anom = f_d(rng_engine);
    double ecc = ecc_easy_d(rng_engine);
    double res = e2f(f2e(true_anom, ecc), ecc);
    REQUIRE(std::sin(res) ==
            Approx(std::sin(true_anom)).epsilon(0.).margin(1e-14));
    REQUIRE(std::cos(res) ==
            Approx(std::cos(true_anom)).epsilon(0.).margin(1e-14));
  }
  // Testing on N random calls (difficult)
  for (auto i = 0u; i < N; ++i) {
    double true_anom = f_d(rng_engine);
    double ecc = ecc_difficult_d(rng_engine);
    double res = e2f(f2e(true_anom, ecc), ecc);
    REQUIRE(std::sin(res) ==
            Approx(std::sin(true_anom)).epsilon(0.).margin(1e-14));
    REQUIRE(std::cos(res) ==
            Approx(std::cos(true_anom)).epsilon(0.).margin(1e-14));
  }
}

TEST_CASE("zeta2e") {
  using Catch::Detail::Approx;
  std::random_device rd;
  //
  // Engines
  //
  std::mt19937 rng_engine(rd());
  //
  // Distributions
  //
  std::uniform_real_distribution<double> ecc_difficult_d(1.01, 1.1);
  std::uniform_real_distribution<double> ecc_easy_d(2., 100.);
  std::uniform_real_distribution<double> f_d(-100, 100.);

  // Testing on N random calls (easy)
  unsigned N = 10000;
  for (auto i = 0u; i < N; ++i) {
    double true_anom = f_d(rng_engine);
    double ecc = ecc_easy_d(rng_engine);
    double res = zeta2f(f2zeta(true_anom, ecc), ecc);
    REQUIRE(std::sin(res) ==
            Approx(std::sin(true_anom)).epsilon(0.).margin(1e-14));
    REQUIRE(std::cos(res) ==
            Approx(std::cos(true_anom)).epsilon(0.).margin(1e-14));
  }
  // Testing on N random calls (difficult)
  for (auto i = 0u; i < N; ++i) {
    double true_anom = f_d(rng_engine);
    double ecc = ecc_difficult_d(rng_engine);
    double res = zeta2f(f2zeta(true_anom, ecc), ecc);
    REQUIRE(std::sin(res) ==
            Approx(std::sin(true_anom)).epsilon(0.).margin(1e-14));
    REQUIRE(std::cos(res) ==
            Approx(std::cos(true_anom)).epsilon(0.).margin(1e-14));
  }
}
