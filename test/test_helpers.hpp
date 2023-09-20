// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_TEST_HELPERS_H
#define kep3_TEST_HELPERS_H

#include <algorithm>
#include <array>
#include <cmath>

namespace kep3_tests {
// This is a float test which, while controversial, will test for abs differences in small numbers, 
// relative otherwise.
inline double floating_point_error(double a, double b) {
  return std::abs(a - b) / std::max(1., std::max(a, b));
}

// This tests how close two vectors are in the euclidean metric. err = r2-r1
inline double floating_point_error_vector(const std::array<double, 3> &r1, const std::array<double, 3> &r2) {
  double R1 = std::sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);
  std::array<double, 3> r12 = {{r2[0]-r1[0], r2[1]-r1[1], r2[2]-r1[2]}};
  double R12 = std::sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
  return R12 / std::max(1., R1);
}
} // namespace kep3_tests

#endif // kep3_TEST_HELPERS_H
