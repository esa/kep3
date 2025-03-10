// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <array> // For std::array
#include <cmath>
#include <utility> // For std::pair

#include <kep3/core_astro/basic_transfers.hpp>

namespace kep3
{
std::pair<double, std::array<double, 2>> hohmann(double r1, double r2, double mu)
{
    double v1 = std::sqrt(mu / r1);
    double v2 = std::sqrt(mu / r2);
    double vt1 = std::sqrt(mu / r1 * (2 * r2 / (r1 + r2)));
    double vt2 = std::sqrt(mu / r2 * (2 * r1 / (r1 + r2)));
    double dv1 = vt1 - v1;
    double dv2 = v2 - vt2;
    double dv_total = dv1 + dv2;
    return {dv_total, {dv1, dv2}};
}

std::pair<double, std::array<double, 3>> bielliptic(double r1, double r2, double rb, double mu)
{
    double v1 = std::sqrt(mu / r1);
    double v2 = std::sqrt(mu / r2);

    double vt1 = std::sqrt(mu / r1 * (2 * rb / (r1 + rb)));
    double vt1a = std::sqrt(mu / rb * (2 * r1 / (r1 + rb)));

    double vt2 = std::sqrt(mu / rb * (2 * r2 / (rb + r2)));
    double vt2a = std::sqrt(mu / r2 * (2 * rb / (rb + r2)));

    double dv1 = vt1 - v1;
    double dv2 = vt2 - vt1a;
    double dv3 = v2 - vt2a;
    double dv_total = dv1 + dv2 + dv3;

    return {dv_total, {dv1, dv2, dv3}};
}
} // namespace kep3