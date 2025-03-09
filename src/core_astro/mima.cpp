// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <cmath>
#include <kep3/core_astro/mima.hpp>

namespace kep3
{

std::pair<double, double> mima(std::array<double, 3> dv1, std::array<double, 3> dv2, double tof, double Tmax, // NOLINT
                               double veff)                                                                   // NOLINT
{
    std::array<double, 3> dv{dv1[0]+dv2[0], dv1[1]+dv2[1], dv1[2]+dv2[2]};
    std::array<double, 3> dv_diff{-dv1[0]+dv2[0], -dv1[1]+dv2[1], -dv1[2]+dv2[2]};
    double ab = dv[0]*dv_diff[0] + dv[1]*dv_diff[1]+dv[2]*dv_diff[2]; // ab = dv@dv_diff
    double aa = dv[0]*dv[0] + dv[1]*dv[1]+dv[2]*dv[2]; // aa = dv@dv
    double bb = dv_diff[0]*dv_diff[0] + dv_diff[1]*dv_diff[1]+dv_diff[2]*dv_diff[2]; // bb = dv_diff@dv_diff
    double ad = std::sqrt(aa + 2*bb + 2*std::sqrt((ab*ab + bb*bb)))/tof;
    double mima = 2 * Tmax / ad / (1 + std::exp(-ad*tof/veff));
    return {mima, ad};
}

} // namespace kep3