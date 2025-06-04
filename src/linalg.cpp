// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <kep3/linalg.hpp>

namespace kep3::linalg
{
mat33 _skew(const mat31 &v)
{
    return {{0., -v(2, 0), v(1, 0)}, {v(2, 0), 0., -v(0, 0)}, {-v(1, 0), v(0, 0), 0.}};
}

mat31 _cross(const mat31 &v1, const mat31 &v2)
{
    return {{v1(1, 0) * v2(2, 0) - v2(1, 0) * v1(2, 0), v2(0, 0) * v1(2, 0) - v1(0, 0) * v2(2, 0),
             v1(0, 0) * v2(1, 0) - v2(0, 0) * v1(1, 0)}};
}

void _normalize(std::array<double, 3> &v1)
{
    const double norm = std::sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
    v1[0] /= norm;
    v1[1] /= norm;
    v1[2] /= norm;
}
} // namespace kep3::linalg
