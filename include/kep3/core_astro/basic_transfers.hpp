// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_IMPULSIVE_TRANSFERS_H
#define kep3_IMPULSIVE_TRANSFERS_H

#include <utility> // For std::pair

#include <kep3/detail/visibility.hpp>

namespace kep3
{
kep3_DLL_PUBLIC std::pair<double, std::pair<double, double>> hohmann(double r1, double r2, double mu);
} // namespace kep3
#endif // kep3_IMPULSIVE_TRANSFERS_H