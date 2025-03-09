// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_MIMA_H
#define kep3_MIMA_H

#include <array>
#include <utility>

#include <kep3/detail/visibility.hpp>

namespace kep3
{

kep3_DLL_PUBLIC std::pair<double, double> mima(std::array<double, 3> dv1, std::array<double, 3> dv2, double tof, double Tmax, double veff);

} // namespace kep3
#endif // kep3_IC2EQ2IC_H