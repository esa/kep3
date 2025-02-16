// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_EQ2PAR2EQ_H
#define kep3_EQ2PAR2EQ_H

#include <array>

#include <kep3/detail/visibility.hpp>

namespace kep3
{

kep3_DLL_PUBLIC std::array<double, 6> eq2par(const std::array<double, 6> &eq, bool retrogade = false);

kep3_DLL_PUBLIC std::array<double, 6> par2eq(const std::array<double, 6> &par, bool retrogade = false);

} // namespace kep3
#endif // kep3_EQ2PAR2EQ_H