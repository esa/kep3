// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_IC2EQ2IC_H
#define kep3_IC2EQ2IC_H

#include <array>

#include <kep3/detail/visibility.hpp>

namespace kep3
{

kep3_DLL_PUBLIC std::array<double, 6> ic2eq(const std::array<std::array<double, 3>, 2> &pos_vel, double mu,
                                            bool retrogade = false);

kep3_DLL_PUBLIC std::array<std::array<double, 3>, 2> eq2ic(const std::array<double, 6> &eq, double mu,
                                                           bool retrogade = false);

} // namespace kep3
#endif // kep3_IC2EQ2IC_H