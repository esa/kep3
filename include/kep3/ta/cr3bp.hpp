// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_TA_CR3BP_H
#define kep3_TA_CR3BP_H

#include <tuple>
#include <vector>

#include <kep3/detail/visibility.hpp>

#include <heyoka/expression.hpp>
#include <heyoka/taylor.hpp>

namespace kep3::ta
{
// Returns the low-thrust dynamics (heyoka API) of the cr3bp. 6 states, 1 parameter: mu.
kep3_DLL_PUBLIC std::vector<std::pair<heyoka::expression, heyoka::expression>> cr3bp_dyn();

// These return const references to function level static variables of type heyoka::taylor_adaptive<double>.
// NOTE: The object retruned are expected to be copied to then be modified.
kep3_DLL_PUBLIC const heyoka::taylor_adaptive<double> &get_ta_cr3bp(double tol);
kep3_DLL_PUBLIC const heyoka::taylor_adaptive<double> &get_ta_cr3bp_var(double tol); // variational (x,y,z,vx,vy,vz,ux,uy,uz) first order

// Methods to access the cache dimensions.
kep3_DLL_PUBLIC size_t get_ta_cr3bp_cache_dim();
kep3_DLL_PUBLIC size_t get_ta_cr3bp_var_cache_dim();
} // namespace kep3::ta

#endif