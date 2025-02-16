// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_TA_PONTRYAGIN_CARTESIAN_H
#define kep3_TA_PONTRYAGIN_CARTESIAN_H

#include <vector>

#include <kep3/detail/visibility.hpp>

#include <heyoka/expression.hpp>
#include <heyoka/taylor.hpp>

namespace kep3::ta
{
// The Pontryagin Cartesian (pc) case.
// In the context of the TPBVP resulting from applying Pontryagin principle (minimum) to the NEP
// low-thrust problem in Cartesian Coordinates:

// Returns dynamics, hamiltonian, switching function control and control direction
// as expressions of the state-costate
kep3_DLL_PUBLIC std::tuple<std::vector<std::pair<heyoka::expression, heyoka::expression>>, heyoka::expression,
                           heyoka::expression, heyoka::expression, std::vector<heyoka::expression>>
pc_expression_factory();

// Returns the dynamics only. Offered for convenience and consistency within the ta namespace.
kep3_DLL_PUBLIC std::vector<std::pair<heyoka::expression, heyoka::expression>> pc_dyn();

// These return const references to function level static variables of type heyoka::taylor_adaptive<double>.
kep3_DLL_PUBLIC const heyoka::taylor_adaptive<double> &get_ta_pc(double tol);
kep3_DLL_PUBLIC const heyoka::taylor_adaptive<double> &
get_ta_pc_var(double tol); // variational (lx,ly,lz,lvx,lvy,lvz,lm,l0) first order

// Methods to access the cache dimensions.
kep3_DLL_PUBLIC size_t get_ta_pc_cache_dim();
kep3_DLL_PUBLIC size_t get_ta_pc_var_cache_dim();

// These return const references to function level static variables of type heyoka::cfunc<double>
// useful in the context of this Pontryagin cartesian OCP
kep3_DLL_PUBLIC const heyoka::cfunc<double> &get_pc_H_cfunc();
kep3_DLL_PUBLIC const heyoka::cfunc<double> &get_pc_SF_cfunc();
kep3_DLL_PUBLIC const heyoka::cfunc<double> &get_pc_u_cfunc();
kep3_DLL_PUBLIC const heyoka::cfunc<double> &get_pc_i_vers_cfunc();

} // namespace kep3::ta

#endif