// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <map>
#include <mutex>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <heyoka/config.hpp>
#include <heyoka/expression.hpp>
#include <heyoka/math/log.hpp>
#include <heyoka/math/pow.hpp>
#include <heyoka/math/sqrt.hpp>
#include <heyoka/math/sum.hpp>
#include <heyoka/taylor.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/ta/pontryagin_cartesian.hpp>

using heyoka::diff;
using heyoka::expression;
using heyoka::log;
using heyoka::make_vars;
using heyoka::par;
using heyoka::pow;
using heyoka::sqrt;
using heyoka::taylor_adaptive;
using heyoka::var_ode_sys;

namespace kep3::ta
{
// All the relevant expression in the TPBVP of the Cartesian Low-Thrust indirect OCP are created here.
std::tuple<std::vector<std::pair<expression, expression>>, expression, expression, expression, std::vector<expression>>
pc_expression_factory()
{
    // The state
    auto [x, y, z, vx, vy, vz, m] = make_vars("x", "y", "z", "vx", "vy", "vz", "m");

    // The costate
    auto [lx, ly, lz, lvx, lvy, lvz, lm] = make_vars("lx", "ly", "lz", "lvx", "lvy", "lvz", "lm");
    // The controls
    auto [u, ix, iy, iz] = make_vars("u", "ix", "iy", "iz");

    // Useful expressions
    auto r3 = pow(x * x + y * y + z * z, 1.5);
    auto lv_norm = sqrt(lvx * lvx + lvy * lvy + lvz * lvz);
    auto rho = 1. - par[1] * lv_norm / m / par[4] - par[2] * lm / par[4];

    // Hamiltonian
    auto H_full = lx * vx + ly * vy + lz * vz + lvx * (par[1] * u / m * ix - (par[0] / r3) * x)
                  + lvy * (par[1] * u / m * iy - (par[0] / r3) * y) + lvz * (par[1] * u / m * iz - (par[0] / r3) * z)
                  - lm * par[2] * u + par[4] * (u - par[3] * log(u * (1. - u)));

    // #Augmented equations of motion
    auto rhs = std::vector<expression>{};
    rhs.reserve(14);
    for (const auto &var : {lx, ly, lz, lvx, lvy, lvz, lm, x, y, z, vx, vy, vz, m}) {
        rhs.push_back(diff(H_full, var));
    }
    // The costate equations have the minus sign (Hamiltonian formalism)
    for (auto j = 7u; j < 14; ++j) {
        rhs[j] = -rhs[j];
    }

    // We apply Pontryagin minimum principle(primer vector and u ^ * = 2eps / (rho + 2eps + sqrt(rho ^ 2 + 4 * eps ^
    // 2)))
    std::map<expression, expression> argmin_H_full = {
        {ix, -lvx / lv_norm},
        {iy, -lvy / lv_norm},
        {iz, -lvz / lv_norm},
        {u, 2. * par[3] / (rho + 2. * par[3] + sqrt(rho * rho + 4. * par[3] * par[3]))},
    };
    rhs = subs(rhs, argmin_H_full);

    // We build the Hamiltonian as a function of the state / co - state only
    //(i.e.no longer of controls)
    auto H = subs(H_full, argmin_H_full);

    // We assemble the system dynamics
    std::vector<std::pair<expression, expression>> sys{};
    sys.reserve(14);
    std::vector<expression> vars{x, y, z, vx, vy, vz, m, lx, ly, lz, lvx, lvy, lvz, lm};
    for (auto i = 0u; i < 14; ++i) {
        sys.emplace_back(vars[i], rhs[i]);
    }
    return {sys, H, rho, argmin_H_full[u], {argmin_H_full[ix], argmin_H_full[iy], argmin_H_full[iz]}};
}

std::vector<std::pair<expression, expression>> pc_dyn()
{
    return std::get<0>(pc_expression_factory());
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_pc_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_pc_cache;

const taylor_adaptive<double> &get_ta_pc(double tol)
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_pc_mutex);

    // Lookup.
    if (auto it = ta_pc_cache.find(tol); it == ta_pc_cache.end()) {
        // Cache miss, create new one.
        auto new_ta = taylor_adaptive<double>{std::get<0>(pc_expression_factory()), heyoka::kw::tol = tol};
        return ta_pc_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        // Cache hit, return existing.
        return it->second;
    }
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_pc_var_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_pc_var_cache;

const taylor_adaptive<double> &get_ta_pc_var(double tol)
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_pc_var_mutex);

    // Lookup.
    if (auto it = ta_pc_var_cache.find(tol); it == ta_pc_var_cache.end()) {
        auto [lx, ly, lz, lvx, lvy, lvz, lm] = make_vars("lx", "ly", "lz", "lvx", "lvy", "lvz", "lm");
        auto vsys = var_ode_sys(std::get<0>(pc_expression_factory()), {lx, ly, lz, lvx, lvy, lvz, lm, par[4]}, 1);
        // Cache miss, create new one.
        auto new_ta = taylor_adaptive<double>{vsys, heyoka::kw::tol = tol, heyoka::kw::compact_mode = true};
        return ta_pc_var_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        // Cache hit, return existing.
        return it->second;
    }
}

size_t get_ta_pc_cache_dim()
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_pc_mutex);
    return ta_pc_cache.size();
}

size_t get_ta_pc_var_cache_dim()
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_pc_mutex);
    return ta_pc_var_cache.size();
}

// Factory functions to help the static variable initialization later
auto pc_H_cfunc_factory()
{
    auto [x, y, z, vx, vy, vz, m, lx, ly, lz, lvx, lvy, lvz, lm]
        = make_vars("x", "y", "z", "vx", "vy", "vz", "m", "lx", "ly", "lz", "lvx", "lvy", "lvz", "lm");
    return heyoka::cfunc<double>({std::get<1>(pc_expression_factory())},
                                 {x, y, z, vx, vy, vz, m, lx, ly, lz, lvx, lvy, lvz, lm});
}
auto pc_SF_cfunc_factory()
{
    auto [x, y, z, vx, vy, vz, m, lx, ly, lz, lvx, lvy, lvz, lm]
        = make_vars("x", "y", "z", "vx", "vy", "vz", "m", "lx", "ly", "lz", "lvx", "lvy", "lvz", "lm");
    return heyoka::cfunc<double>({std::get<2>(pc_expression_factory())},
                                 {x, y, z, vx, vy, vz, m, lx, ly, lz, lvx, lvy, lvz, lm});
}
auto pc_u_cfunc_factory()
{
    auto [x, y, z, vx, vy, vz, m, lx, ly, lz, lvx, lvy, lvz, lm]
        = make_vars("x", "y", "z", "vx", "vy", "vz", "m", "lx", "ly", "lz", "lvx", "lvy", "lvz", "lm");
    return heyoka::cfunc<double>({std::get<3>(pc_expression_factory())},
                                 {x, y, z, vx, vy, vz, m, lx, ly, lz, lvx, lvy, lvz, lm});
}
auto pc_i_vers_cfunc_factory()
{
    auto [lvx, lvy, lvz] = make_vars("lvx", "lvy", "lvz");
    return heyoka::cfunc<double>({std::get<4>(pc_expression_factory())}, {lvx, lvy, lvz});
}

// Function-level static variable: it is initialised the first time the function is invoked
// and the initialisation is guaranteed to be thread-safe (the compiler will put the necessary
// mutex and locking).
const heyoka::cfunc<double> &get_pc_H_cfunc()
{
    static const auto pc_H_cfunc = pc_H_cfunc_factory();
    return pc_H_cfunc;
}

const heyoka::cfunc<double> &get_pc_SF_cfunc()
{
    static const auto pc_SF_cfunc = pc_SF_cfunc_factory();
    return pc_SF_cfunc;
}
const heyoka::cfunc<double> &get_pc_u_cfunc()
{
    static const auto pc_u_cfunc = pc_u_cfunc_factory();
    return pc_u_cfunc;
}
const heyoka::cfunc<double> &get_pc_i_vers_cfunc()
{
    static const auto pc_i_vers_cfunc = pc_i_vers_cfunc_factory();
    return pc_i_vers_cfunc;
}

} // namespace kep3::ta