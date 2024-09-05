// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <fmt/core.h>

#include <heyoka/taylor.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/propagate_taylor.hpp>
#include <optional>

namespace kep3
{

/// Taylor propagation for the Stark problem
/**
 * This function propagates an initial Cartesian state for a time t assuming a central body and a keplerian motion
 * perturbed by a constant inertial thrust. Heyoka Taylor adaptive integration is used as a basic numerical technique.
 * All units systems can be used, as long as the input parameters are all expressed in the same system. The information
 * on the state transition matrix can be optionally asked, in which case also the sensitivities to the throttles are
 * returned.
 */

using heyoka::taylor_outcome;
using ta::get_ta_lt_kepler;
using ta::get_ta_lt_kepler_var;

std::tuple<std::array<std::array<double, 3>, 2>, double,
           std::optional<std::pair<std::array<double, 49>, std::array<double, 21>>>>
propagate_taylor(const std::array<std::array<double, 3>, 2> &pos_vel, double m, std::array<double, 3> thrust,
                 // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                 double tof, double mu, double veff, double tol, bool stm)
{
    // Assemble the return value
    if (!stm) {
        std::array<std::array<double, 3>, 2> retval_posvel{};

        auto ta_cached = get_ta_lt_kepler(tol);
        // We copy the const cached one
        heyoka::taylor_adaptive<double> ta(ta_cached);
        // We set the parameters mu, veff and the thrust magnitude.
        *(ta.get_pars_data()) = mu;
        *(ta.get_pars_data() + 1) = veff;
        std::copy(thrust.begin(), thrust.end(), ta.get_state_data() + 2);
        // Set the Taylor Integration initial conditions
        ta.set_time(0.);
        std::copy(pos_vel[0].begin(), pos_vel[0].end(), ta.get_state_data());
        std::copy(pos_vel[1].begin(), pos_vel[1].end(), ta.get_state_data() + 3);
        *(ta.get_state_data() + 6) = m;
        // ... and integrate
        auto out = ta.propagate_until(tof);
        if (std::get<0>(out) != taylor_outcome::time_limit) {
            throw std::domain_error("propagate_taylor: failiure to reach the final time requested");
        }
        // We now copy the result into the various return values
        std::copy(ta.get_state().begin(), ta.get_state().begin() + 3, retval_posvel[0].begin());
        std::copy(ta.get_state().begin() + 3, ta.get_state().begin() + 6, retval_posvel[1].begin());
        return {retval_posvel, ta.get_state()[6], std::nullopt};
    }
    // The stm has been requested, using the variational integrator
    std::array<std::array<double, 3>, 2> retval_posvel{};
    std::array<double, 49> dxdx{};
    std::array<double, 21> dxdu{};

    auto ta_cached = get_ta_lt_kepler_var(tol);
    // We copy the const cached one
    heyoka::taylor_adaptive<double> ta(ta_cached);
    // We set the parameters mu, veff and the thrust magnitude.
    *(ta.get_pars_data()) = mu;
    *(ta.get_pars_data() + 1) = veff;
    std::copy(thrust.begin(), thrust.end(), ta.get_state_data() + 2);
    // Set the Taylor Integration initial conditions
    ta.set_time(0.);
    std::copy(pos_vel[0].begin(), pos_vel[0].end(), ta.get_state_data());
    std::copy(pos_vel[1].begin(), pos_vel[1].end(), ta.get_state_data() + 3);
    *(ta.get_state_data() + 6) = m;
    // ... and integrate
    auto out = ta.propagate_until(tof);
    if (std::get<0>(out) != taylor_outcome::time_limit) {
        throw std::domain_error("propagate_taylor (with stm): failiure to reach the final time requested");
    }
    // We now copy the result into the various return values
    std::copy(ta.get_state().begin(), ta.get_state().begin() + 3, retval_posvel[0].begin());
    std::copy(ta.get_state().begin() + 3, ta.get_state().begin() + 6, retval_posvel[1].begin());
    for (auto i = 0; i < 7; ++i) {
        std::copy(ta.get_state().begin() + 7 + 10l * i, ta.get_state().begin() + 7 + 10l * i + 7,
                  dxdx.begin() + 7l * i);
        std::copy(ta.get_state().begin() + 14 + 10l * i, ta.get_state().begin() + 14 + 10l * i + 3,
                  dxdu.begin() + 3l * i);
    }
    return {retval_posvel, ta.get_state()[6], std::make_pair(dxdx, dxdu)};
}

} // namespace kep3