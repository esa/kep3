// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_CONSTANTS_H
#define kep3_CONSTANTS_H

#include <boost/math/constants/constants.hpp>

namespace kep3
{

enum elements_type {
    KEP_M,  // Keplerian Osculating (with Mean Anomaly)
    KEP_F,  // Keplerian Osculating (with True Anomaly)
    MEQ,    // Modified Equinoctial Elements
    MEQ_R,  // Modified Equinoctial Elements (retrogade)
    POSVEL, // position and Velocity
};

inline constexpr double pi = boost::math::constants::pi<double>();
inline constexpr double half_pi = boost::math::constants::half_pi<double>();

inline constexpr double AU = 149597870700.0;                 // Astronomical Unit (m)
inline constexpr double CAVENDISH = 73.6687e-11;             // Cavendish constant (N M^2 / kg^2)
inline constexpr double MU_SUN = 1.32712440018e20;           // Sun's gravitational parameter (m^3/s^2 kg)
inline constexpr double MU_EARTH = 398600441800000.0;        // Earth's gravitational parameter (m^3/s^2 kg)
inline constexpr double EARTH_VELOCITY = 29784.691831696804; // Earth's velocity. (m/s)
inline constexpr double EARTH_J2 = 1.08262668E-03;
inline constexpr double EARTH_RADIUS = 6378137; // Earth's radius (m)
inline constexpr double DEG2RAD = (pi / 180.0);
inline constexpr double RAD2DEG = (180.0 / pi);
inline constexpr double DAY2SEC = 86400.0;
inline constexpr double SEC2DAY = (1. / DAY2SEC);
inline constexpr double DAY2YEAR = (1. / 365.25);
inline constexpr double G0 = 9.80665; // Acceleration at Earth's surface (m/s^2)

} // namespace kep3

#endif // kep3_CONSTANTS_H