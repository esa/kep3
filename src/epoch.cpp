// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <ratio>
#include <string>
#include <type_traits>

#include "kep3/core_astro/convert_julian_dates.hpp"
#include "kep3/epoch.hpp"

namespace kep3 {

/**
 * @brief Constructs a default epoch .
 */
epoch::epoch() = default;

/**
 * @brief Constructs an epoch from a non-Gregorian date.
 *
 * @param[in] epoch_in A double indicating the number of days
                        since day 0 in the specified calendar.
 * @param[in] epoch_type epoch::julian_type
 */
epoch::epoch(const double epoch_in, const julian_type epoch_type)
    : tp{make_tp(epoch_in, epoch_type)} {}

/**
 * @brief Constructs an epoch from offsets relative to 0 MJD2000.
 *
 * @param[in] y The number of years.
 * @param[in] d The number of days.
 * @param[in] h The number of hours.
 * @param[in] min The number of minutes.
 * @param[in] s The number of seconds.
 * @param[in] ms The number of milliseconds.
 * @param[in] us The number of microseconds.
 */
epoch::epoch(const int y, const int d, const int h, const int min, const int s,
             const int ms, const int us)
    : tp{make_tp(y, d, h, min, s, ms, us)} {}

/**
 * @brief Constructs an epoch from a const reference to a time point.
 *
 * @param[in] time_point Self-explanatory.
 */
epoch::epoch(const kep_clock::time_point &time_point) : tp{time_point} {}

/**
 * @brief Constructs an epoch from an rvalue reference to a time point.
 *
 * @param[in] time_point Self-explanatory.
 */
epoch::epoch(kep_clock::time_point &&time_point) : tp{time_point} {}

kep_clock::time_point epoch::make_tp(const int y, const int d, const int h,
                                     const int min, const int s, const int ms,
                                     const int us)

{
  return kep_clock::time_point{} + chr::years(y) + chr::days(d) +
         chr::hours(h) + chr::minutes(min) + chr::seconds(s) +
         chr::milliseconds(ms) + chr::microseconds(us);
}

kep_clock::time_point epoch::make_tp(const double epoch_in,
                                     const julian_type epoch_type) {
  switch (epoch_type) {
  case julian_type::MJD2000:
    return epoch::tp_from_days(epoch_in);
  case julian_type::MJD:
    return mjd2mjd2000(epoch::tp_from_days(epoch_in));
  case julian_type::JD:
    return jd2mjd2000(epoch::tp_from_days(epoch_in));
  default:
    throw;
  }
}

/**
 * @brief Creates time point from the number of days since 0 MJD2000.
 *
 * @return A time point
 */
constexpr kep_clock::time_point epoch::tp_from_days(const double days) {
  return kep_clock::time_point{} +
         chr::duration_cast<kep_clock::duration>(
             chr::duration<double, std::ratio<86400>>(days));
}

/**
 * @brief Returns a time point formatted as a date/time string
 *
 * @param tp The time point.
 * @return A formatted date/time string.
 */
auto epoch::as_utc_string(const kep_clock::time_point &tp) {
  auto t = kep_clock::to_time_t(tp);
  return std::put_time(gmtime(&t), "%FT%T");
}

/**
 * @brief Streams out a date in the format 2000-Jan-01 00:12:30
 *
 * @param[in] s Stream to which the epoch will be sent.
 * @param[in] ep Epoch to be sent to the stream.
 *
 * @return Reference to s.
 */
std::ostream &operator<<(std::ostream &s, const epoch &ep) {
  s << epoch::as_utc_string(ep.tp);
  return s;
}

bool operator>(const epoch &c1, const epoch &c2) { return c1.tp > c2.tp; }
bool operator<(const epoch &c1, const epoch &c2) { return c1.tp < c2.tp; }
bool operator>=(const epoch &c1, const epoch &c2) { return c1.tp >= c2.tp; }
bool operator<=(const epoch &c1, const epoch &c2) { return c1.tp <= c2.tp; }
bool operator==(const epoch &c1, const epoch &c2) { return c1.tp == c2.tp; }
bool operator!=(const epoch &c1, const epoch &c2) { return c1.tp != c2.tp; }

kep_clock::duration operator-(const epoch &lhs, const epoch &rhs) {
  return lhs.tp - rhs.tp;
}

} // namespace kep3
