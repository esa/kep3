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
#include <fmt/core.h>
#include <fmt/chrono.h>
#include <iomanip>
#include <iostream>
#include <ratio>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>


#include "kep3/core_astro/convert_julian_dates.hpp"
#include "kep3/epoch.hpp"

namespace kep3
{

/**
 * @brief Constructs a default epoch .
 */
epoch::epoch()
    : tp{kep_clock::ref_epoch} {}

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
epoch::epoch(const uint y, const uint mon, const uint d, const uint h,
             const uint min, const uint s, const uint ms, const uint us)
    : tp{make_tp(y, mon, d, h, min, s, ms, us)} {}

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
epoch::epoch(kep_clock::time_point &&time_point) : tp{std::move(time_point)} {}

kep_clock::time_point epoch::make_tp(const uint y, const uint mon, const uint d, const uint h,
                                     const uint min, const uint s, const uint ms, const uint us)

{
    return kep_clock::time_point{} + chr::sys_days(chr::year_month_day{chr::year(y) / chr::month(mon) / chr::day(d)} + chr::months{0}).time_since_epoch() +
         chr::hours(h) + chr::minutes(min) + chr::seconds(s) +
         chr::milliseconds(ms) + chr::microseconds(us);
}

kep_clock::time_point epoch::make_tp(const double epoch_in,
                                     const julian_type epoch_type) {
  switch (epoch_type) {
  case julian_type::MJD2000:
    return epoch::tp_from_days(epoch_in);
  case julian_type::MJD:
    return epoch::tp_from_days(epoch_in) - 4453401600s;
  case julian_type::JD:
    return epoch::tp_from_days(epoch_in) - 211813444800s;
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
  return kep_clock::ref_epoch +
         chr::duration_cast<kep_clock::duration>(
             chr::duration<double, std::ratio<86400>>(days));
}

/**
 * @brief Returns a time point formatted as a date/time string
 * in the in the format 2000-12-31T12:34:56.123456.
 *
 * @param tp The time point.
 * @return A formatted date/time string.
 */
std::string epoch::as_utc_string(const kep_clock::time_point& tp) {
    std::stringstream iss;
    const auto tse{tp.time_since_epoch()};
    const auto dp{std::chrono::duration_cast<std::chrono::days>(tse)};
    const auto hms{std::chrono::duration_cast<std::chrono::seconds>(tse - dp)};
    const auto us{std::chrono::duration_cast<std::chrono::microseconds>(tse - dp - hms)};
    iss << fmt::format("{:%F}", chr::sys_days(dp)) << "T" << fmt::format("{:%T}", hms) << "." << fmt::format("{:06}", us.count());
    return iss.str();
}

std::string epoch::as_utc_string() const {
    return epoch::as_utc_string(tp);
}

/**
 * @brief Streams out an epoch as a UTC string.
 *
 * @param[in] s Stream to which the epoch will be sent.
 * @param[in] ep Epoch to be sent to the stream.
 *
 * @return Reference to s.
 */
std::ostream &operator<<(std::ostream &s, const epoch &ep) {
  s << ep.as_utc_string();
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
