// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <chrono>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <boost/numeric/conversion/cast.hpp>

#include <fmt/chrono.h>
#include <fmt/core.h>

#include "kep3/epoch.hpp"

namespace kep3
{

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
epoch::epoch(double epoch_in, julian_type epoch_type)
{
    switch (epoch_type) {
        case julian_type::MJD2000:
            m_tp = epoch::tp_from_days(epoch_in);
            break;
        case julian_type::MJD:
            m_tp = epoch::tp_from_days(epoch_in) - mjd_offset;
            break;
        case julian_type::JD:
            m_tp = epoch::tp_from_days(epoch_in) - jd_offset;
            break;
        default:
            throw std::invalid_argument(
                fmt::format("An unsupported julian_type enumerator with value {} was used in the epoch constructor",
                            static_cast<std::underlying_type_t<julian_type>>(epoch_type)));
    }
}

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
epoch::epoch(const std::int32_t y, const std::uint32_t mon, const std::uint32_t d, const std::int32_t h, // NOLINT
             const std::int32_t min, const std::int32_t s, const std::int32_t ms, const std::int32_t us)
    : m_tp{make_tp(y, mon, d, h, min, s, ms, us)}
{
}

// Epoch constructor from string
epoch::epoch(const std::string &in, string_format sf)
{
    // NOTE: only ISO format supported so far.
    if (sf != string_format::ISO) {
        throw std::invalid_argument(
            fmt::format("An unsupported string_format enumerator with value {} was used in the epoch constructor",
                        static_cast<std::underlying_type_t<string_format>>(sf)));
    }

    // We assume: 1980-10-17T11:36:21.121841 and allow crops such as 1980-10.
    constexpr std::array<decltype(in.size()), 11> allowed_lenghts{7, 10, 13, 16, 19, 21, 22, 23, 24, 25, 26};
    auto len = in.size();
    auto foo = std::find(std::begin(allowed_lenghts), std::end(allowed_lenghts), len); // NOLINT
    if (foo == std::end(allowed_lenghts)) {
        throw std::logic_error(
            "Malformed input string when constructing an epoch. Must be 'YYYY-MM-DD HH:MM:SS:XXXXXX'. "
            "D,H,M,S and X can be missing incrementally.");
    }
    unsigned d = 1u;
    int h = 0, min = 0, s = 0, us = 0;
    int y = std::stoi(in.substr(0, 4));
    auto mon = boost::numeric_cast<unsigned>(std::stoi(in.substr(5, 2)));
    if (len >= 10) {
        d = boost::numeric_cast<unsigned>(std::stoi(in.substr(8, 2)));
        if (len >= 13) {
            h = std::stoi(in.substr(11, 2));
            if (len >= 16) {
                min = std::stoi(in.substr(14, 2));
                if (len >= 19) {
                    s = std::stoi(in.substr(17, 2));
                    if (len >= 21) {
                        std::string rest = in.substr(20);
                        us = std::stoi(rest);
                        for (decltype(rest.size()) i = 0; i < 6u - rest.size(); ++i) {
                            us *= 10;
                        }
                    }
                }
            }
        }
    }
    m_tp = make_tp(y, mon, d, h, min, s, 0, us);
}

/**
 * @brief Constructs an epoch from a time point.
 *
 * @param[in] time_point Self-explanatory.
 */
epoch::epoch(time_point tp) : m_tp{tp} {}

time_point epoch::make_tp(const std::int32_t y, const std::uint32_t mon, const std::uint32_t d, const std::int32_t h,
                          const std::int32_t min, const std::int32_t s, const std::int32_t ms, const std::int32_t us)

{
    return /*time_point{}
           +*/
        static_cast<time_point>(
            std::chrono::sys_days(
                std::chrono::year_month_day{std::chrono::year(y) / std::chrono::month(mon) / std::chrono::day(d)}
                + std::chrono::months{0})
                .time_since_epoch()
            + std::chrono::hours(h) + std::chrono::minutes(min) + std::chrono::seconds(s)
            + std::chrono::milliseconds(ms) + std::chrono::microseconds(us));
}

/**
 * @brief Creates time point from the number of days since 0 MJD2000.
 *
 * @return A time point
 */
time_point epoch::tp_from_days(double days)
{
    return y2k + std::chrono::duration_cast<microseconds>(std::chrono::duration<double, seconds_day_ratio>(days));
}

/**
 * @brief Returns a time point formatted as a date/time string
 * in the in the format 2000-12-31T12:34:56.123456.
 *
 * @param tp The time point.
 * @return A formatted date/time string.
 */
std::string epoch::as_utc_string(const time_point &tp)
{
    std::stringstream iss;
    const auto tse{tp.time_since_epoch()};
    const auto dp{std::chrono::floor<std::chrono::days>(tse)};
    const auto hms{std::chrono::floor<std::chrono::seconds>(tse - dp)};
    const auto us{std::chrono::floor<std::chrono::microseconds>(tse - dp - hms)};
    iss << fmt::format("{:%F}", std::chrono::sys_days(dp)) << "T" << fmt::format("{:%T}", hms) << "."
        << fmt::format("{:06}", us.count());
    return iss.str();
}

std::string epoch::as_utc_string() const
{
    return epoch::as_utc_string(m_tp);
}

/**
 * @brief Streams out an epoch as a UTC string.
 *
 * @param[in] s Stream to which the epoch will be sent.
 * @param[in] ep Epoch to be sent to the stream.
 *
 * @return Reference to s.
 */
std::ostream &operator<<(std::ostream &s, const epoch &ep)
{
    s << ep.as_utc_string();
    return s;
}

microseconds operator-(const epoch &lhs, const epoch &rhs)
{
    return lhs.m_tp - rhs.m_tp;
}

time_point epoch::get_tp() const
{
    return m_tp;
}
epoch utc_now()
{
    auto tp = time_point{std::chrono::duration_cast<microseconds>(std::chrono::system_clock::now().time_since_epoch())};

    return epoch(tp);
}

} // namespace kep3
