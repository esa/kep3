// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_EPOCH_HPP
#define kep3_EPOCH_HPP

#include <chrono>
#include <cstdint>
#include <ctime>
#include <fmt/ostream.h>
#include <iostream>

#include <kep3/detail/s11n.hpp>
#include <kep3/detail/visibility.hpp>
#include <ratio>
#include <type_traits>
#include <utility>

/// Keplerian Toolbox
/**
 * This namespace contains astrodynamics and space flight mechanics routines
 * that are related to keplerian motions or models.
 */

namespace kep3
{
using namespace std::literals;
namespace chr = std::chrono;
using uint = unsigned int;

template <typename T>
struct is_duration : std::false_type {
};

template <typename Rep, typename Period>
struct is_duration<chr::duration<Rep, Period>> : std::true_type {
};

template <typename T>
using enable_if_is_duration = std::enable_if_t<is_duration<T>::value>;

struct kep_clock : public chr::system_clock {

    /**
     * @brief Custom clock.
     * Used for constructing epochs with a custom reference point (1 Jan 2000).
     * By defining a custom clock, we avoid the overflow
     * that std::chrono::system_clock suffers at +/- 292 years.
     * To do that, we lower the resolution to 1 us (microsecond),
     * which gives the clock a range of +/- 292 thousand years.
     *
     * NOTE: Adding durations of less than 1 us to an epoch (defined below)
     * would not be registered.
     *
     * NOTE: As of C++20, the standard guarantees that std::chrono::system_clock
     * uses the UNIX time reference point, which is midnight on 1 January 1970
     * (1970-01-01T00:00:00). We correct for that here in order to bring the
     * reference point forward to midnight on 1 January 2000
     * (2000-01-01T00:00:00), which is 0 MJD2000.
     */
    using rep = int_fast64_t;
    // Resolution of (1 / 1'000'000)s = 1 us
    using period = std::ratio<1, 1'000'000>;
    using duration = chr::duration<rep, period>;
    using time_point = chr::time_point<kep_clock, duration>;
    static constexpr bool is_steady = false;
    // Number of seconds from midnight on 1 Jan 1970 to midnight on 1 Jan 2000.
    static constexpr chr::seconds y2k_offset{946684800s};

    static constexpr time_point ref_epoch{kep_clock::time_point{} + y2k_offset};

    static constexpr std::time_t to_time_t(const time_point &t) noexcept
    {
        return static_cast<std::time_t>(chr::duration_cast<chr::seconds>(t.time_since_epoch() + y2k_offset).count());
    }

    static constexpr time_point from_time_t(std::time_t t) noexcept
    {
        return chr::time_point_cast<duration>(time_point(chr::seconds(t) - y2k_offset));
    }
};

/// epoch class.
/**
 * This class defines and contains a non-Gregorian date (i.e., a date expressed
 * in Julian format). The date is defined in MJD2000 format as a
 * kep_clock::time_point private member. Types of non-Gregorian dates supported:
 *      - Julian Date (JD): the number of days passed since January 1, 4713 BC
 * at 12:00 (noon).
 *      - Modified Julian Date (MJD): the number of days passed since November
 * 17, 1858 at 00:00 (midnight).
 *      - Modified Julian Date 2000 (MJD2000): the number of days passed since
 * Juanuary 1, 2000 at 00:00 (midnight).
 */
class kep3_DLL_PUBLIC epoch
{
public:
    enum class julian_type { MJD2000, MJD, JD };

    /** Constructors */
    // Default constructor
    epoch();

    // Constructor from a julian date (as a floating-point value)
    explicit epoch(double epoch_in, julian_type epoch_type = julian_type::MJD2000);

    /**
     * Constructs an epoch from a std::chrono::duration.
     * The reference point is assumed to be MJD2000.
     * \param[in] time The time as a duration.
     */
    template <class Duration, class = enable_if_is_duration<Duration>>
    explicit epoch(const Duration &duration) : tp{kep_clock::time_point{} + duration}
    {
    }

    // Constructor from duration&&)
    template <class Duration, class = enable_if_is_duration<Duration>>
    explicit epoch(Duration &&duration) : tp{kep_clock::time_point{} + std::forward(duration)}
    {
    }

    // Constructor for datetime broken down into its constituents.
    explicit epoch(int y, uint mon, uint d, int h = 0, int min = 0, int s = 0, int ms = 0, int us = 0);

    /* Computing non-Gregorian dates */

    /**
     * @return Number of days since 0 JD (including fractional days).
     */
    [[nodiscard]] constexpr double jd() const
    {
        return chr::duration<double, std::ratio<86400>>(tp.time_since_epoch() - kep_clock::y2k_offset + 211813444800s)
            .count();
    }

    /**
     * @return Number of days since 0 MJD (including fractional days).
     */
    [[nodiscard]] constexpr double mjd() const
    {
        return chr::duration<double, std::ratio<86400>>(tp.time_since_epoch() - kep_clock::y2k_offset + 4453401600s)
            .count();
    }

    /**
     * @return Number of days since 0 MJD2000 (including fractional days).
     */
    [[nodiscard]] constexpr double mjd2000() const
    {
        return chr::duration<double, std::ratio<86400>>(tp.time_since_epoch() - kep_clock::y2k_offset).count();
    }

    /* Helper functions for constructors */
    static kep_clock::time_point make_tp(int y, uint mon, uint d, int h = 0, int min = 0,
                                         int s = 0, int ms = 0, int us = 0);

    static kep_clock::time_point make_tp(double epoch_in, julian_type epoch_type);

    // Conversions
    static constexpr kep_clock::time_point tp_from_days(double days);

    // Duration conversions
    static constexpr double as_sec(kep_clock::duration &&d)
    {
        return std::chrono::duration<double, std::chrono::seconds::period>(d).count();
    }

    // Printing
    [[nodiscard]] std::string as_utc_string() const;
    static std::string as_utc_string(const kep_clock::time_point &tp);

    /** operator overloads for sum and diff (epoch-days) and comparison
     * operators
     * **/

    kep3_DLL_PUBLIC friend std::ostream &operator<<(std::ostream &s, epoch const &epoch_in);

    template <class Duration, class = enable_if_is_duration<Duration>>
    epoch &operator+=(const Duration &duration)
    {
        tp += chr::duration_cast<kep_clock::duration>(duration);
        return *this;
    }

    template <class Duration, class = enable_if_is_duration<Duration>>
    epoch &operator-=(const Duration &duration)
    {
        tp -= chr::duration_cast<kep_clock::duration>(duration);
        return *this;
    }

    kep3_DLL_PUBLIC friend bool operator>(const epoch &c1, const epoch &c2);
    kep3_DLL_PUBLIC friend bool operator<(const epoch &c1, const epoch &c2);
    kep3_DLL_PUBLIC friend bool operator>=(const epoch &c1, const epoch &c2);
    kep3_DLL_PUBLIC friend bool operator<=(const epoch &c1, const epoch &c2);
    kep3_DLL_PUBLIC friend bool operator==(const epoch &c1, const epoch &c2);
    kep3_DLL_PUBLIC friend bool operator!=(const epoch &c1, const epoch &c2);

    template <class Duration, class = enable_if_is_duration<Duration>>
    epoch operator+(const Duration &duration)
    {
        return epoch(tp + chr::duration_cast<kep_clock::duration>(duration));
    }

    template <class Duration, class = enable_if_is_duration<Duration>>
    epoch operator-(const Duration &duration)
    {
        return epoch(tp - chr::duration_cast<kep_clock::duration>(duration));
    }

    static constexpr auto days(const double value)
    {
        return chr::duration_cast<kep_clock::duration>(chr::duration<double, std::ratio<86400>>(value));
    }

    static constexpr auto sec(const double value)
    {
        return chr::duration_cast<kep_clock::duration>(chr::duration<double, std::ratio<1>>(value));
    }

    kep3_DLL_PUBLIC friend kep_clock::duration operator-(const epoch &lhs, const epoch &rhs);

private:
    // Constructor for const time_point&)
    explicit epoch(const kep_clock::time_point &time_point);

    // Constructor for const time_point&&)
    explicit epoch(kep_clock::time_point &&time_point);

    // Serialization code
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const uint)
    {
        ar &boost::serialization::make_binary_object(&tp, sizeof(tp));
    }
    // Serialization code (END)

    // Time point relative to 1 Jan 2000 (MJD2000)
    kep_clock::time_point tp;
};

kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &s, const epoch &epoch_in);

} // end of namespace kep3

template <>
struct fmt::formatter<kep3::epoch> : fmt::ostream_formatter {
};

#endif // kep3_EPOCH_HPP
