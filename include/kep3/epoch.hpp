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

#include <kep3/detail/duration.hpp>
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
namespace chr = std::chrono;

template <typename T>
concept Duration = detail::is_duration<T>::value;



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
     * (1970-01-01T00:00:00).
     */
    using rep = int_fast64_t;
    // Resolution of (1 / 1'000'000)s = 1 us
    using period = std::ratio<1, 1'000'000>;
    using duration = chr::duration<rep, period>;
    using time_point = chr::time_point<kep_clock, duration>;
    static constexpr bool is_steady = false;
    // Number of seconds from midnight on 1 Jan 1970 to midnight on 1 Jan 2000.
    static constexpr chr::seconds y2k_offset{chr::seconds{946684800}};

    static constexpr time_point y2k{kep_clock::time_point{} + y2k_offset};

    static constexpr std::time_t to_time_t(const time_point &t) noexcept
    {
        return static_cast<std::time_t>(chr::duration_cast<chr::seconds>(t.time_since_epoch()).count());
    }

    static constexpr time_point from_time_t(std::time_t t) noexcept
    {
        return chr::time_point_cast<duration>(time_point(chr::seconds(t)));
    }

    static time_point utc_now() noexcept;
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
    enum class string_format { ISO };

    /** Constructors */
    // Default constructor
    epoch();

    // Constructor from a julian date (as a floating-point value)
    explicit epoch(double epoch_in, julian_type epoch_type = julian_type::MJD2000);

    // Conxtructor from string
    explicit epoch(const std::string &, string_format = epoch::string_format::ISO);

    // Constructor for const time_point&)
    explicit epoch(const kep_clock::time_point &time_point);

    // Constructor for const time_point&&)
    explicit epoch(kep_clock::time_point &&time_point);

    /**
     * Constructs an epoch from a std::chrono::duration.
     * The reference point is assumed to be MJD2000.
     * \param[in] time The time as a duration.
     */
    template <Duration D>
    explicit epoch(const D &duration) : m_tp{kep_clock::time_point{} + duration}
    {
    }

    // Constructor from duration&&
    //template <Duration D>
    //explicit epoch(D &&duration) : m_tp{kep_clock::time_point{} + std::forward(duration)}
    //{
    //}

    // Constructor for datetime broken down into its constituents.
    explicit epoch(std::int32_t y, std::uint32_t mon, std::uint32_t d, std::int32_t h = 0, std::int32_t min = 0, std::int32_t s = 0, std::int32_t ms = 0, std::int32_t us = 0);

    /* Computing non-Gregorian dates */

    /**
     * @return Number of days since 0 JD (including fractional days).
     */
    [[nodiscard]] constexpr double jd() const
    {
        return chr::duration<double, std::ratio<86400>>(m_tp.time_since_epoch() - kep_clock::y2k_offset
                                                        + chr::seconds{211813444800})
            .count();
    }

    /**
     * @return Number of days since 0 MJD (including fractional days).
     */
    [[nodiscard]] constexpr double mjd() const
    {
        return chr::duration<double, std::ratio<86400>>(m_tp.time_since_epoch() - kep_clock::y2k_offset
                                                        + chr::seconds{4453401600})
            .count();
    }

    /**
     * @return Number of days since 0 MJD2000 (including fractional days).
     */
    [[nodiscard]] constexpr double mjd2000() const
    {
        return chr::duration<double, std::ratio<86400>>(m_tp.time_since_epoch() - kep_clock::y2k_offset).count();
    }

    /* Helper functions for constructors */
    static kep_clock::time_point make_tp(std::int32_t y, std::uint32_t mon, std::uint32_t d, std::int32_t h = 0, std::int32_t min = 0, std::int32_t s = 0, std::int32_t ms = 0, std::int32_t us = 0);

    static kep_clock::time_point make_tp(double epoch_in, julian_type epoch_type);

    // Conversions
    static constexpr kep_clock::time_point tp_from_days(double days);

    // Duration conversions
    static constexpr double as_sec(kep_clock::duration &&d)
    {
        return chr::duration<double, chr::seconds::period>(d).count();
    }

    // Printing
    [[nodiscard]] std::string as_utc_string() const;
    static std::string as_utc_string(const kep_clock::time_point &tp);

    /** operator overloads for sum and diff (epoch-days) and comparison
     * operators
     * **/

    kep3_DLL_PUBLIC friend std::ostream &operator<<(std::ostream &s, epoch const &epoch_in);

    template <Duration D>
    epoch &operator+=(const D &duration)
    {
        m_tp += chr::duration_cast<kep_clock::duration>(duration);
        return *this;
    }

    template <Duration D>
    epoch &operator-=(const D &duration)
    {
        m_tp -= chr::duration_cast<kep_clock::duration>(duration);
        return *this;
    }

    kep3_DLL_PUBLIC friend bool operator>(const epoch &c1, const epoch &c2);
    kep3_DLL_PUBLIC friend bool operator<(const epoch &c1, const epoch &c2);
    kep3_DLL_PUBLIC friend bool operator>=(const epoch &c1, const epoch &c2);
    kep3_DLL_PUBLIC friend bool operator<=(const epoch &c1, const epoch &c2);
    kep3_DLL_PUBLIC friend bool operator==(const epoch &c1, const epoch &c2);
    kep3_DLL_PUBLIC friend bool operator!=(const epoch &c1, const epoch &c2);

    template <Duration D>
    epoch operator+(const D &duration)
    {
        return epoch(m_tp + chr::duration_cast<kep_clock::duration>(duration));
    }

    template <Duration D>
    epoch operator-(const D &duration)
    {
        return epoch(m_tp - chr::duration_cast<kep_clock::duration>(duration));
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

    [[nodiscard ]]kep_clock::time_point get_tp() const;

private:
    // Serialization code
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive&ar, const unsigned) const
    {
        const auto count{ m_tp.time_since_epoch().count() };
        ar&count;
    }
    template<class Archive>
    void load(Archive&ar, const unsigned)
    {
        kep3::kep_clock::rep count{0};
        ar&count;
        m_tp = kep_clock::time_point{std::chrono::microseconds(count)};
    }

    template<class Archive>
    void serialize(
        Archive & ar,
        const unsigned int file_version
    ) {
        boost::serialization::split_member(ar, *this, file_version);
    }

    // Time point relative to 1 Jan 2000 (MJD2000)
    kep_clock::time_point m_tp;
};

kep3_DLL_PUBLIC epoch utc_now();

kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &s, const epoch &epoch_in);

} // end of namespace kep3

template <>
struct fmt::formatter<kep3::epoch> : fmt::ostream_formatter {
};


namespace boost::serialization
{
    template<class Archive>
    void save(Archive&ar, const std::chrono::microseconds&us, const unsigned)
    {
        auto rep{reinterpret_cast<kep3::kep_clock::rep>(us.count())};
        ar & rep;
    }
    template<class Archive>
    void load(Archive&ar, std::chrono::microseconds&us, const unsigned)
    {
        kep3::kep_clock::rep rep{0};
        ar & rep;
        us = std::chrono::microseconds{rep};
    }
}  // namespace boost::serialization


BOOST_SERIALIZATION_SPLIT_FREE(std::chrono::microseconds)

#endif // kep3_EPOCH_HPP
