// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EPOCH_HPP
#define EPOCH_HPP

#include <iostream>

#include <fmt/ostream.h>

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/lexical_cast.hpp>

#include <chrono>
#include <ctime>

#include <kep3/detail/s11n.hpp>
#include <kep3/detail/visibility.hpp>

/// Keplerian Toolbox
/**
 * This namespace contains astrodynamics and space flight mechanics routines
 * that are related to keplerian motions or models.
 */
namespace kep3
{

    using namespace std::literals;
    namespace chr = std::chrono;

    // struct MJD2KClock
    // {

    //     using Clock = chr::system_clock;
    //     using rep = long long int;
    //     using period = typename Clock::period;
    //     using duration = typename Clock::duration;
    //     using time_point = typename chr::time_point<MJD2KClock>;

    //     static constexpr bool is_steady = Clock::is_steady;
    //     static constexpr typename Clock::time_point ref_epoch =
    //         chr::sys_days{ chr::year(2000) / chr::month(1) / chr::day(1) };

    //     static time_point now() { return time_point(Clock::now() - ref_epoch); }

    //     static time_t to_time_t(const time_point& t)
    //     {
    //         return Clock::to_time_t(ref_epoch + t.time_since_epoch());
    //     }

    //     static time_point from_time_t(time_t t)
    //     {
    //         return time_point(Clock::from_time_t(t) - ref_epoch);
    //     }
    // };
    struct MJD2KClock: public chr::system_clock
    {

        using rep = long long int;
        using period = std::ratio<1>;
        using duration = chr::duration<rep, period>;
        using time_point = chr::time_point<MJD2KClock>;
        static constexpr bool is_steady = false;
        static constexpr chr::seconds unix_diff{ 946684800s };

        static constexpr time_point ref_epoch{ MJD2KClock::time_point{} + unix_diff };

        static std::time_t to_time_t(const time_point& t) noexcept
        {
            return std::time_t(chr::duration_cast<duration>(t.time_since_epoch() + unix_diff).count());
        }

        static time_point from_time_t(std::time_t t) noexcept
        {
            return time_point(chr::seconds(t) - unix_diff);
        }
    };

    using namespace std::literals;
    /// epoch class.
    /**
     * This class defines and contains a non-gregorian date (i.e. a date expressed
     * in julian form). It also provides the user with an interface to boost
     * gregorian dates (see boost documentation at
     * http://www.boost.org/doc/libs/1_42_0/doc/html/date_time.html)
     * using the posix time.
     * The date is defined in MJD2000 (double) as a
     * private member
     *
     * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
     */
    class kep3_DLL_PUBLIC epoch
    {
    public:
        /** Types of non gregorian dates supported. Julian Date (JD) is the number of
         * days passed since January 1, 4713 BC at noon. Modified Julian Date (MJD) is
         * the number of days passed since November 17, 1858 at 00:00 am. The Modified
         * Julian Date 2000 (MJD2000) is the number of days passed since Juanuary 1,
         * 2000 at 00:00am.
         */
        enum class julian_type { MJD2000, MJD, JD };

        /** Constructors */
        explicit epoch(const double epoch_in = 0.0, const julian_type epoch_type = julian_type::MJD2000);
        explicit epoch(const int year, const int month, const int day, const int hour = 0, const int minute = 0, const int second = 0);
        explicit epoch(const boost::posix_time::ptime& posix_time);

        /** Computing non-gregorian dates */
        [[nodiscard]] MJD2KClock::time_point mjd2000() const;
        [[nodiscard]] MJD2KClock::time_point jd() const;
        [[nodiscard]] MJD2KClock::time_point mjd() const;

        /** Interface to boost::posix_time::ptime */
        [[nodiscard]] boost::posix_time::ptime get_posix_time() const;

        void set(const int year, const int month, const int day, const int hour = 0,
                 const int minute = 0, const int second = 0);
        MJD2KClock::time_point make_tm(const int year, const int month, const int day,
                                       const int hour = 0, const int minute = 0, const int second = 0);

        // Calendar conversions
        MJD2KClock::time_point tp_from_seconds(const double seconds) const;
        MJD2KClock::time_point tp_from_days(const double days) const;
        void set_tp_mjd2000(const double epoch_in);
        void set_tp_mjd(const double epoch_in);
        void set_tp_jd(const double epoch_in);

        // Duration conversions
        MJD2KClock::rep as_sec();
        double as_days();

        // Printing
        static const char* as_utc_string(const MJD2KClock::time_point&);
        static std::time_t as_gmtime(const MJD2KClock::time_point&);

        /** operators overloads for sum diff (epoch-days) and the comparison
         * operators
         * **/

        std::chrono::seconds day2sec(const double days) const;
        kep3_DLL_PUBLIC epoch epoch_from_string(const std::string& date);
        kep3_DLL_PUBLIC epoch epoch_from_iso_string(const std::string& date);

        kep3_DLL_PUBLIC friend std::ostream& operator<<(std::ostream& s, epoch const& epoch_in);
        template<typename T, std::chrono::__enable_if_is_duration<T>>
        kep3_DLL_PUBLIC friend epoch& operator+=(epoch& ep, T rhs);
        template<typename T, std::chrono::__enable_if_is_duration<T>>
        kep3_DLL_PUBLIC friend epoch& operator-=(epoch& ep, T rhs);
        kep3_DLL_PUBLIC friend bool operator>(const epoch& c1, const epoch& c2);
        kep3_DLL_PUBLIC friend bool operator<(const epoch& c1, const epoch& c2);
        kep3_DLL_PUBLIC friend bool operator>=(const epoch& c1, const epoch& c2);
        kep3_DLL_PUBLIC friend bool operator<=(const epoch& c1, const epoch& c2);
        kep3_DLL_PUBLIC friend bool operator==(const epoch& c1, const epoch& c2);
        kep3_DLL_PUBLIC friend bool operator!=(const epoch& c1, const epoch& c2);
        template<typename T, std::chrono::__enable_if_is_duration<T>>
        kep3_DLL_PUBLIC friend epoch operator+(const epoch& lhs, T rhs);
        template<typename T, std::chrono::__enable_if_is_duration<T>>
        kep3_DLL_PUBLIC friend epoch operator-(const epoch& lhs, T rhs);
        kep3_DLL_PUBLIC friend MJD2KClock::rep operator-(const epoch lhs, const epoch rhs);

    private:
        // Serialization code
        friend class boost::serialization::access;
        template <class Archive> void serialize(Archive& ar, const unsigned int)
        {
            ar& m_mjd2000;
        }
        // Serialization code (END)

        /// the modified julian date 2000 stored in a double
        double m_mjd2000;
        MJD2KClock::time_point tp;
    };
    // kep3_DLL_PUBLIC bool operator>(const epoch& c1, const epoch& c2);
    // kep3_DLL_PUBLIC bool operator<(const epoch& c1, const epoch& c2);
    // kep3_DLL_PUBLIC bool operator>=(const epoch& c1, const epoch& c2);
    // kep3_DLL_PUBLIC bool operator<=(const epoch& c1, const epoch& c2);
    // kep3_DLL_PUBLIC bool operator==(const epoch& c1, const epoch& c2);
    // kep3_DLL_PUBLIC bool operator!=(const epoch& c1, const epoch& c2);
    // template<typename T, std::chrono::__enable_if_is_duration<T>>
    // kep3_DLL_PUBLIC epoch operator+(const epoch& lhs, T rhs);
    // template<typename T, std::chrono::__enable_if_is_duration<T>>
    // kep3_DLL_PUBLIC epoch operator-(const epoch& lhs, T rhs);
    // kep3_DLL_PUBLIC MJD2KClock::rep operator-(const epoch lhs, const epoch rhs);

} // end of namespace kep3

template <> struct fmt::formatter<kep3::epoch>: ostream_formatter {};

#endif // EPOCH_HPP
