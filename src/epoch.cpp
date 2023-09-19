// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/date_time/gregorian/gregorian.hpp>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>

#include "kep3/core_astro/convert_julian_dates.hpp"
#include "kep3/epoch.hpp"

// This is set to the precision of the boost date library (microseconds is
// default, nanoseconds can be set when compiling boosts. Note that the code has
// not been checked in that case)
#define BOOST_DATE_PRECISION 1e-6

namespace kep3
{
    using ::boost::gregorian::date;
    using ::boost::gregorian::greg_day;
    using ::boost::gregorian::greg_month;
    using ::boost::gregorian::greg_year;
    using ::boost::posix_time::from_iso_string;
    using ::boost::posix_time::ptime;
    using ::boost::posix_time::time_duration;
    using ::boost::posix_time::time_from_string;

    using namespace std::literals;

    /// Constructor.
    /**
     * Constructs an epoch from a non-gregorian date.
     * \param[in] epoch_in A double indicating the non-gregorian date
     * \param[in] epoch_type One of [epoch::julian_type::MJD2000,
     * epoch::julian_type::MJD, epoch::julian_type::JD]
     */
    epoch::epoch(const double epoch_in, const julian_type epoch_type)
        :
        m_mjd2000(epoch_in)
    {
        switch (epoch_type)
        {
        case julian_type::MJD2000:
            set_tp_mjd2000(epoch_in);
            break;
        case julian_type::MJD:
            set_tp_mjd(epoch_in);
            break;
        case julian_type::JD:
            set_tp_jd(epoch_in);
            break;
        default:
            throw;
        }
    }

    /// Constructor.
    /**
     * Constructs an epoch from a gregorian date. The time of the day is assumed to
     * be midnight. \param[in] year The gregorian year \param[in] month The month of
     * the year \param[in] day The day of the month
     */
    epoch::epoch(const int year, const int month, const int day,
                 const int hour, const int minute, const int second)
        :
        m_mjd2000(0),
        tp(make_tm(year, month, day, hour, minute, second))
    {
    }

    /// Constructor.
    /**
     * Constructs an epoch from a boost ptime object (posix time)
     * \param[in] posix_time The posix_time
     */
    epoch::epoch(const boost::posix_time::ptime& posix_time)
    {

        time_duration dt = posix_time - ptime(date(2000, 1, 1));
        bool flag = false;
        if (dt.is_negative())
        {
            flag = true;
            dt = dt.invert_sign();
        }
        const double fr_secs =
            static_cast<double>(dt.fractional_seconds()) * BOOST_DATE_PRECISION;
        m_mjd2000 = static_cast<double>(dt.hours()) / 24.0 +
            static_cast<double>(dt.minutes()) / 1440.0 +
            (static_cast<double>(dt.seconds()) + fr_secs) / 86400.0;

        // NOTE: Why not this?
        // m_mjd2000 = dt.total_seconds() / 86400.0;
        if (flag)
        {
            m_mjd2000 = -m_mjd2000;
        }

        tp = tp_from_seconds(flag ? -dt.total_seconds() : dt.total_seconds());
        // std::cout << "Posix time: " << posix_time << "\n";
        // std::cout << "Time point : " << epoch::as_utc_string(tp) << "\n";
    }

    MJD2KClock::time_point epoch::make_tm(const int year, const int month, const int day,
                                          const int hr, const int min, const int sec)
    {
        std::tm t{};
        t.tm_year = year - 1900;
        t.tm_mon = month - 1;
        t.tm_mday = day;
        t.tm_hour = hr;
        t.tm_min = min;
        t.tm_sec = sec;
        auto tm{ std::mktime(&t) };
        auto utctm{ std::gmtime(&tm) };

        return MJD2KClock::from_time_t(std::mktime(utctm));
    }

    /// jd getter.
    /**
     * Returns the julian date
     *
     * @return double containing the julian date
     *
     */
     // double epoch::jd() const { return mjd20002jd(m_mjd2000); }
    MJD2KClock::time_point epoch::jd() const { return mjd20002jd(tp); }

    /// mjd getter.
    /**
     * Returns the modified julian date
     *
     * @return double containing the modified julian date
     *
     */
    MJD2KClock::time_point epoch::mjd() const { return mjd20002mjd(tp); }

    /// mjd2000 getter.
    /**
     * Gets the modified julian date 2000
     * @return const reference to mjd2000
     */
    MJD2KClock::time_point epoch::mjd2000() const { return tp; }

    /// Extracts the posix time
    /**
     * Returns the posix_time representation of the epoch. The method evaluates
     * from the mjd2000 the number of days, months, seconds and
     * micro/nano seconds passed since the 1st of January 2000 and uses this
     * information to build the posix time
     *
     * @return ptime containing the posix time
     *
     */
    ptime epoch::get_posix_time() const
    {
        long hrs = 0., min = 0., sec = 0., fsec = 0.;
        bool flag = false;
        double copy = m_mjd2000;
        if (copy < 0)
        {
            copy = -copy;
            flag = true;
        }
        hrs = static_cast<long>(copy * 24);
        min = static_cast<long>((copy * 24 - static_cast<double>(hrs)) * 60);
        sec = static_cast<long>((((copy * 24 - static_cast<double>(hrs)) * 60) -
            static_cast<double>(min)) * 60);
        const double dblfsec = ((((copy * 24 - static_cast<double>(hrs)) * 60) -
            static_cast<double>(min)) * 60) -
            static_cast<double>(sec);
        std::ostringstream fsecstr;
        fsecstr << std::setiosflags(std::ios::fixed)
            << std::setprecision(
                         -static_cast<int>(std::log10(BOOST_DATE_PRECISION)))
            << dblfsec;
        fsec = boost::lexical_cast<long>(fsecstr.str().substr(
            2, static_cast<unsigned long>(-std::log10(BOOST_DATE_PRECISION) + 1)));
        ptime retval;
        if (flag)
        {
            retval = ptime(date(2000, 1, 1), time_duration(-hrs, -min, -sec, -fsec));
        }
        else
        {
            retval = ptime(date(2000, 1, 1), time_duration(hrs, min, sec, fsec));
        }

        ptime temp_pt = boost::posix_time::from_time_t(MJD2KClock::to_time_t(tp));

        return temp_pt;
    }

    /// Sets the epoch from a posix time
    /**
     * Sets the epoch to a particular posix_time.
     *
     * \param[in] posix_time containing the posix time
     *
     */
    void epoch::set(const int year, const int month, const int day,
                    const int hour, const int minute, const int second)
    {
        tp = make_tm(year, month, day, hour, minute, second);
    }

    MJD2KClock::time_point epoch::tp_from_seconds(const double seconds) const
    {
        return MJD2KClock::time_point{} + chr::round<chr::seconds>(chr::duration<double, std::ratio<1>>(seconds));
    }

    MJD2KClock::time_point epoch::tp_from_days(const double days) const
    {
        return MJD2KClock::time_point{} + day2sec(days);
    }

    std::chrono::seconds epoch::day2sec(const double days) const
    {
        return chr::round<chr::seconds>(chr::duration<double, std::ratio<86400>>(days));
    }

    void epoch::set_tp_mjd2000(const double days)
    {
        tp = tp_from_days(days);
        std::cout << days << " days in MJD2000: " << epoch::as_utc_string(tp);
    }
    void epoch::set_tp_mjd(const double days)
    {
        tp = mjd2mjd2000(tp_from_days(days));
        std::cout << days << " days in MJD: " << epoch::as_utc_string(tp);
    }
    void epoch::set_tp_jd(const double days)
    {
        tp = jd2mjd2000(tp_from_days(days));
        std::cout << days << " in JD: " << epoch::as_utc_string(tp);
    }

    MJD2KClock::rep epoch::as_sec()
    {
        return tp.time_since_epoch().count();
    }
    double epoch::as_days()
    {
        return std::chrono::duration<double, std::ratio<86400>>(tp.time_since_epoch()).count();
    }

    std::time_t epoch::as_gmtime(const MJD2KClock::time_point& tp)
    {
        auto tt{ MJD2KClock::to_time_t(tp) };
        auto utc = std::gmtime(&tt);
        return std::mktime(utc);
    }

    const char* epoch::as_utc_string(const MJD2KClock::time_point& tp)
    {
        auto tm = epoch::as_gmtime(tp);
        return std::ctime(&tm);
    }

    template<typename T, std::chrono::__enable_if_is_duration<T>>
    epoch& operator+=(epoch& ep, T rhs)
    {
        /* addition of rhs to *this takes place here */
        ep.tp += rhs;
        return ep; // return the result by reference
    }

    template<typename T, std::chrono::__enable_if_is_duration<T>>
    epoch& operator-=(epoch& ep, T rhs)
    {
        /* addition of rhs to *this takes place here */
        ep.tp -= rhs;
        return ep; // return the result by reference
    }

    /// Returns an epoch constructed from a delimited string containing a date
    /**
     *  Builds an epoch from a delimited string. Excess digits in fractional seconds
     * will be dropped. Ex: "1:02:03.123456999" => "1:02:03.123456". This behavior
     * depends on the precision defined in astro_constant.h used to compile
     *
     * Example:
     * 	std::string ts("2002-01-20 23:59:54.003");
     * 	epoch e(epoch_from_string(ts))
     *
     */
    epoch epoch_from_string(const std::string& date)
    {
        return epoch(time_from_string(date));
    }

    /// Returns an epoch constructed from a non delimited iso string containing a
    /// date
    /**
     *  Builds an epoch from a non delimited iso string containing a date.
     *
     * Example:
     * 	std::string ts("20020131T235959");
     * 	epoch e(epoch_from_iso_string(ts))
     *
     */
    epoch epoch_from_iso_string(const std::string& date)
    {
        return epoch(from_iso_string(date));
    }

    /// Overload the stream operator for kep_toolbox::epoch
    /**
     * Streams out a date in the format 2000-Jan-01 00:12:30.123457
     *
     * \param[in] s stream to which the epoch will be sent
     * \param[in] now epoch to be sent to stream
     *
     * \return reference to s
     *
     */
    std::ostream& operator<<(std::ostream& s, const epoch& ep)
    {
        s << epoch::as_utc_string(ep.tp);
        return s;
    }
    bool operator>(const epoch& c1, const epoch& c2)
    {
        return c1.tp > c2.tp;
    }
    bool operator<(const epoch& c1, const epoch& c2)
    {
        return c1.tp < c2.tp;
    }
    bool operator>=(const epoch& c1, const epoch& c2)
    {
        return c1.tp >= c2.tp;
    }
    bool operator<=(const epoch& c1, const epoch& c2)
    {
        return c1.tp <= c2.tp;
    }
    bool operator==(const epoch& c1, const epoch& c2)
    {
        return c1.tp == c2.tp;
    }
    bool operator!=(const epoch& c1, const epoch& c2)
    {
        return c1.tp != c2.tp;
    }

    // template<typename T, std::chrono::__enable_if_is_duration<T>>
    // epoch operator+(const epoch& lhs, T rhs)
    // {
    //     lhs += rhs; // reuse compound assignment
    //     return lhs; // return the result by value (uses move constructor)
    // }
    // template<typename T, std::chrono::__enable_if_is_duration<T>>
    // epoch operator-(const epoch& lhs, T rhs)
    // {
    //     lhs -= rhs; // reuse compound assignment
    //     return lhs; // return the result by value (uses move constructor)
    // }
    MJD2KClock::rep operator-(const epoch lhs, const epoch rhs)
    {
        return (lhs - rhs); // return the result by value (uses move constructor)
    }

} // namespace kep3
