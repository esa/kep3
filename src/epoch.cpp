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

namespace kep3
{

    /// Constructor.
    /**
    * Constructs a reference epoch.
    */
    epoch::epoch()
        : tp{}
    {}

    /// Constructor.
    /**
    * Constructs an epoch from a non-gregorian date.
    * \param[in] epoch_in A double indicating the non-gregorian date
    * \param[in] epoch_type epoch::julian_type
    */
    epoch::epoch(const double epoch_in, const julian_type epoch_type)
      : tp{ make_tp( epoch_in, epoch_type ) }

    {
    }

    /// Constructor.
    /**
     * Constructs an epoch from a gregorian date. The time of the day is assumed to
     * be midnight. \param[in] yr The gregorian year \param[in] mon The month of the
     * year [0-11] \param[in] day The day of the month [1-31] \param[in] hr The hour
     * of the day [0-23] \param[in] min The minutes [0-59] \param[in] s The seconds
     * [0-59] \param[in] ms The milliseconds [0-999] \param[in] us The milliseconds
     * [0-999]
     */
    epoch::epoch(const int y, const int d, const int h, const int min, const int s, const int ms, const int us)
        : tp{ make_tp( y, d, h, min, s, ms, us ) }
    {
    }


    epoch::epoch( const kep_clock::time_point& time_point )
        : tp{ time_point }
    {
    }

    epoch::epoch( kep_clock::time_point&& time_point )
        : tp{ std::move( time_point ) }
    {
    }

    kep_clock::time_point
    epoch::make_tp(const int y, const int d, const int h, const int min,
                   const int s, const int ms, const int us)

    {
        return kep_clock::time_point{} + chr::years( y ) + chr::days( d ) + chr::hours( h ) + chr::minutes( min ) + chr::seconds( s ) +
               chr::milliseconds( ms ) + chr::microseconds( us );
    }

    kep_clock::time_point epoch::make_tp( const double epoch_in, const julian_type epoch_type )
    {
        switch ( epoch_type )
        {
            case julian_type::MJD2000:
                return epoch::tp_from_days( epoch_in );
            case julian_type::MJD:
                return mjd2mjd2000( epoch::tp_from_days( epoch_in ) );
            case julian_type::JD:
                return jd2mjd2000( epoch::tp_from_days( epoch_in ) );
            default:
                throw;
        }
    }
    // /// jd getter.
    // /**
    //  * Returns the julian date
    //  *
    //  * @return double containing the julian date
    //  *
    //  */
    // constexpr double epoch::jd() const
    // {
    //     return chr::duration<double, std::ratio<86400>>(tp.time_since_epoch() + 211813444800s).count();
    // }

    // /// mjd getter.
    // /**
    //  * Returns the modified julian date
    //  *
    //  * @return double containing the modified julian date
    //  *
    //  */
    // constexpr double epoch::mjd() const
    // {
    //     return chr::duration<double, std::ratio<86400>>(tp.time_since_epoch() + 4453401600s).count();
    // }

    // /// mjd2000 getter.
    // /**
    //  * Gets the modified julian date 2000
    //  * @return const reference to mjd2000
    //  */
    // constexpr double epoch::mjd2000() const
    // {
    //     return chr::duration<double, std::ratio<86400>>(tp.time_since_epoch()).count();
    // }

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

    // kep_clock::duration epoch::tp_from_seconds(const double seconds) const {
    //   return chr::duration_cast<kep_clock::duration>(chr::seconds(seconds));
    // }

    constexpr kep_clock::time_point epoch::tp_from_days( const double days )
    {
        return kep_clock::time_point{} +
               chr::duration_cast<kep_clock::duration>(
                   chr::duration<double, std::ratio<86400>>( days ) );
    }

    //    void epoch::set_tp_mjd2000( const double days )
    //    {
    //        //        tp = tp_from_days( days );
    //        //        std::cout << days << " days in MJD2000: " << epoch::as_utc_string( tp );
    //    }
    //    void epoch::set_tp_mjd( const double days )
    //    {
    //        //        tp = mjd2mjd2000( tp_from_days( days ) );
    //        //        std::cout << days << " days in MJD: " << epoch::as_utc_string( tp );
    //    }
    //    void epoch::set_tp_jd( const double days )
    //    {
    //        //        tp = jd2mjd2000( tp_from_days( days ) );
    //        //        std::cout << days << " in JD: " << epoch::as_utc_string( tp );
    //    }

    // double epoch::as_days() {
    //   return std::chrono::duration<double, std::seconds::ratio>(
    //              tp.time_since_epoch())
    //       .count();
    // }

    // std::time_t epoch::as_gmtime(const kep_clock::time_point &tp) {
    //   auto tt{kep_clock::to_time_t(tp)};
    //   auto utc = std::gmtime(&tt);
    //   return std::mktime(utc);
    // }

    auto epoch::as_utc_string( const kep_clock::time_point& tp )
    {
        auto t = kep_clock::to_time_t( tp );
        return std::put_time( gmtime( &t ), "%FT%T");
    }

    //    template <llint Num, llint Den>
    //    epoch& epoch::operator+=( dur<Num, Den>&& rhs )
    //    {
    //        /* addition of rhs to *this takes place here */
    //        tp += chr::duration_cast<kep_clock::duration>( rhs );
    //        return *this;
    //    }

    //    template <llint Num, llint Den>
    //    epoch& epoch::operator-=( dur<Num, Den>&& rhs )
    //    {
    //        /* addition of rhs to *this takes place here */
    //        tp -= chr::duration_cast<kep_clock::duration>( rhs );
    //        return *this;
    //    }

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
    // epoch epoch_from_string(const std::string &date) {
    //   return epoch(time_from_string(date));
    // }

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
    // epoch epoch_from_iso_string(const std::string &date) {
    //   return epoch(from_iso_string(date));
    // }

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
    std::ostream& operator<<( std::ostream& s, const epoch& ep )
    {
        s << epoch::as_utc_string(ep.tp);
        return s;
    }

    bool operator>( const epoch& c1, const epoch& c2 )
    {
        return c1.tp > c2.tp;
    }
    bool operator<( const epoch& c1, const epoch& c2 )
    {
        return c1.tp < c2.tp;
    }
    bool operator>=( const epoch& c1, const epoch& c2 )
    {
        return c1.tp >= c2.tp;
    }
    bool operator<=( const epoch& c1, const epoch& c2 )
    {
        return c1.tp <= c2.tp;
    }
    bool operator==( const epoch& c1, const epoch& c2 )
    {
        return c1.tp == c2.tp;
    }
    bool operator!=( const epoch& c1, const epoch& c2 )
    {
        return c1.tp != c2.tp;
    }

    //    template <llint Num, llint Den>
    //    epoch operator+( const epoch& lhs, dur<Num, Den>&& rhs )
    //    {
    //        lhs.tp +=
    //            chr::duration_cast<kep_clock::duration>( rhs ); // reuse compound assignment
    //        return lhs;                                         // return the result by value (uses move constructor)
    //    }
    //    template <llint Num, llint Den>
    //    epoch operator-( const epoch& lhs, dur<Num, Den>&& rhs )
    //    {
    //        lhs.tp -=
    //            chr::duration_cast<kep_clock::duration>( rhs ); // reuse compound assignment
    //        return lhs;                                         // return the result by value (uses move constructor)
    //    }

    kep_clock::duration operator-( const epoch& lhs, const epoch& rhs )
    {
        return lhs.tp - rhs.tp;
    }

} // namespace kep3
