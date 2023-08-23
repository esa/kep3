// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/date_time/gregorian/gregorian.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <kep3/core_astro/convert_julian_dates.hpp>
#include <kep3/epoch.hpp>

// This is set to the precision of the boost date library (microseconds is
// default, nanoseconds can be set when compiling boosts. Note that the code has
// not been checked in that case)
#define BOOST_DATE_PRECISION 1e-6

namespace kep3 {
using namespace boost::gregorian;
using namespace boost::posix_time;

/// Constructor.
/**
 * Constructs an epoch from a non-gregorian date.
 * \param[in] epoch_in A double indicating the non-gregorian date
 * \param[in] epoch_type One of [epoch::MJD2000, epoch::MJD, epoch::JD]
 */
epoch::epoch(const double &epoch_in, julian_type epoch_type)
    : m_mjd2000(epoch_in) {
  switch (epoch_type) {
  case MJD2000:
    break;
  case MJD:
    m_mjd2000 = mjd2mjd2000(epoch_in);
    break;
  case JD:
    m_mjd2000 = jd2mjd2000(epoch_in);
    break;
  }
}

/// Constructor.
/**
 * Constructs an epoch from a gregorian date. The time of the day is assumed to
 * be midnight. \param[in] year The gregorian year \param[in] month The month of
 * the year \param[in] day The day of the month
 */
epoch::epoch(const greg_day &day, const greg_month &month,
             const greg_year &year) {
  set_posix_time(ptime(date(year, month, day)));
}

/// Constructor.
/**
 * Constructs an epoch from a boost ptime object (posix time)
 * \param[in] posix_time The posix_time
 */
epoch::epoch(const boost::posix_time::ptime &posix_time) {
  time_duration dt = posix_time - ptime(date(2000, 1, 1));
  bool flag = false;
  if (dt.is_negative()) {
    flag = true;
    dt = dt.invert_sign();
  }
  double fr_secs =
      static_cast<double>(dt.fractional_seconds()) * BOOST_DATE_PRECISION;
  m_mjd2000 = static_cast<double>(dt.hours()) / 24.0 +
              static_cast<double>(dt.minutes()) / 1440.0 +
              (static_cast<double>(dt.seconds()) + fr_secs) / 86400.0;
  if (flag)
    m_mjd2000 = -m_mjd2000;
}

/// jd getter.
/**
 * Returns the julian date
 *
 * @return double containing the julian date
 *
 */
double epoch::jd() const { return mjd20002jd(m_mjd2000); }

/// mjd getter.
/**
 * Returns the modified julian date
 *
 * @return double containing the modified julian date
 *
 */
double epoch::mjd() const { return mjd20002mjd(m_mjd2000); }

/// mjd2000 getter.
/**
 * Gets the modified julian date 2000
 * @return const reference to mjd2000
 */
double epoch::mjd2000() const { return m_mjd2000; }

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
ptime epoch::get_posix_time() const {
  long hrs, min, sec, fsec;
  bool flag = false;
  double copy = m_mjd2000;
  if (copy < 0) {
    copy = -copy;
    flag = true;
  }
  hrs = static_cast<long>(copy * 24);
  min = static_cast<long>((copy * 24 - static_cast<double>(hrs)) * 60);
  sec = static_cast<long>((((copy * 24 - static_cast<double>(hrs)) * 60) -
                           static_cast<double>(min)) *
                          60);
  double dblfsec = ((((copy * 24 - static_cast<double>(hrs)) * 60) -
                     static_cast<double>(min)) *
                    60) -
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
    retval = ptime(date(2000, 1, 1), time_duration(-hrs, -min, -sec, -fsec));
  else
    retval = ptime(date(2000, 1, 1), time_duration(hrs, min, sec, fsec));
  return retval;
}

/// Sets the epoch from a posix time
/**
 * Sets the epoch to a particular posix_time.
 *
 * \param[in] posix_time containing the posix time
 *
 */
void epoch::set_posix_time(const boost::posix_time::ptime &posix_time) {

  m_mjd2000 = epoch(posix_time).mjd2000();
}

epoch &epoch::operator+=(double rhs) {
  /* addition of rhs to *this takes place here */
  m_mjd2000 += rhs;
  return *this; // return the result by reference
}

epoch &epoch::operator-=(double rhs) {
  /* addition of rhs to *this takes place here */
  m_mjd2000 -= rhs;
  return *this; // return the result by reference
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
epoch epoch_from_string(const std::string date) {
  return epoch(
      boost::posix_time::ptime(boost::posix_time::time_from_string(date)));
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
epoch epoch_from_iso_string(const std::string date) {
  return epoch(
      boost::posix_time::ptime(boost::posix_time::from_iso_string(date)));
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
std::ostream &operator<<(std::ostream &s, const epoch &now) {
  s << now.get_posix_time();
  return s;
}
bool operator>(const epoch &c1, const epoch &c2) {
  return (c1.m_mjd2000 > c2.m_mjd2000) ? true : false;
}
bool operator<(const epoch &c1, const epoch &c2) {
  return (c1.m_mjd2000 < c2.m_mjd2000) ? true : false;
}
bool operator>=(const epoch &c1, const epoch &c2) {
  return (c1.m_mjd2000 >= c2.m_mjd2000) ? true : false;
}
bool operator<=(const epoch &c1, const epoch &c2) {
  return (c1.m_mjd2000 <= c2.m_mjd2000) ? true : false;
}
bool operator==(const epoch &c1, const epoch &c2) {
  return (c1.m_mjd2000 == c2.m_mjd2000) ? true : false;
}
bool operator!=(const epoch &c1, const epoch &c2) {
  return (c1.m_mjd2000 != c2.m_mjd2000) ? true : false;
}
epoch operator+(epoch lhs, double rhs) {
  lhs += rhs; // reuse compound assignment
  return lhs; // return the result by value (uses move constructor)
}
epoch operator-(epoch lhs, double rhs) {
  lhs -= rhs; // reuse compound assignment
  return lhs; // return the result by value (uses move constructor)
}
double operator-(const epoch &lhs, const epoch &rhs) {
  return lhs.mjd2000() -
         rhs.mjd2000(); // return the result by value (uses move constructor)
}

} // namespace kep3