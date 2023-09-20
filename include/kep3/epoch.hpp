// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_EPOCH_H
#define kep3_EPOCH_H

#include <iostream>

#include <fmt/ostream.h>

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/lexical_cast.hpp>

#include <kep3/detail/s11n.hpp>
#include <kep3/detail/visibility.hpp>

/// Keplerian Toolbox
/**
 * This namespace contains astrodynamics and space flight mechanics routines
 * that are related to keplerian motions or models.
 */
namespace kep3 {

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
class kep3_DLL_PUBLIC epoch {
public:
  /** Types of non gregorian dates supported. Julian Date (JD) is the number of
   * days passed since January 1, 4713 BC at noon. Modified Julian Date (MJD) is
   * the number of days passed since November 17, 1858 at 00:00 am. The Modified
   * Julian Date 2000 (MJD2000) is the number of days passed since Juanuary 1,
   * 2000 at 00:00am.
   */
  enum julian_type { MJD2000, MJD, JD };

  /** Constructors */
  explicit epoch(const double & = 0., julian_type = MJD2000);
  epoch(const boost::gregorian::greg_day &day,
        const boost::gregorian::greg_month &month,
        const boost::gregorian::greg_year &year);
  explicit epoch(const boost::posix_time::ptime &posix_time);

  /** Computing non-gregorian dates */
  [[nodiscard]] double mjd2000() const;
  [[nodiscard]] double jd() const;
  [[nodiscard]] double mjd() const;

  /** Interface to boost::posix_time::ptime */
  [[nodiscard]] boost::posix_time::ptime get_posix_time() const;
  void set_posix_time(const boost::posix_time::ptime &);

  /** operators overloads for sum diff (epoch-days) and the comparison operators
   * **/
  epoch &operator+=(double rhs);
  epoch &operator-=(double rhs);
  kep3_DLL_PUBLIC friend bool operator>(const epoch &c1, const epoch &c2);
  kep3_DLL_PUBLIC friend bool operator<(const epoch &c1, const epoch &c2);
  kep3_DLL_PUBLIC friend bool operator>=(const epoch &c1, const epoch &c2);
  kep3_DLL_PUBLIC friend bool operator<=(const epoch &c1, const epoch &c2);
  kep3_DLL_PUBLIC friend bool operator==(const epoch &c1, const epoch &c2);
  kep3_DLL_PUBLIC friend bool operator!=(const epoch &c1, const epoch &c2);
  kep3_DLL_PUBLIC friend epoch operator+(epoch lhs, double rhs);
  kep3_DLL_PUBLIC friend epoch operator-(epoch lhs, double rhs);
  kep3_DLL_PUBLIC friend double operator-(const epoch &lhs, const epoch &rhs);

private:
  // Serialization code
  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, const unsigned int) {
    ar &m_mjd2000;
  }
  // Serialization code (END)

  /// the modified julian date 2000 stored in a double
  double m_mjd2000;
};

kep3_DLL_PUBLIC epoch epoch_from_string(const std::string &date);
kep3_DLL_PUBLIC epoch epoch_from_iso_string(const std::string &date);

kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &s,
                                         const epoch &epoch_in);


} // end of namespace kep3

template <> struct fmt::formatter<kep3::epoch> : ostream_formatter {};

#endif // kep3_EPOCH_H
