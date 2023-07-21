/*****************************************************************************
 *   Copyright (C) 2004-2018 The pykep development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://gitter.im/esa/pykep                                             *
 *   https://github.com/esa/pykep                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#ifndef kep3_EPOCH_H
#define kep3_EPOCH_H

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

#include <detail/s11n.hpp>
#include <detail/visibility.hpp>

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
 * To achieve higher performance the date is defined in MJD2000 (double) as a
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
  epoch(const double & = 0., julian_type = MJD2000);
  epoch(const boost::gregorian::greg_day &day,
        const boost::gregorian::greg_month &month,
        const boost::gregorian::greg_year &year);
  epoch(const boost::posix_time::ptime &posix_time);

  /** Computing non-gregorian dates */
  double mjd2000() const;
  double jd() const;
  double mjd() const;

  /** Interface to boost::posix_time::ptime */
  boost::posix_time::ptime get_posix_time() const;
  void set_posix_time(const boost::posix_time::ptime &);

  /** operators overloads for sum diff (epoch-days) and the comparison operators
   * **/
  epoch &operator+=(double rhs);
  epoch &operator-=(double rhs);
  friend epoch operator+(epoch lhs, double rhs);
  friend epoch operator-(epoch lhs, double rhs);
  friend double operator-(const epoch &lhs, const epoch &rhs);
  friend bool operator>(const epoch &c1, const epoch &c2);
  friend bool operator<(const epoch &c1, const epoch &c2);
  friend bool operator>=(const epoch &c1, const epoch &c2);
  friend bool operator<=(const epoch &c1, const epoch &c2);
  friend bool operator==(const epoch &c1, const epoch &c2);
  friend bool operator!=(const epoch &c1, const epoch &c2);

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

kep3_DLL_PUBLIC epoch epoch_from_string(const std::string date);
kep3_DLL_PUBLIC epoch epoch_from_iso_string(const std::string date);

kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &s,
                                         const epoch &epoch_in);
kep3_DLL_PUBLIC bool operator>(const epoch &c1, const epoch &c2);
kep3_DLL_PUBLIC bool operator<(const epoch &c1, const epoch &c2);
kep3_DLL_PUBLIC bool operator>=(const epoch &c1, const epoch &c2);
kep3_DLL_PUBLIC bool operator<=(const epoch &c1, const epoch &c2);
kep3_DLL_PUBLIC bool operator==(const epoch &c1, const epoch &c2);
kep3_DLL_PUBLIC bool operator!=(const epoch &c1, const epoch &c2);
kep3_DLL_PUBLIC epoch operator+(epoch lhs, double rhs);
kep3_DLL_PUBLIC epoch operator-(epoch lhs, double rhs);
kep3_DLL_PUBLIC double operator-(const epoch &lhs, const epoch &rhs);



} // end of namespace kep3

#endif // kep3_EPOCH_H
