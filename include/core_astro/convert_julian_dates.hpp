/*****************************************************************************
 *   Copyright (C) 2023 The pykep development team,                     *
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

#ifndef kep3_CONVERT_JULIAN_DATES_H
#define kep3_CONVERT_JULIAN_DATES_H

namespace kep3
{
inline double jd2mjd(double in)
{
    return (in - 2400000.5);
}
inline double jd2mjd2000(double in)
{
    return (in - 2451544.5);
}
inline double mjd2jd(double in)
{
    return (in + 2400000.5);
}
inline double mjd2mjd2000(double in)
{
    return (in - 51544);
}
inline double mjd20002jd(double in)
{
    return (in + 2451544.5);
}
inline double mjd20002mjd(double in)
{
    return (in + 51544);
}
}

#endif // kep3_CONVERT_JULIAN_DATES_H
