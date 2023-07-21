/*****************************************************************************
 *   Copyright (C) 2023 The pykep development team,                          *
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

#ifndef kep3_kep3LER_EQUATIONS_H
#define kep3_kep3LER_EQUATIONS_H

#include <cmath>

namespace kep3 {
// In terms of the eccentric anomaly (E)
// -------------------------------------------
inline double kepE(double E, double M, double ecc) {
  return (E - ecc * std::sin(E) - M);
}
// Its first derivative
inline double d_kepE(double E, double ecc) { return (1 - ecc * std::cos(E)); }
// And its second derivative
inline double dd_kepE(double E, double ecc) { return ecc * std::sin(E); }
// -------------------------------------------

// In terms of the hyperbolic anomaly (E)
// -------------------------------------------
inline double kepH(double H, double M, double ecc) {
  return (ecc * std::sinh(H) - H - M);
}
// Its first derivative
inline double d_kepH(double H, double ecc) { return (ecc * std::cosh(H) - 1); }
// And its second derivative
inline double dd_kepH(double H, double ecc) { return ecc * std::sinh(H); }
// -------------------------------------------

// In terms of the eccentric anomaly difference (DE)
// -------------------------------------------
inline double kepDE(double &DE, double DM, double sigma0, double sqrta,
                    double a, double R) {
  return ((-DM + DE + sigma0 / sqrta * (1 - std::cos(DE)) -
           (1 - R / a) * std::sin(DE)));
}

inline double d_kepDE(double DE, double sigma0, double sqrta, double a,
                      double R) {
  return ((1 + sigma0 / sqrta * std::sin(DE) - (1 - R / a) * std::cos(DE)));
}

// In terms of the hyperbolic anomaly difference (DH)
// -------------------------------------------
inline double kepDH(double DH, double DN, double sigma0, double sqrta, double a,
                    double R) {
  return (-DN - DH + sigma0 / sqrta * (std::cosh(DH) - 1) +
          (1 - R / a) * std::sinh(DH));
}

inline double d_kepDH(double DH, double sigma0, double sqrta, double a,
                      double R) {
  return (-1. + sigma0 / sqrta * std::sinh(DH) + (1 - R / a) * std::cosh(DH));
}
} // namespace kep3
#endif // kep3_kep3LER_EQUATIONS_H
