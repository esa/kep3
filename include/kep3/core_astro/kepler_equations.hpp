// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
