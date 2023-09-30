// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <string>

#include "docstrings.hpp"

namespace pykep {

std::string core_module_doc() {
  return R"(core is the Pykep module that contains most of its core routines efficiently coded in c++
)";
}

std::string m2e_doc() {
  return R"(m2e(M, ecc)
    
    Converts from Mean to Eccentric anomaly. Requires ecc < 1.

    Args:
      M (float): the Mean anomaly (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the Eccentric anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> M = 1.2
      >>> ecc = 0.1
      >>> pk.m2e(M, ecc)
      1.1929459841174401
)";
}

std::string e2m_doc() {
  return R"(e2m(E, ecc)
    
    Converts from Eccentric to Mean anomaly. Requires ecc < 1.

    Args:
      E (float): the Eccentric anomaly (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the Mean anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> E = 0.5
      >>> ecc = 0.1
      >>> pk.e2m(E, ecc)
      0.4520574461395797
)";
}

std::string e2f_doc() {
  return R"(e2f(E, ecc)
    
    Converts from eccentric to true anomaly. Requires ecc < 1.

    Args:
      E (float): the Eccentric anomaly (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> E = 0.5
      >>> ecc = 0.1
      >>> pk.e2f(E, ecc)
      0.5502639747136633
)";
}

std::string f2e_doc() {
  return R"(f2e(f, ecc)
    
    Converts from True to Eccentric anomaly. Requires ecc < 1.

    Args:
      f (float): the True anomaly (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the Eccentric anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> f = 1.2
      >>> ecc = 0.1
      >>> pk.f2e(f, ecc)
      1.1082931139529482
)";
}

std::string f2m_doc() {
  return R"(f2m(f, ecc)
    
    Converts from True to Mean anomaly. Requires ecc < 1.

    Args:
      f (float): the True anomaly (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the Mean anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> f = -0.34
      >>> ecc = 0.67
      >>> pk.f2m(f, ecc)
      -0.8380280766377411
)";
}

std::string m2f_doc() {
  return R"(m2f(M, ecc)
    
    Converts from Mean to True anomaly. Requires ecc < 1.

    Args:
      M (float): the Mean anomaly (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> M = 0.32
      >>> ecc = 0.65
      >>> pk.m2f(M, ecc)
      1.4497431281728277
)";
}

std::string h2n_doc() {
  return R"(h2n(H, ecc)
    
    Converts from Hyperbolic to Hyperbolic Mean anomaly. Requires ecc > 1.

    Args:
      H (float): the Hyperbolic anomaly (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the Hyperbolic Mean anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> H = 1.2
      >>> ecc = 10.32
      >>> pk.h2n(H, ecc)
      14.377641187853621
)";
}

std::string n2h_doc() {
  return R"(n2h(N, ecc)
    
    Converts from Hyperbolic Mean to Hyperbolic anomaly. Requires ecc > 1.

    Args:
      N (float): the Hyperbolic Mean anomaly (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the Hyperbolic anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> N = 1.2
      >>> ecc = 10.32
      >>> pk.n2h(N, ecc)
      0.12836469743916526
)";
}

std::string h2f_doc() {
  return R"(h2f(H, ecc)
    
    Converts from Hyperbolic to True anomaly. Requires ecc > 1.

    Args:
      H (float): the Hyperbolic anomaly (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> H = 10.32
      >>> ecc = 4.5
      >>> pk.h2f(H, ecc)
      1.7948251330114304
)";
}

std::string f2h_doc() {
  return R"(f2h(f, ecc)
    
    Converts from True to Hyperbolic anomaly. Requires ecc > 1.

    Args:
      f (float): the True anomaly (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the Hyperbolic anomaly

    Examples:
      >>> import pykep as pk
      >>> f = 1.2
      >>> ecc = 0.1
      >>> pk.f2h(f, ecc)
      1.1929459841174401
)";
}

std::string f2n_doc() {
  return R"(f2n(f, ecc)
    
    Converts from True to Hyperbolic Mean anomaly. Requires ecc > 1.

    Args:
      f (float): the True anomaly (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the Hyperbolic Mean anomaly

    Examples:
      >>> import pykep as pk
      >>> f = 1.2
      >>> ecc = 5.7
      >>> pk.f2n(f, ecc)
      8.421335633880908
)";
}

std::string n2f_doc() {
  return R"(n2f(N, ecc)
    
    Converts from Hyperbolic Mean to True anomaly. Requires ecc > 1.

    Args:
      N (float): the Hyperbolic Mean anomaly (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the True anomaly

    Examples:
      >>> import pykep as pk
      >>> N = 10.32
      >>> ecc = 13.45
      >>> pk.m2e(M, ecc)
      0.7373697968359353
)";
}

std::string zeta2f_doc() {
  return R"(zeta2f(zeta, ecc)
    
    Converts from Gudermannian to True anomaly. Requires ecc > 1.

    See Battin: "An Introduction to the Mathematics and Methods of Astrodynamics" for a 
    definition of zeta and the treatment of the resulting equations.

    Args:
      zeta (float): the Gudermannian (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the True anomaly

    Examples:
      >>> import pykep as pk
      >>> zeta = 8.2
      >>> ecc = 2.2
      >>> pk.zeta2f(zeta, ecc)
      2.3290929552114266
)";
}

std::string f2zeta_doc() {
  return R"(f2zeta(f, ecc)
    
    Converts from True anomaly to Gudermannian. Requires ecc > 1.

    Args:
      f (float): the True anomaly (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the Gudermannian 

    Examples:
      >>> import pykep as pk
      >>> f = 0.5
      >>> ecc = 3.3
      >>> pk.f2zeta(f, ecc)
      0.36923933496389816
)";
}
} // namespace pykep