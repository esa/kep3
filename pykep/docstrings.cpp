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
    
    Converts from mean to eccentric anomaly. Requires ecc < 1.

    Args:
      M (float): the Mean anomaly (rad.)
      ecc (float): the eccentricity

    Returns:
      float: the Eccentric anomaly

    Examples:
      >>> import pykep as pk
      >>> M = 1.2
      >>> ecc = 0.1
      >>> pk.m2e(M, ecc)
      1.1929459841174401
)";
}
} // namespace pykep