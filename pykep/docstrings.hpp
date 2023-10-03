// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef PYKEP_DOCSTRINGS_HPP
#define PYKEP_DOCSTRINGS_HPP

#include <string>

namespace pykep {
// Modules
std::string core_module_doc();

// Anomaly conversions
std::string m2e_doc();
std::string e2m_doc();
std::string m2f_doc();
std::string f2m_doc();
std::string e2f_doc();
std::string f2e_doc();

std::string n2h_doc();
std::string h2n_doc();
std::string n2f_doc();
std::string f2n_doc();
std::string h2f_doc();
std::string f2h_doc();

std::string zeta2f_doc();
std::string f2zeta_doc();

// Anomaly conversions (vectorized)
std::string m2e_v_doc();
std::string e2m_v_doc();
std::string m2f_v_doc();
std::string f2m_v_doc();
std::string e2f_v_doc();
std::string f2e_v_doc();

std::string n2h_v_doc();
std::string h2n_v_doc();
std::string n2f_v_doc();
std::string f2n_v_doc();
std::string h2f_v_doc();
std::string f2h_v_doc();

std::string zeta2f_v_doc();
std::string f2zeta_v_doc();

} // namespace pykep

#endif