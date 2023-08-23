// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_CONFIG_HPP
#define kep3_CONFIG_HPP

// NOTE: include this so that we can
// detect _LIBCPP_VERSION below.
#include <ciso646>

// Start of defines instantiated by CMake.
// clang-format off
#define kep3_VERSION "0.0.1"
#define kep3_VERSION_MAJOR 0
#define kep3_VERSION_MINOR 0
#define kep3_VERSION_PATCH 1

// clang-format on
// End of defines instantiated by CMake.

#if defined(__clang__) && defined(_LIBCPP_VERSION)

// When using clang + libc++, prefer the name-based
// extract() implementation for UDx classes. See
// the explanation in typeid_name_extract.hpp.
#define kep3_PREFER_TYPEID_NAME_EXTRACT

#endif

#endif
