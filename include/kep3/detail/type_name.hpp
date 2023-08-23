// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_DETAIL_TYPE_NAME_HPP
#define kep3_DETAIL_TYPE_NAME_HPP

#include <string>
#include <type_traits>
#include <typeinfo>

#include <kep3/detail/visibility.hpp>

namespace kep3
{

namespace detail
{

kep3_DLL_PUBLIC std::string demangle_from_typeid(const char *);

// Determine the name of the type T at runtime.
template <typename T>
inline std::string type_name()
{
    // Get the demangled name without cvref.
    auto ret
        = demangle_from_typeid(typeid(typename std::remove_cv<typename std::remove_reference<T>::type>::type).name());

    // Redecorate it with cv qualifiers.
    constexpr unsigned flag = unsigned(std::is_const<typename std::remove_reference<T>::type>::value)
                              + (unsigned(std::is_volatile<typename std::remove_reference<T>::type>::value) << 1);
    switch (flag) {
        case 0u:
            // NOTE: handle this explicitly to keep compiler warnings at bay.
            break;
        case 1u:
            ret += " const";
            break;
        case 2u:
            ret += " volatile";
            break;
        case 3u:
            ret += " const volatile";
    }

    // Re-add the reference, if necessary.
    if (std::is_lvalue_reference<T>::value) {
        ret += " &";
    } else if (std::is_rvalue_reference<T>::value) {
        ret += " &&";
    }

    return ret;
}

} // namespace detail

} // namespace kep3

#endif