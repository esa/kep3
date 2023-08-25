// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_DETAIL_TYPEID_NAME_EXTRACT_HPP
#define kep3_DETAIL_TYPEID_NAME_EXTRACT_HPP

#include <cstring>
#include <type_traits>
#include <typeinfo>

#include <kep3/detail/type_traits.hpp>

namespace kep3::detail {

// This is an implementation of the extract() functionality
// for UDx classes based on the name() of the UDx C++ type,
// as returned by typeid().name(). This is needed
// because the dynamic_cast() used in the
// usual extract() implementations can fail on some
// compiler/platform/stdlib implementations
// when crossing boundaries between dlopened()
// modules. See:
// https://github.com/pybind/pybind11/issues/912#issuecomment-310157016
// https://bugs.llvm.org/show_bug.cgi?id=33542
template <typename T, typename C>
inline typename std::conditional<std::is_const<C>::value, const T *, T *>::type
typeid_name_extract(C &class_inst) {
  // NOTE: typeid() strips away both reference and cv qualifiers. Thus,
  // if T is cv-qualified or a reference type, return nullptr preemptively
  // (in any case, extraction cannot be successful in such cases).
  if (!std::is_same<T, uncvref_t<T>>::value || std::is_reference<T>::value) {
    return nullptr;
  }

  if (std::strcmp(class_inst.get_type_index().name(), typeid(T).name()) != 0) {
    // The names differ, return null.
    return nullptr;
  } else {
    // The names match, cast to the correct type and return.
    return static_cast<typename std::conditional<std::is_const<C>::value,
                                                 const T *, T *>::type>(
        class_inst.get_ptr());
  }
}

} // namespace kep3::detail

#endif
