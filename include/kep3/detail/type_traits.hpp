// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_DETAIL_TYPE_TRAITS_HPP
#define kep3_DETAIL_TYPE_TRAITS_HPP

#include <initializer_list>
#include <limits>
#include <type_traits>
#include <vector>


namespace kep3
{

namespace detail
{

template <typename T>
using uncvref_t = std::remove_cv_t<std::remove_reference_t<T>>;

template <typename, typename...>
inline constexpr bool always_false_v = false;

// http://en.cppreference.com/w/cpp/experimental/is_detected
template <class Default, class AlwaysVoid, template <class...> class Op, class... Args>
struct detector {
    using value_t = std::false_type;
    using type = Default;
};

template <class Default, template <class...> class Op, class... Args>
struct detector<Default, std::void_t<Op<Args...>>, Op, Args...> {
    using value_t = std::true_type;
    using type = Op<Args...>;
};

// http://en.cppreference.com/w/cpp/experimental/nonesuch
struct nonesuch {
    nonesuch() = delete;
    ~nonesuch() = delete;
    nonesuch(nonesuch const &) = delete;
    nonesuch(nonesuch &&) noexcept = delete;
    void operator=(nonesuch const &) = delete;
    void operator=(nonesuch &&) noexcept = delete;
};

template <template <class...> class Op, class... Args>
using is_detected = typename detector<nonesuch, void, Op, Args...>::value_t;

template <template <class...> class Op, class... Args>
using detected_t = typename detector<nonesuch, void, Op, Args...>::type;

template <template <class...> class Op, class... Args>
inline constexpr bool is_detected_v = is_detected<Op, Args...>::value;

template <typename T>
inline constexpr bool is_supported_fp_v = is_supported_fp<T>::value;

// Detect vector type.
template <typename>
struct is_any_vector : std::false_type {
};

template <typename T>
struct is_any_vector<std::vector<T>> : std::true_type {
};

template <typename T>
inline constexpr bool is_any_vector_v = is_any_vector<T>::value;

// Detect initializer_list type.
template <typename>
struct is_any_ilist : std::false_type {
};

template <typename T>
struct is_any_ilist<std::initializer_list<T>> : std::true_type {
};

template <typename T>
inline constexpr bool is_any_ilist_v = is_any_ilist<T>::value;


} // namespace detail

} // namespace kep3

#endif