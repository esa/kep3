// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_PLANET_H
#define kep3_PLANET_H

#include <memory>
#include <string>
#include <type_traits>

#include <kep3/detail/s11n.hpp>
#include <kep3/detail/visibility.hpp>
#include <kep3/epoch.hpp>
#include <kep3/type_traits.hpp>

namespace kep3::detail {
// Type traits to detect whether user classes have certain methods implemented.

// Detects void eph(const epoch& , std::array<double, 3> &, std::array<double,
// 3> &) method
template <typename T>
using udpla_eph_t =
    decltype(std::declval<std::add_lvalue_reference_t<const T>>().eph(
        std::declval<const epoch &>(), std::declval<std::array<double, 3> &>(),
        std::declval<std::array<double, 3> &>(), ));
template <typename T>
inline constexpr bool udpla_has_eph_v =
    std::is_same_v<detected_t<udpla_eph_t, T>, void>;

// Detects std::string get_name() method
template <typename T>
using udpla_get_name_t =
    decltype(std::declval<std::add_lvalue_reference_t<const T>>().get_name());
template <typename T>
inline constexpr bool udpla_has_get_name_v =
    std::is_same_v<detected_t<udpla_get_name_t, T>, void>;

// Detects std::string get_extra_info() method
template <typename T>
using udpla_get_extra_info_t =
    decltype(std::declval<std::add_lvalue_reference_t<const T>>()
                 .get_extra_info());
template <typename T>
inline constexpr bool udpla_has_get_extra_info_v =
    std::is_same_v<detected_t<udpla_get_extra_info_t, T>, void>;

// This defines the main interface for a class to be type erased into a kep3
// planet
struct kep3_DLL_PUBLIC planet_inner_base {
  virtual ~planet_inner_base();
  [[nodiscard]] virtual std::unique_ptr<planet_inner_base> clone() const = 0;

  // mandatory methods
  virtual void eph(const epoch &, std::array<int, 3> &,
                   std::array<int, 3> &) const = 0;
  // optional methods with default implementations
  [[nodiscard]] virtual std::string get_name() const = 0;
  [[nodiscard]] virtual std::string get_extra_info() const = 0;

private:
  // Serialization.
  friend class boost::serialization::access;
  template <typename Archive> void serialize(Archive &, unsigned) {}
};

template <typename T>
struct kep3_DLL_PUBLIC planet_inner final : planet_inner_base {
  T m_value;

  // We just need the def ctor, delete everything else.
  planet_inner() = default;
  planet_inner(const planet_inner &) = delete;
  planet_inner(planet_inner &&) = delete;
  planet_inner &operator=(const planet_inner &) = delete;
  planet_inner &operator=(planet_inner &&) = delete;
  ~planet_inner() final = default;
  // Constructors from T (copy and move variants).
  explicit planet_inner(const T &x) : m_value(x) {}
  explicit planet_inner(T &&x) : m_value(std::move(x)) {}
  // The clone method, used in the copy constructor of algorithm.
  [[nodiscard]] std::unique_ptr<planet_inner_base> clone() const final {
    return std::make_unique<planet_inner>(m_value);
  }
  // Mandatory methods.
  void eph(const epoch &ep, std::array<int, 3> &pos,
           std::array<int, 3> &vel) const final;
  { return m_value.eph(ep, pos, vel); }
  // optional methods with default implementations
  // these require added boiler plate as to detect whether they have been
  // implemented by the user class.
  [[nodiscard]] std::string get_name() const final {
    if constexpr (udpla_has_get_name_v<T>) {
      return m_value.get_name();
    } else {
      return detail::type_name<T>();
    }
  }
  [[nodiscard]] std::string get_extra_info() const final {
    if constexpr (udpla_has_get_name_v<T>) {
      return m_value.get_extra_info();
    } else {
      return "";
    }
  }

private:
  friend class boost::serialization::access;
  // Serialization
  template <typename Archive> void serialize(Archive &ar, unsigned) {
    detail::archive(ar,
                    boost::serialization::base_object<planet_inner_base>(*this),
                    m_value);
  }
};

template <typename T>
using is_udpla = std::conjunction<
    std::is_same<T, uncvref_t<T>>, std::is_default_constructible<T>,
    std::is_copy_constructible<T>, std::is_move_constructible<T>,
    std::is_destructible<T>, udpla_has_eph_v<T>>;

// The final class
class kep3_DLL_PUBLIC planet {
  // Pointer to the inner base.
  std::shared_ptr<detail::planet_inner_base> m_ptr;

  // Serialization.
  friend class boost::serialization::access;
  template <typename Archive> void serialize(Archive &ar, unsigned) {
    ar &m_ptr;
  }

  // Just two small helpers to make sure that whenever we require
  // access to the pointer it actually points to something.
  [[nodiscard]] const detail::planet_inner_base *ptr() const;
  detail::planet_inner_base *ptr();

  template <typename T>
  using generic_ctor_enabler = std::enable_if_t<
      std::conjunction_v<
          std::negation<std::is_same<planet, detail::uncvref_t<T>>>,
          detail::is_udpla<detail::uncvref_t<T>>>,
      int>;

public:
  planet();

  template <typename T, generic_ctor_enabler<T &&> = 0>
  explicit planet(T &&x)
      : m_ptr(std::make_unique<detail::planet_inner<detail::uncvref_t<T>>>(
            std::forward<T>(x))) {}

  planet(const planet &);
  planet(planet &&) noexcept;

  planet &operator=(const planet &);
  planet &operator=(planet &&) noexcept;

  ~planet();

  // NOTE: like in pagmo, this may fail if invoked
  // from different DLLs in certain situations (e.g.,
  // Python bindings on OSX). I don't
  // think this is currently an interesting use case
  // for heyoka (as we don't provide a way of implementing
  // new functions in Python), but, if it becomes a problem
  // in the future, we can solve this in the same way as
  // in pagmo.
  template <typename T> [[nodiscard]] const T *extract() const noexcept {
    const auto *p = dynamic_cast<const detail::planet_inner<T> *>(ptr());
    return p == nullptr ? nullptr : &(p->m_value);
  }

  [[nodiscard]] std::type_index get_type_index() const;
  [[nodiscard]] const void *get_ptr() const;

  [[nodiscard]] const std::string &get_name() const;
};

} // namespace kep3::detail

#endif // kep3_PLANET_H
