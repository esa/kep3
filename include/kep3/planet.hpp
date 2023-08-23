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
#include <typeindex>

#include <kep3/detail/s11n.hpp>
#include <kep3/detail/type_name.hpp>
#include <kep3/detail/type_traits.hpp>
#include <kep3/detail/visibility.hpp>
#include <kep3/epoch.hpp>

namespace kep3::detail {
// Type traits to detect whether user classes have certain methods implemented.

// udpla_has_eph<T> detects whether T has the method:
// std::array<std::array<double, 3>, 2> eph(const epoch&)
// method
template <typename T>
using udpla_eph_t =
    decltype(std::declval<std::add_lvalue_reference_t<const T>>().eph(
        std::declval<const epoch &>()));

template <typename T>
using udpla_has_eph = std::is_same<detected_t<udpla_eph_t, T>, std::array<std::array<double, 3>, 2>>;

// udpla_has_get_name_v<T> is True if T has the method:
// std::string get_name()
template <typename T>
using udpla_get_name_t =
    decltype(std::declval<std::add_lvalue_reference_t<const T>>().get_name());
template <typename T>
inline constexpr bool udpla_has_get_name_v =
    std::is_same_v<detected_t<udpla_get_name_t, T>, void>;

// udpla_has_get_extra_info_v<T> is True if T has the method:
// std::string get_name()
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
  planet_inner_base(); // why do we need this to be default constructible?
  planet_inner_base(const planet_inner_base &) = delete;
  planet_inner_base(planet_inner_base &&) noexcept = delete;
  planet_inner_base &operator=(const planet_inner_base &) = delete;
  planet_inner_base &operator=(planet_inner_base &&) noexcept = delete;
  virtual ~planet_inner_base();

  [[nodiscard]] virtual std::unique_ptr<planet_inner_base> clone() const = 0;

  [[nodiscard]] virtual std::type_index get_type_index() const = 0;
  [[nodiscard]] virtual const void *get_ptr() const = 0;
  virtual void *get_ptr() = 0;

  // mandatory methods
  [[nodiscard]] virtual std::array<std::array<double, 3>, 2> eph(const epoch &) const = 0;
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

  // Get the type at runtime.
  [[nodiscard]] std::type_index get_type_index() const final {
    return typeid(T);
  }

  // Raw getters for the internal instance.
  [[nodiscard]] const void *get_ptr() const final { return &m_value; }
  void *get_ptr() final { return &m_value; }

  // Mandatory methods.
   [[nodiscard]] std::array<std::array<double, 3>, 2> eph(const epoch &ep) const final {
    return m_value.eph(ep);
  }
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
  // Serialization.
  friend class boost::serialization::access;
  template <typename Archive> void serialize(Archive &ar, unsigned) {
    ar &boost::serialization::base_object<planet_inner_base>(*this);
    ar &m_value;
  }
};
template <typename T>
using is_udpla = std::conjunction<
    std::is_same<T, detail::uncvref_t<T>>, std::is_default_constructible<T>,
    std::is_copy_constructible<T>, std::is_move_constructible<T>,
    std::is_destructible<T>, udpla_has_eph<T>>;

struct null_udpla {
  null_udpla();
  static std::array<std::array<double, 3>, 2> eph(const epoch&);
};
} // namespace kep3::detail

namespace kep3 {
// The final class
class kep3_DLL_PUBLIC planet {
  // Pointer to the inner base.
  std::shared_ptr<detail::planet_inner_base> m_ptr;

  // Serialization.
  friend class boost::serialization::access;
  template <typename Archive> void serialize(Archive &ar, unsigned) {
    ar &m_ptr;
  }

  // Just two small helpers to make sure via assetions that whenever we require
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
  // Default constructor
  planet();
  // Constructor from the UDPLA
  template <typename T, generic_ctor_enabler<T &&> = 0>
  explicit planet(T &&x)
      : m_ptr(std::make_unique<detail::planet_inner<detail::uncvref_t<T>>>(
            std::forward<T>(x))) {}
  // Copy constructor
  planet(const planet &);
  // Move ctor.
  planet(planet &&) noexcept;
  // Move assignment.
  planet &operator=(planet &&) noexcept;
  // Copy assignment.
  planet &operator=(const planet &);
  // Default destructor
  ~planet() = default;
  // Assignment from a user-defined planet of type \p T
  template <typename T, generic_ctor_enabler<T> = 0> planet &operator=(T &&x) {
    return (*this) = planet(std::forward<T>(x));
  }
  // Extract a const pointer to the UDPLA.
  template <typename T> const T *extract() const noexcept {
#if defined(kep3_PREFER_TYPEID_NAME_EXTRACT)
    return detail::typeid_name_extract<T>(*this);
#else
    auto p = dynamic_cast<const detail::planet_inner<T> *>(ptr());
    return p == nullptr ? nullptr : &(p->m_value);
#endif
  }
  /// Extract a pointer to the UDA.
  template <typename T> T *extract() noexcept {
#if defined(kep3_PREFER_TYPEID_NAME_EXTRACT)
    return detail::typeid_name_extract<T>(*this);
#else
    auto p = dynamic_cast<detail::planet_inner<T> *>(ptr());
    return p == nullptr ? nullptr : &(p->m_value);
#endif
  }
  // Checks the user-defined algorithm type at run-time.
  template <typename T> [[nodiscard]] bool is() const noexcept {
    return extract<T>() != nullptr;
  }
  // Check if the planet is valid (i.e. has not been moved from)
  [[nodiscard]] bool is_valid() const;
  // Gets the type of the UDPLA at runtime.
  [[nodiscard]] std::type_index get_type_index() const;
  /// Gets a const pointer to the UDPLA.
  [[nodiscard]] const void *get_ptr() const;
  // Gets a mutable pointer to the UDPLA.
  void *get_ptr();

  std::array<std::array<double, 3>, 2> eph(const epoch &);
  [[nodiscard]] std::string get_name() const;
  [[nodiscard]] std::string get_extra_info() const;
};

} // namespace kep3

#endif // kep3_PLANET_H
