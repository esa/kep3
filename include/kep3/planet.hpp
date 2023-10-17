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

#include <concepts>
#include <string>
#include <typeinfo>

#include <boost/core/demangle.hpp>

#include <fmt/core.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/detail/s11n.hpp>
#include <kep3/detail/visibility.hpp>
#include <kep3/epoch.hpp>
#include <kep3/exceptions.hpp>

#define TANUKI_WITH_BOOST_S11N

#include <kep3/detail/tanuki.hpp>

#undef TANUKI_WITH_BOOST_S11N

namespace kep3::detail
{

// Concepts to detect whether user classes have certain methods implemented.
#define KEP3_UDPLA_CONCEPT_HAS_GET(name, type)                                                                         \
    template <typename T>                                                                                              \
    concept udpla_has_get_##name = requires(const T &p) {                                                              \
        {                                                                                                              \
            p.get_##name()                                                                                             \
        } -> std::same_as<type>;                                                                                       \
    }

KEP3_UDPLA_CONCEPT_HAS_GET(mu_central_body, double);
KEP3_UDPLA_CONCEPT_HAS_GET(mu_self, double);
KEP3_UDPLA_CONCEPT_HAS_GET(radius, double);
KEP3_UDPLA_CONCEPT_HAS_GET(safe_radius, double);
KEP3_UDPLA_CONCEPT_HAS_GET(name, std::string);
KEP3_UDPLA_CONCEPT_HAS_GET(extra_info, std::string);

#undef KEP3_UDPLA_CONCEPT_HAS_GET

template <typename T>
concept udpla_has_eph = requires(const T &p, const epoch &e) {
    {
        p.eph(e)
    } -> std::same_as<std::array<std::array<double, 3>, 2>>;
};

template <typename T>
concept udpla_has_period = requires(const T &p, const epoch &e) {
    {
        p.period(e)
    } -> std::same_as<double>;
};

template <typename T>
concept udpla_has_elements = requires(const T &p, const epoch &e, kep3::elements_type el_t) {
    {
        p.elements(e, el_t)
    } -> std::same_as<std::array<double, 6>>;
};

// Concept detecting if the type T can be used as a udpla.
template <typename T>
concept any_udpla = std::default_initializable<T> && std::copy_constructible<T> && std::move_constructible<T>
                    && std::destructible<T> && udpla_has_eph<T>;

template <typename, typename>
struct planet_iface {
};

// Planet interface.
template <>
// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)
struct planet_iface<void, void> {
    virtual ~planet_iface() = default;

    [[nodiscard]] virtual double get_mu_central_body() const = 0;
    [[nodiscard]] virtual double get_mu_self() const = 0;
    [[nodiscard]] virtual double get_radius() const = 0;
    [[nodiscard]] virtual double get_safe_radius() const = 0;
    [[nodiscard]] virtual std::string get_name() const = 0;
    [[nodiscard]] virtual std::string get_extra_info() const = 0;
    [[nodiscard]] virtual std::array<std::array<double, 3>, 2> eph(const epoch &) const = 0;
    // NOLINTNEXTLINE(google-default-arguments)
    [[nodiscard]] virtual double period(const epoch & = kep3::epoch()) const = 0;
    // NOLINTNEXTLINE(google-default-arguments)
    [[nodiscard]] virtual std::array<double, 6> elements(const kep3::epoch & = kep3::epoch(),
                                                         kep3::elements_type = kep3::elements_type::KEP_F) const
        = 0;
};

// Helper macro to implement getters in the planet interface implementation.
#define KEP3_UDPLA_IMPLEMENT_GET(name, type, def_value)                                                                \
    [[nodiscard]] type get_##name() const final                                                                        \
    {                                                                                                                  \
        if constexpr (udpla_has_get_##name<T>) {                                                                       \
            return this->value().get_##name();                                                                         \
        } else {                                                                                                       \
            return def_value;                                                                                          \
        }                                                                                                              \
    }

// NOTE: implement this in the cpp in order to avoid
// instantiating the same code over and over.
kep3_DLL_PUBLIC double period_from_energy(const std::array<double, 3> &, const std::array<double, 3> &, double);
kep3_DLL_PUBLIC std::array<double, 6> elements_from_posvel(const std::array<std::array<double, 3>, 2> &, double,
                                                           kep3::elements_type);

// Planet interface implementation.
template <typename Holder, typename T>
    requires any_udpla<T>
struct planet_iface<Holder, T> : planet_iface<void, void>, tanuki::iface_impl_helper<Holder, T, planet_iface> {
    KEP3_UDPLA_IMPLEMENT_GET(mu_central_body, double, -1)
    KEP3_UDPLA_IMPLEMENT_GET(mu_self, double, -1)
    KEP3_UDPLA_IMPLEMENT_GET(radius, double, -1)
    KEP3_UDPLA_IMPLEMENT_GET(safe_radius, double, -1)
    KEP3_UDPLA_IMPLEMENT_GET(name, std::string, boost::core::demangle(typeid(T).name()))
    KEP3_UDPLA_IMPLEMENT_GET(extra_info, std::string, "")

    [[nodiscard]] std::array<std::array<double, 3>, 2> eph(const epoch &ep) const final
    {
        return this->value().eph(ep);
    }

    // NOLINTNEXTLINE(google-default-arguments)
    [[nodiscard]] double period(const epoch &ep = kep3::epoch()) const final
    {
        // If the user provides an efficient way to compute the period, then use it
        if constexpr (udpla_has_period<T>) {
            return this->value().period(ep);
        } else if constexpr (udpla_has_get_mu_central_body<T>) {
            // If the user provides the central body parameter, then compute the
            // period from the energy at epoch
            auto [r, v] = eph(ep);
            double mu = get_mu_central_body();
            return period_from_energy(r, v, mu);
        } else {
            // There is no way to compute a period for this planet
            throw not_implemented_error(fmt::format("A period nor a central body mu has been declared for '{}', "
                                                    "impossible to provide a default implementation",
                                                    get_name()));
        }
    }

    // NOLINTNEXTLINE(google-default-arguments)
    [[nodiscard]] std::array<double, 6> elements(const kep3::epoch &ep = kep3::epoch(),
                                                 kep3::elements_type el_type = kep3::elements_type::KEP_F) const final
    {
        // If the user provides an efficient way to compute the orbital elements, then use it.
        if constexpr (udpla_has_elements<T>) {
            return this->value().elements(ep, el_type);
        } else if constexpr (udpla_has_get_mu_central_body<T>) {
            // If the user provides the central body parameter, then compute the
            // elements using posvel computed at ep and converted.
            auto pos_vel = eph(ep);
            double mu = get_mu_central_body();
            return elements_from_posvel(pos_vel, mu, el_type);
        } else {
            // There is no way to compute osculating elements for this planet
            throw not_implemented_error(
                fmt::format("A central body mu has not been declared for '{}', "
                            "impossible to provide a default implementation to compute the osculating elements",
                            get_name()));
        }
    }
};

#undef KEP3_UDPLA_IMPLEMENT_GET

// The udpla used in the default constructor of planet.
struct kep3_DLL_PUBLIC null_udpla {
    null_udpla() = default;
    static std::array<std::array<double, 3>, 2> eph(const epoch &);

private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &, unsigned){};
};

} // namespace kep3::detail

// Make a def-constructed planet serialisable.
TANUKI_S11N_WRAP_EXPORT_KEY(kep3::detail::null_udpla, kep3::detail::planet_iface)

namespace tanuki
{

// Implement the reference interface for the planet class.
template <typename Wrap>
struct ref_iface<Wrap, kep3::detail::planet_iface> {
    TANUKI_REF_IFACE_MEMFUN(get_mu_central_body)
    TANUKI_REF_IFACE_MEMFUN(get_mu_self)
    TANUKI_REF_IFACE_MEMFUN(get_radius)
    TANUKI_REF_IFACE_MEMFUN(get_safe_radius)
    TANUKI_REF_IFACE_MEMFUN(get_name)
    TANUKI_REF_IFACE_MEMFUN(get_extra_info)
    TANUKI_REF_IFACE_MEMFUN(eph)
    TANUKI_REF_IFACE_MEMFUN(period)
    TANUKI_REF_IFACE_MEMFUN(elements)

    // Implement the extract functionality.
    template <typename T>
    T *extract()
    {
        return value_ptr<T>(*static_cast<Wrap *>(this));
    }
};

} // namespace tanuki

namespace kep3
{

#if defined(__GNUC__)

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"

#endif

// Definition of the planet class.
using planet = tanuki::wrap<detail::planet_iface, tanuki::config<detail::null_udpla>{.pointer_interface = false}>;

#if defined(__GNUC__)

#pragma GCC diagnostic pop

#endif

namespace detail
{

// Streaming operator for planet.
kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &, const planet &);

} // namespace detail

} // namespace kep3

#endif // kep3_PLANET_H
