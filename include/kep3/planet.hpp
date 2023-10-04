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

#include <cmath>
#include <concepts>
#include <limits>
#include <string>
#include <typeinfo>

#include <boost/core/demangle.hpp>

#include <fmt/core.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/detail/s11n.hpp>
#include <kep3/detail/type_traits.hpp>
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

template <typename>
struct planet_iface;

// Planet interface.
template <>
// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)
struct planet_iface<void> {
    virtual ~planet_iface() = default;

    [[nodiscard]] virtual double get_mu_central_body() const = 0;
    [[nodiscard]] virtual double get_mu_self() const = 0;
    [[nodiscard]] virtual double get_radius() const = 0;
    [[nodiscard]] virtual double get_safe_radius() const = 0;
    [[nodiscard]] virtual std::string get_name() const = 0;
    [[nodiscard]] virtual std::string get_extra_info() const = 0;
    [[nodiscard]] virtual std::array<std::array<double, 3>, 2> eph(const epoch &) const = 0;
    [[nodiscard]] virtual double period(const epoch &) const = 0;
};

#define KEP3_UDPLA_IMPLEMENT_GET(name, type, def_value)                                                                \
    [[nodiscard]] type get_##name() const final                                                                        \
    {                                                                                                                  \
        if constexpr (udpla_has_get_##name<typename Holder::value_type>) {                                             \
            return this->value().get_##name();                                                                         \
        } else {                                                                                                       \
            return def_value;                                                                                          \
        }                                                                                                              \
    }

// Planet interface implementation.
template <typename Holder>
struct planet_iface : planet_iface<void>, tanuki::iface_impl_helper<Holder> {
    KEP3_UDPLA_IMPLEMENT_GET(mu_central_body, double, -1)
    KEP3_UDPLA_IMPLEMENT_GET(mu_self, double, -1)
    KEP3_UDPLA_IMPLEMENT_GET(radius, double, -1)
    KEP3_UDPLA_IMPLEMENT_GET(safe_radius, double, -1)
    KEP3_UDPLA_IMPLEMENT_GET(name, std::string, boost::core::demangle(typeid(typename Holder::value_type).name()))
    KEP3_UDPLA_IMPLEMENT_GET(extra_info, std::string, "")

    [[nodiscard]] std::array<std::array<double, 3>, 2> eph(const epoch &ep) const final
    {
        return this->value().eph(ep);
    }

    [[nodiscard]] double period(const epoch &ep) const final
    {
        using value_type = typename Holder::value_type;

        // If the user provides an efficient way to compute the period, then use it
        if constexpr (udpla_has_period<value_type>) {
            return this->value().period(ep);
        } else if constexpr (udpla_has_get_mu_central_body<value_type>) {
            // If the user provides the central body parameter, then compute the
            // period from the energy at epoch
            auto [r, v] = eph(ep);
            double mu = get_mu_central_body();
            double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
            double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            double en = v2 / 2. - mu / R;
            if (en > 0) {
                // If the energy is positive we have an hyperbolae and we return nan
                return std::numeric_limits<double>::quiet_NaN();
            } else {
                double a = -mu / 2. / en;
                return kep3::pi * 2. * std::sqrt(a * a * a / mu);
            }
        } else {
            // There is no way to compute a period for this planet
            throw not_implemented_error(fmt::format("A period nor a central body mu has been declared for '{}', "
                                                    "impossible to provide a default implementation",
                                                    get_name()));
        }
    }
};

#undef KEP3_UDPLA_IMPLEMENT_GET

template <typename T>
concept any_udpla = udpla_has_eph<T>;

struct kep3_DLL_PUBLIC null_udpla {
    null_udpla() = default;
    static std::array<std::array<double, 3>, 2> eph(const epoch &);

private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &, unsigned){};
};
} // namespace kep3::detail

TANUKI_S11N_WRAP_EXPORT_KEY(kep3::detail::null_udpla, kep3::detail::planet_iface)

namespace tanuki
{

template <typename T>
inline constexpr bool is_wrappable<T, kep3::detail::planet_iface> = kep3::detail::any_udpla<T>;

template <typename Wrap>
struct ref_iface<Wrap, kep3::detail::planet_iface> {
    TANUKI_REF_IFACE_MEMFUN(get_mu_central_body)
    TANUKI_REF_IFACE_MEMFUN(get_mu_self)
    TANUKI_REF_IFACE_MEMFUN(get_radius)
    TANUKI_REF_IFACE_MEMFUN(get_safe_radius)
    TANUKI_REF_IFACE_MEMFUN(get_name)
    TANUKI_REF_IFACE_MEMFUN(get_extra_info)
    TANUKI_REF_IFACE_MEMFUN(eph)

    [[nodiscard]] double period(const kep3::epoch &ep = kep3::epoch()) const
    {
        return iface_ptr(*static_cast<const Wrap *>(this))->period(ep);
    }
};

} // namespace tanuki

namespace kep3
{

#if defined(__GNUC__)

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"

#endif

using planet = tanuki::wrap<detail::planet_iface, tanuki::config<detail::null_udpla>{.pointer_interface = false}>;

#if defined(__GNUC__)

#pragma GCC diagnostic pop

#endif

namespace detail
{

// Streaming operator for algorithm.
kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &, const planet &);

} // namespace detail

} // namespace kep3

#endif // kep3_PLANET_H
