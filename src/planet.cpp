// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <fmt/core.h>

#include <kep3/exceptions.hpp>
#include <kep3/planet.hpp>

namespace kep3::detail
{

std::array<std::array<double, 3>, 2> null_udpla::eph(const epoch &)
{
    std::array<double, 3> pos = {1., 0., 0.};
    std::array<double, 3> vel = {0., 1., 0.};
    return {pos, vel};
};
} // namespace kep3::detail

namespace kep3
{
planet::planet() : planet(detail::null_udpla{}){};
planet::planet(const planet &other) : m_ptr(other.m_ptr->clone()){};
planet::planet(planet &&other) noexcept : m_ptr(std::move(other.m_ptr)){};

planet &planet::operator=(planet &&other) noexcept
{
    if (this != &other) {
        m_ptr = std::move(other.m_ptr);
    }
    return *this;
}

planet &planet::operator=(const planet &other)
{
    return *this = planet(other);
}

bool planet::is_valid() const
{
    return static_cast<bool>(m_ptr);
}

std::type_index planet::get_type_index() const
{
    return ptr()->get_type_index();
}

const void *planet::get_ptr() const
{
    return ptr()->get_ptr();
}

void *planet::get_ptr()
{
    return ptr()->get_ptr();
}

detail::planet_inner_base const *planet::ptr() const
{
    assert(m_ptr.get() != nullptr);
    return m_ptr.get();
}

detail::planet_inner_base *planet::ptr()
{
    assert(m_ptr.get() != nullptr);
    return m_ptr.get();
}

std::array<std::array<double, 3>, 2> planet::eph(const epoch &ep)
{
    return ptr()->eph(ep);
}

double planet::period(const epoch &ep) const
{
    return ptr()->period(ep);
}

std::string planet::get_name() const
{
    return ptr()->get_name();
}

std::string planet::get_extra_info() const
{
    return ptr()->get_extra_info();
}

double planet::get_mu_central_body() const
{
    return ptr()->get_mu_central_body();
}

double planet::get_mu_self() const
{
    return ptr()->get_mu_self();
    ;
}

double planet::get_radius() const
{
    return ptr()->get_radius();
}

double planet::get_safe_radius() const
{
    return ptr()->get_safe_radius();
    ;
}

std::ostream &operator<<(std::ostream &os, const planet &p)
{
    os << "Planet name: " << p.get_name();
    os << "\nC++ class name: " << detail::demangle_from_typeid(p.get_type_index().name()) << '\n';
    os << fmt::format("mu central body (-1 if not defined): {}\n", p.get_mu_central_body());
    os << fmt::format("mu body (-1 if not defined): {}\n", p.get_mu_self());
    os << fmt::format("radius body (-1 if not defined): {}\n", p.get_radius());
    os << fmt::format("safe body radius (-1 if not defined): {}\n", p.get_safe_radius());

    const auto extra_str = p.get_extra_info();
    if (!extra_str.empty()) {
        os << "\nExtra info:\n" << extra_str;
    }
    return os;
}
} // namespace kep3