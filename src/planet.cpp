// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/core/demangle.hpp>

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

namespace kep3::detail
{

std::ostream &operator<<(std::ostream &os, const planet &p)
{
    os << "Planet name: " << p.get_name();
    os << "\nC++ class name: " << boost::core::demangle(value_type_index(p).name()) << '\n';
    const auto extra_str = p.get_extra_info();
    if (!extra_str.empty()) {
        os << "\nExtra info:\n" << extra_str;
    }
    return os;
}

} // namespace kep3::detail

// NOLINTNEXTLINE
TANUKI_S11N_WRAP_EXPORT_IMPLEMENT(kep3::detail::null_udpla, kep3::detail::planet_iface)
