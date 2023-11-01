// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <cmath>
#include <limits>

#include <boost/core/demangle.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/core_astro/eq2par2eq.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/exceptions.hpp>
#include <kep3/planet.hpp>

namespace kep3::detail
{

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
double period_from_energy(const std::array<double, 3> &r, const std::array<double, 3> &v, double mu)
{
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
}

std::array<double, 6> elements_from_posvel(const std::array<std::array<double, 3>, 2> &pos_vel, double mu,
                                           kep3::elements_type el_type)
{
    std::array<double, 6> retval = kep3::ic2par(pos_vel, mu);
    switch (el_type) {
        case kep3::elements_type::KEP_F:
            break;
        case kep3::elements_type::KEP_M:
            if (retval[0] < 1) {
                throw std::logic_error("Mean anomaly is only available for ellipses.");
            }
            retval[5] = kep3::f2m(retval[5], retval[1]);
            break;
        case kep3::elements_type::MEQ:
            retval = kep3::par2eq(retval, false);
            break;
        case kep3::elements_type::MEQ_R:
            retval = kep3::par2eq(retval, true);
            break;
        // LCOV_EXCL_START
        default:
            throw std::logic_error("You should not go here!");
            // LCOV_EXCL_END
    }
    return retval;
}

std::array<std::array<double, 3>, 2> null_udpla::eph(double)
{
    std::array<double, 3> pos = {1., 0., 0.};
    std::array<double, 3> vel = {0., 1., 0.};
    return {pos, vel};
}

} // namespace kep3::detail

namespace kep3::detail
{

std::ostream &operator<<(std::ostream &os, const planet &p)
{
    os << "Planet name: " << p.get_name();
    os << "\nC++ class name: " << boost::core::demangle(value_type_index(p).name()) << '\n';
    os << fmt::format("\nmu central body (-1 if not defined): {}\n", p.get_mu_central_body());
    os << fmt::format("mu body (-1 if not defined): {}\n", p.get_mu_self());
    os << fmt::format("radius body (-1 if not defined): {}\n", p.get_radius());
    os << fmt::format("safe body radius (-1 if not defined): {}\n", p.get_safe_radius());

    const auto extra_str = p.get_extra_info();
    if (!extra_str.empty()) {
        os << "\nExtra info:\n" << extra_str;
    }
    return os;
}

} // namespace kep3::detail

// NOLINTNEXTLINE
TANUKI_S11N_WRAP_EXPORT_IMPLEMENT(kep3::detail::null_udpla, kep3::detail::planet_iface)
