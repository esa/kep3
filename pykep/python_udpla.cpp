// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <memory>
#include <string>
#include <type_traits>
#include <typeindex>

#include <fmt/core.h>
#include <fmt/std.h>
#include <kep3/planet.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "common_utils.hpp"
#include "python_udpla.hpp"

namespace pykep
{

python_udpla::python_udpla() = default;

python_udpla::python_udpla(py::object obj) : m_obj(std::move(obj))
{
    // Forbid the use of a pykep.planet as a UDPLA.
    if (pykep::type(m_obj).equal(py::module::import("pykep").attr("planet"))) {
        pykep::py_throw(PyExc_TypeError,
                        ("a pykep.planet cannot be used as a UDPLA for another pykep.problem (if you need to copy a "
                         "planet please use the standard Python copy()/deepcopy() functions)"));
    }
    //// Check that o is an instance of a class, and not just a type.
    check_not_type(m_obj, "planet");
    //// Check the presence of the mandatory methods
    check_mandatory_method(m_obj, "eph", "planet");
};

// Mandatory methods
[[nodiscard]] std::array<std::array<double, 3>, 2> python_udpla::eph(const kep3::epoch &ep) const
{
    auto bf = pykep::callable_attribute(m_obj, "eph");
    if (bf.is_none()) {
        pykep::py_throw(PyExc_NotImplementedError, ("the eph() method has been invoked, but it is not implemented "
                                                    "in the user-defined Python planet '"
                                                    + pykep::str(m_obj) + "' of type '" + pykep::str(pykep::type(m_obj))
                                                    + "': the method is either not present or not callable")
                                                       .c_str());
    }
    return py::cast<std::array<std::array<double, 3>, 2>>(m_obj.attr("eph")(ep));
}

// Optional methods
[[nodiscard]] std::string python_udpla::get_name() const
{
    return getter_wrapper<std::string>(m_obj, "get_name", pykep::str(pykep::type(m_obj)));
}
[[nodiscard]] std::string python_udpla::get_extra_info() const
{
    return getter_wrapper<std::string>(m_obj, "get_extra_info", "");
}
[[nodiscard]] double python_udpla::get_mu_central_body() const
{
    return getter_wrapper<double>(m_obj, "get_mu_central_body", -1);
}
[[nodiscard]] double python_udpla::get_mu_self() const
{
    return getter_wrapper<double>(m_obj, "get_mu_self", -1);
}
[[nodiscard]] double python_udpla::get_radius() const
{
    return getter_wrapper<double>(m_obj, "get_radius", -1);
}
[[nodiscard]] double python_udpla::get_safe_radius() const
{
    return getter_wrapper<double>(m_obj, "get_safe_radius", -1);
}
[[nodiscard]] double python_udpla::period(const kep3::epoch &ep) const
{
    auto method_period = pykep::callable_attribute(m_obj, "period");
    auto method_mu = pykep::callable_attribute(m_obj, "get_mu_central_body");

    if (method_period.is_none()) {
        if (method_mu.is_none()) {
            pykep::py_throw(
                PyExc_NotImplementedError,
                ("the period() method has been invoked, but "
                 "in the user-defined Python planet '"
                 + pykep::str(m_obj) + "' of type '" + pykep::str(pykep::type(m_obj))
                 + "': the methods period() or get_mu_central_body() are either not present or not callable")
                    .c_str());
        } else {
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
        }
    }
    return py::cast<double>(m_obj.attr("period")());
}

} // namespace pykep
