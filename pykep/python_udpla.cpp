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

#include <pybind11/pybind11.h>

#include <kep3/planet.hpp>

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
    return py::cast<std::array<std::array<double, 3>, 2>>(m_obj.attr("eph")(ep));
}

// Optional methods
[[nodiscard]] std::string python_udpla::get_name() const
{
    auto bf = pykep::callable_attribute(m_obj, "get_name");
    if (bf.is_none()) {
        pykep::py_throw(PyExc_NotImplementedError, ("the get_name() method has been invoked, but it is not implemented "
                                                    "in the user-defined Python planet '"
                                                    + pykep::str(m_obj) + "' of type '" + pykep::str(pykep::type(m_obj))
                                                    + "': the method is either not present or not callable")
                                                       .c_str());
    }
    return py::cast<std::string>(m_obj.attr("get_name")());
}
[[nodiscard]] std::string python_udpla::get_extra_info() const
{
    auto bf = pykep::callable_attribute(m_obj, "get_extra_info");
    if (bf.is_none()) {
        pykep::py_throw(PyExc_NotImplementedError,
                        ("the get_extra_info() method has been invoked, but it is not implemented "
                         "in the user-defined Python planet '"
                         + pykep::str(m_obj) + "' of type '" + pykep::str(pykep::type(m_obj))
                         + "': the method is either not present or not callable")
                            .c_str());
    }
    return py::cast<std::string>(m_obj.attr("get_extra_info")());
}
[[nodiscard]] double python_udpla::get_mu_central_body() const
{
    auto bf = pykep::callable_attribute(m_obj, "get_mu_central_body");
    if (bf.is_none()) {
        pykep::py_throw(PyExc_NotImplementedError,
                        ("the get_mu_central_body() method has been invoked, but it is not implemented "
                         "in the user-defined Python planet '"
                         + pykep::str(m_obj) + "' of type '" + pykep::str(pykep::type(m_obj))
                         + "': the method is either not present or not callable")
                            .c_str());
    }
    return py::cast<double>(m_obj.attr("get_mu_central_body")());
}
[[nodiscard]] double python_udpla::get_mu_self() const
{
    auto bf = pykep::callable_attribute(m_obj, "get_mu_self");
    if (bf.is_none()) {
        pykep::py_throw(PyExc_NotImplementedError,
                        ("the get_mu_self() method has been invoked, but it is not implemented "
                         "in the user-defined Python planet '"
                         + pykep::str(m_obj) + "' of type '" + pykep::str(pykep::type(m_obj))
                         + "': the method is either not present or not callable")
                            .c_str());
    }
    return py::cast<double>(m_obj.attr("get_mu_self")());
}
[[nodiscard]] double python_udpla::get_radius() const
{
    auto bf = pykep::callable_attribute(m_obj, "get_radius");
    if (bf.is_none()) {
        pykep::py_throw(PyExc_NotImplementedError,
                        ("the get_radius() method has been invoked, but it is not implemented "
                         "in the user-defined Python planet '"
                         + pykep::str(m_obj) + "' of type '" + pykep::str(pykep::type(m_obj))
                         + "': the method is either not present or not callable")
                            .c_str());
    }
    return py::cast<double>(m_obj.attr("get_radius")());
}
[[nodiscard]] double python_udpla::get_safe_radius() const
{
    auto bf = pykep::callable_attribute(m_obj, "get_safe_radius");
    if (bf.is_none()) {
        pykep::py_throw(PyExc_NotImplementedError,
                        ("the get_safe_radius() method has been invoked, but it is not implemented "
                         "in the user-defined Python planet '"
                         + pykep::str(m_obj) + "' of type '" + pykep::str(pykep::type(m_obj))
                         + "': the method is either not present or not callable")
                            .c_str());
    }
    return py::cast<double>(m_obj.attr("get_safe_radius")());
}
[[nodiscard]] double python_udpla::period() const
{
    auto bf = pykep::callable_attribute(m_obj, "period");
    if (bf.is_none()) {
        pykep::py_throw(PyExc_NotImplementedError, ("the period() method has been invoked, but it is not implemented "
                                                    "in the user-defined Python planet '"
                                                    + pykep::str(m_obj) + "' of type '" + pykep::str(pykep::type(m_obj))
                                                    + "': the method is either not present or not callable")
                                                       .c_str());
    }
    return py::cast<double>(m_obj.attr("period")());
}
} // namespace pykep
