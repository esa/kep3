// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef PYKEP_PLANET_HPP
#define PYKEP_PLANET_HPP

#include <memory>
#include <string>
#include <type_traits>
#include <typeindex>

#include <pybind11/pybind11.h>

#include <kep3/planet.hpp>

namespace pykep
{
namespace py = pybind11;
struct python_udpla {
    py::object m_obj;

    python_udpla() = default;
    explicit python_udpla(py::object obj) : m_obj(std::move(obj)) {};

    [[nodiscard]] std::array<std::array<double, 3>, 2> eph(const kep3::epoch &) const
    {
        return py::cast<std::array<std::array<double, 3>, 2>>(m_obj.attr("eph")());
    }

    [[nodiscard]] std::string get_name() const
    {
        return py::cast<std::string>(m_obj.attr("get_name")());
    }


};
} // namespace pykep

#endif