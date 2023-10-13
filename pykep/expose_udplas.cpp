// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <string>

#include <pybind11/pybind11.h>

#include <kep3/planet.hpp>
#include <kep3/planets/jpl_lp.hpp>
#include <kep3/planets/keplerian.hpp>

#include "common_utils.hpp"
#include "docstrings.hpp"
#include "expose_udplas.hpp"

namespace pykep
{
namespace py = pybind11;

// Split algorithm exposition functions.
void expose_all_udplas(py::module &udpla_module, py::class_<kep3::planet> &planet_class) // NOLINT
{
    // null udpla
    auto null_udpla = pykep::expose_one_udpla<kep3::detail::null_udpla>(
        udpla_module, planet_class, "_null_udpla", "A moot udpla used as default to construct a planet.");
    // Constructor.
    null_udpla.def(py::init<>());

    // keplerian udpla
    auto keplerian_udpla
        = pykep::expose_one_udpla<kep3::udpla::keplerian>(udpla_module, planet_class, "_keplerian", "keplerian udpla");
    // Constructors.
    keplerian_udpla
        .def(py::init<const kep3::epoch &, const std::array<double, 6> &, double, std::string, std::array<double, 3>,
                      kep3::elements_type>(),
             py::arg("ep"), py::arg("elem"), py::arg("mu_central_body"), py::arg("name") = "unknown",
             py::arg("added_params") = std::array<double, 3>({-1, -1, -1}),
             py::arg("elem_type") = kep3::elements_type::KEP_F)
        .def(py::init<const kep3::epoch &, const std::array<std::array<double, 3>, 2> &, double, std::string,
                      std::array<double, 3>>(),
             py::arg("ep"), py::arg("posvel"), py::arg("mu_central_body"), py::arg("name") = "unknown",
             py::arg("added_params") = std::array<double, 3>({-1, -1, -1}))
        // repr().
        .def("__repr__", &pykep::ostream_repr<kep3::udpla::keplerian>)
        // other methods
        .def_property_readonly("ref_epoch", &kep3::udpla::keplerian::get_ref_epoch);

    // jpl_lp_udpla udpla
    auto jpl_lp_udpla = pykep::expose_one_udpla<kep3::udpla::jpl_lp>(
        udpla_module, planet_class, "_jpl_lp", "Low precision ephemerides from JPL for the solar system planets");
    // Constructors.
    jpl_lp_udpla
        .def(py::init<std::string>(), py::arg("name") = "earth")
        // repr().
        .def("__repr__", &pykep::ostream_repr<kep3::udpla::jpl_lp>);
}

} // namespace pykep
