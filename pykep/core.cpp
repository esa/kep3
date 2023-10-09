// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include <string>

#include <fmt/chrono.h>
#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/epoch.hpp>
#include <kep3/planet.hpp>
#include <kep3/planets/keplerian.hpp>
#include <pybind11/chrono.h>
#include <pybind11/detail/common.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "common_utils.hpp"
#include "docstrings.hpp"
#include "expose_udplas.hpp"
#include "python_udpla.hpp"

namespace py = pybind11;
namespace pk = pykep;

PYBIND11_MODULE(core, m)
{
    py::options options;
    options.disable_function_signatures();
    m.doc() = pk::core_module_doc();

    // We expose various global constants:
    m.attr("AU") = py::float_(kep3::AU);
    m.attr("CAVENDISH") = py::float_(kep3::CAVENDISH);
    m.attr("MU_SUN") = py::float_(kep3::MU_SUN);
    m.attr("MU_EARTH") = py::float_(kep3::MU_EARTH);
    m.attr("EARTH_VELOCITY") = py::float_(kep3::EARTH_VELOCITY);
    m.attr("EARTH_J2") = py::float_(kep3::EARTH_J2);
    m.attr("EARTH_RADIUS") = py::float_(kep3::EARTH_RADIUS);
    m.attr("RAD2DEG") = py::float_(kep3::RAD2DEG);
    m.attr("DAY2SEC") = py::float_(kep3::DAY2SEC);
    m.attr("SEC2DAY") = py::float_(kep3::SEC2DAY);
    m.attr("DAY2YEAR") = py::float_(kep3::DAY2YEAR);
    m.attr("G0") = py::float_(kep3::G0);

    // We expose here global enums:
    py::enum_<kep3::elements_type>(m, "elements_type", "")
        .value("KEP_M", kep3::KEP_M, "Keplerian Elements a,e,i,W,w,M (Mean anomaly)")
        .value("KEP_F", kep3::KEP_F, "Keplerian Elements a,e,i,W,w,f (True anomaly)")
        .value("MEQ", kep3::MEQ, "Modified Equinoctial Elements p,f,g,h,k,L (Mean Longitude)")
        .value("MEQ_R", kep3::MEQ_R,
               "Modified Equinoctial Elements (retrograde) p,f,g,h,k,L (Mean "
               "Longitude)")
        .value("POSVEL", kep3::POSVEL, "Position and Velocity")
        .export_values();

    // We expose the various anomaly conversions
    m.def("m2e", &kep3::m2e, pk::m2e_doc().c_str());
    m.def("e2m", &kep3::e2m, pk::e2m_doc().c_str());
    m.def("m2f", &kep3::m2f, pk::m2f_doc().c_str());
    m.def("f2m", &kep3::f2m, pk::f2m_doc().c_str());
    m.def("e2f", &kep3::e2f, pk::e2f_doc().c_str());
    m.def("f2e", &kep3::f2e, pk::f2e_doc().c_str());
    m.def("n2h", &kep3::n2h, pk::n2h_doc().c_str());
    m.def("h2n", &kep3::h2n, pk::h2n_doc().c_str());
    m.def("n2f", &kep3::n2f, pk::n2f_doc().c_str());
    m.def("f2n", &kep3::f2n, pk::f2n_doc().c_str());
    m.def("h2f", &kep3::h2f, pk::h2f_doc().c_str());
    m.def("f2h", &kep3::f2h, pk::f2h_doc().c_str());
    m.def("zeta2f", &kep3::zeta2f, pk::zeta2f_doc().c_str());
    m.def("f2zeta", &kep3::f2zeta, pk::f2zeta_doc().c_str());

    // And their vectorized versions
    m.def("m2e_v", py::vectorize(kep3::m2e), pk::m2e_v_doc().c_str());
    m.def("e2m_v", py::vectorize(kep3::e2m), pk::e2m_v_doc().c_str());
    m.def("m2f_v", py::vectorize(kep3::m2f), pk::m2f_v_doc().c_str());
    m.def("f2m_v", py::vectorize(kep3::f2m), pk::f2m_v_doc().c_str());
    m.def("e2f_v", py::vectorize(kep3::e2f), pk::e2f_v_doc().c_str());
    m.def("f2e_v", py::vectorize(kep3::f2e), pk::f2e_v_doc().c_str());
    m.def("n2h_v", py::vectorize(kep3::n2h), pk::n2h_v_doc().c_str());
    m.def("h2n_v", py::vectorize(kep3::h2n), pk::h2n_v_doc().c_str());
    m.def("n2f_v", py::vectorize(kep3::n2f), pk::n2f_v_doc().c_str());
    m.def("f2n_v", py::vectorize(kep3::f2n), pk::f2n_v_doc().c_str());
    m.def("h2f_v", py::vectorize(kep3::h2f), pk::h2f_v_doc().c_str());
    m.def("f2h_v", py::vectorize(kep3::f2h), pk::f2h_v_doc().c_str());
    m.def("zeta2f_v", py::vectorize(kep3::zeta2f), pk::zeta2f_v_doc().c_str());
    m.def("f2zeta_v", py::vectorize(kep3::f2zeta), pk::f2zeta_v_doc().c_str());

    // Class epoch
    py::class_<kep3::epoch> epoch_class(m, "epoch");

    py::enum_<kep3::epoch::julian_type>(epoch_class, "julian_type")
        .value("MJD2000", kep3::epoch::julian_type::MJD2000, "Modified Julian Date 2000")
        .value("MJD", kep3::epoch::julian_type::MJD, "Modified Julian Date")
        .value("JD", kep3::epoch::julian_type::JD, "Julian Date");

    py::enum_<kep3::epoch::string_format>(epoch_class, "string_format")
        .value("ISO", kep3::epoch::string_format::ISO, "ISO 8601 format for dates");

    // This must go after the enum class registration
    epoch_class
        // Construtor from julian floats/int
        .def(py::init<double, kep3::epoch::julian_type>(), py::arg("when"),
             py::arg("julian_type") = kep3::epoch::julian_type::MJD2000)
        .def(py::init<int, kep3::epoch::julian_type>(), py::arg("when"),
             py::arg("julian_type") = kep3::epoch::julian_type::MJD2000)
        // Constructor from string
        .def(py::init<std::string, kep3::epoch::string_format>(), py::arg("when"),
             py::arg("string_format") = kep3::epoch::string_format::ISO)
        // Constructor from datetime py::object
        .def(py::init([](const py::object &in) {
                 // We check that `in` is a datetimeobject
                 py::object Datetime = py::module_::import("datetime").attr("datetime");
                 if (!py::isinstance(in, Datetime)) {
                     pykep::py_throw(PyExc_TypeError, ("it seems you are trying to construct kep3::epoch object from a "
                                                       "python object that is not of type datetime"));
                 }
                 // We collect its info
                 int y = in.attr("year").cast<int>();
                 auto m = in.attr("month").cast<unsigned>();
                 auto d = in.attr("day").cast<unsigned>();
                 int h = in.attr("hour").cast<int>();
                 int min = in.attr("minute").cast<int>();
                 int s = in.attr("second").cast<int>();
                 int us = in.attr("microsecond").cast<int>();
                 return kep3::epoch(y, m, d, h, min, s, 0, us);
             }),
             py::arg("when"))
        // repr()
        .def("__repr__", &pykep::ostream_repr<kep3::epoch>)
        // Copy and deepcopy.
        .def("__copy__", &pykep::generic_copy_wrapper<kep3::epoch>)
        .def("__deepcopy__", &pykep::generic_deepcopy_wrapper<kep3::epoch>)
        // Pickle support.
        .def(py::pickle(&pykep::pickle_getstate_wrapper<kep3::epoch>, &pykep::pickle_setstate_wrapper<kep3::epoch>))
        // julian dates
        .def("mjd2000", &kep3::epoch::mjd2000)
        .def("mjd", &kep3::epoch::mjd)
        .def("jd", &kep3::epoch::jd)
        // comparison operators
        .def("__lt__", [](const kep3::epoch &ep1, const kep3::epoch &ep2) { return ep1 < ep2; })
        .def("__gt__", [](const kep3::epoch &ep1, const kep3::epoch &ep2) { return ep1 > ep2; })
        .def("__le__", [](const kep3::epoch &ep1, const kep3::epoch &ep2) { return ep1 <= ep2; })
        .def("__ge__", [](const kep3::epoch &ep1, const kep3::epoch &ep2) { return ep1 >= ep2; })
        .def("__eq__", [](const kep3::epoch &ep1, const kep3::epoch &ep2) { return ep1 == ep2; })
        .def("__ne__", [](const kep3::epoch &ep1, const kep3::epoch &ep2) { return ep1 != ep2; })
        // math
        .def("__add__",
             [](kep3::epoch ep, double dt) { return ep + std::chrono::duration<double, std::ratio<86400>>(dt); })
        .def("__add__", [](kep3::epoch ep, std::chrono::duration<double, std::ratio<1>> dt) { return ep + dt; })
        .def("__sub__",
             [](kep3::epoch ep, double dt) { return ep - std::chrono::duration<double, std::ratio<86400>>(dt); })
        .def("__sub__", [](kep3::epoch ep, std::chrono::duration<double, std::ratio<1>> dt) { return ep - dt; });

    // Epoch related utils
    m.def("utc_now", &kep3::utc_now);

    // Class planet (type erasure machinery here)
    py::class_<kep3::planet> planet_class(m, "planet", py::dynamic_attr{});
    // Constructor.
    // Expose extract.
    planet_class.def("_cpp_extract", &pykep::generic_cpp_extract<kep3::planet, kep3::udpla::keplerian>,
                     py::return_value_policy::reference_internal);
    // repr().
    planet_class.def("__repr__", &pykep::ostream_repr<kep3::planet>);
    // Copy and deepcopy.
    planet_class.def("__copy__", &pykep::generic_copy_wrapper<kep3::planet>);
    planet_class.def("__deepcopy__", &pykep::generic_deepcopy_wrapper<kep3::planet>);
    // UDPLA extraction for python stuff.
    planet_class.def("_py_extract", &pykep::generic_py_extract<kep3::planet>);
    // Pickle support.
    planet_class.def(
        py::pickle(&pykep::pickle_getstate_wrapper<kep3::planet>, &pykep::pickle_setstate_wrapper<kep3::planet>));
    // Planet methods.
    planet_class.def(
        "eph", [](const kep3::planet &pl, const kep3::epoch &ep) { return pl.eph(ep); }, py::arg("ep"));

#define PYKEP3_EXPOSE_PLANET_GETTER(name)                                                                              \
    planet_class.def("get_" #name, [](const kep3::planet &pl) { return pl.get_##name(); })

    PYKEP3_EXPOSE_PLANET_GETTER(name);
    PYKEP3_EXPOSE_PLANET_GETTER(extra_info);
    PYKEP3_EXPOSE_PLANET_GETTER(mu_central_body);
    PYKEP3_EXPOSE_PLANET_GETTER(mu_self);
    PYKEP3_EXPOSE_PLANET_GETTER(radius);
    PYKEP3_EXPOSE_PLANET_GETTER(safe_radius);

#undef PYKEP3_EXPOSE_PLANET_GETTER

    planet_class.def(
        "period", [](const kep3::planet &pl, const kep3::epoch &ep) { return pl.period(ep); },
        py::arg("ep") = kep3::epoch{});

    // We now expose the cpp udplas. They will also add a constructor and the extract machinery to the planet_class
    // UDPLA module
    auto udpla_module = m.def_submodule("udpla", "User defined planets that can construct a pykep.planet");
    pykep::expose_all_udplas(udpla_module, planet_class);

    // Finalize (this constructor must be the last one else overload will fail with all the others)
    planet_class.def(py::init([](const py::object &o) { return kep3::planet{pk::python_udpla(o)}; }), py::arg("udpla"));
}
