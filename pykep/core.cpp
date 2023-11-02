// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include <pybind11/pytypes.h>
#include <string>

#include <fmt/chrono.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/core_astro/eq2par2eq.hpp>
#include <kep3/core_astro/ic2eq2ic.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/epoch.hpp>
#include <kep3/lambert_problem.hpp>
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
    py::enum_<kep3::elements_type>(m, "el_type", "")
        .value("KEP_M", kep3::KEP_M, "Keplerian Elements :math:`[a,e,i,\\Omega,\\omega,M]` (Mean anomaly)")
        .value("KEP_F", kep3::KEP_F, "Keplerian Elements :math:`[a,e,i,\\Omega,\\omega,f]` (True anomaly)")
        .value("MEQ", kep3::MEQ, "Modified Equinoctial Elements :math:`[p,f,g,h,k,L]` (Mean Longitude)")
        .value("MEQ_R", kep3::MEQ_R,
               "Modified Equinoctial Elements (retrograde) :math:`[p,f,g,h,k,L]` (Mean "
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

    // Eposing element conversions
    m.def("ic2par", &kep3::ic2par);
    m.def("par2ic", &kep3::par2ic);
    m.def("ic2eq", &kep3::ic2eq);
    m.def("eq2ic", &kep3::eq2ic);
    m.def("par2eq", &kep3::par2eq);
    m.def("eq2par", &kep3::eq2par);

    // Class epoch
    py::class_<kep3::epoch> epoch_class(m, "epoch", "Represents a specific point in time.");

    py::enum_<kep3::epoch::julian_type>(epoch_class, "julian_type")
        .value("MJD2000", kep3::epoch::julian_type::MJD2000, "Modified Julian Date 2000.")
        .value("MJD", kep3::epoch::julian_type::MJD, "Modified Julian Date.")
        .value("JD", kep3::epoch::julian_type::JD, "Julian Date.");

    py::enum_<kep3::epoch::string_format>(epoch_class, "string_format")
        .value("ISO", kep3::epoch::string_format::ISO, "ISO 8601 format for dates.");

    // This must go after the enum class registration
    epoch_class
        // Construtor from julian floats/int
        .def(py::init<double, kep3::epoch::julian_type>(), py::arg("when"),
             py::arg("julian_type") = kep3::epoch::julian_type::MJD2000, pk::epoch_from_float_doc().c_str())
        .def(py::init<int, kep3::epoch::julian_type>(), py::arg("when"),
             py::arg("julian_type") = kep3::epoch::julian_type::MJD2000)
        // Constructor from string
        .def(py::init<std::string, kep3::epoch::string_format>(), py::arg("when"),
             py::arg("string_format") = kep3::epoch::string_format::ISO, pk::epoch_from_string_doc().c_str())
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
             py::arg("when"), pk::epoch_from_datetime_doc().c_str())
        // repr()
        .def("__repr__", &pykep::ostream_repr<kep3::epoch>)
        // Copy and deepcopy.
        .def("__copy__", &pykep::generic_copy_wrapper<kep3::epoch>)
        .def("__deepcopy__", &pykep::generic_deepcopy_wrapper<kep3::epoch>)
        // Pickle support.
        .def(py::pickle(&pykep::pickle_getstate_wrapper<kep3::epoch>, &pykep::pickle_setstate_wrapper<kep3::epoch>))
        // now().
        .def_static("now", &kep3::epoch::now, "Returns a pykep.epoch with the current UTC date.")
        // julian dates
        .def_property_readonly("mjd2000", &kep3::epoch::mjd2000, "The Modified Julian Date 2000")
        .def_property_readonly("mjd", &kep3::epoch::mjd, "The Modified Julian Date")
        .def_property_readonly("jd", &kep3::epoch::jd, "The Julian Date")
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

    // Class planet (type erasure machinery here)
    py::class_<kep3::planet> planet_class(m, "planet", py::dynamic_attr{}, pykep::planet_docstring().c_str());
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
        "eph",
        [](const kep3::planet &pl, const std::variant<double, kep3::epoch> &when) {
            return std::visit([&](const auto &v) { return pl.eph(v); }, when);
        },
        py::arg("when"), pykep::planet_eph_docstring().c_str());
    // Vectorized version. Note that the udpla method flattens everything but planet returns a non flat array.
    planet_class.def(
        "eph_v",
        [](const kep3::planet &pl, const std::vector<double> &eps) {
            std::vector<double> res = pl.eph_v(eps);
            // We create a capsule for the py::array_t to manage ownership change.
            auto vec_ptr = std::make_unique<std::vector<double>>(std::move(res));

            py::capsule vec_caps(vec_ptr.get(), [](void *ptr) {
                std::unique_ptr<std::vector<double>> vptr(static_cast<std::vector<double> *>(ptr));
            });

            // NOTE: at this point, the capsule has been created successfully (including
            // the registration of the destructor). We can thus release ownership from vec_ptr,
            // as now the capsule is responsible for destroying its contents. If the capsule constructor
            // throws, the destructor function is not registered/invoked, and the destructor
            // of vec_ptr will take care of cleaning up.
            auto *ptr = vec_ptr.release();

            return py::array_t<double>(py::array::ShapeContainer{boost::numeric_cast<py::ssize_t>(eps.size()),
                                                                 static_cast<py::ssize_t>(6)}, // shape
                                       ptr->data(), std::move(vec_caps));
        },
        py::arg("when"), pykep::planet_eph_v_docstring().c_str());

#define PYKEP3_EXPOSE_PLANET_GETTER(name)                                                                              \
    planet_class.def(                                                                                                  \
        "get_" #name, [](const kep3::planet &pl) { return pl.get_##name(); },                                          \
        pykep::planet_get_##name##_docstring().c_str());

    PYKEP3_EXPOSE_PLANET_GETTER(name);
    PYKEP3_EXPOSE_PLANET_GETTER(extra_info);
    PYKEP3_EXPOSE_PLANET_GETTER(mu_central_body);
    PYKEP3_EXPOSE_PLANET_GETTER(mu_self);
    PYKEP3_EXPOSE_PLANET_GETTER(radius);
    PYKEP3_EXPOSE_PLANET_GETTER(safe_radius);

#undef PYKEP3_EXPOSE_PLANET_GETTER

    planet_class.def(
        "period",
        [](const kep3::planet &pl, const std::variant<double, kep3::epoch> &when) {
            return std::visit([&](const auto &v) { return pl.period(v); }, when);
        },
        py::arg("when") = 0., pykep::planet_period_docstring().c_str());

    planet_class.def(
        "elements",
        [](const kep3::planet &pl, const std::variant<double, kep3::epoch> &when, kep3::elements_type el_ty) {
            return std::visit([&](const auto &v) { return pl.elements(v, el_ty); }, when);
        },
        py::arg("when") = 0., py::arg("el_type") = kep3::elements_type::KEP_F,
        pykep::planet_elements_docstring().c_str());

    // We now expose the cpp udplas. They will also add a constructor and the extract machinery to the planet_class
    // UDPLA module
    pykep::expose_all_udplas(m, planet_class);

    // Finalize (this constructor must be the last one of planet_class: else overload will fail with all the others)
    planet_class.def(py::init([](const py::object &o) { return kep3::planet{pk::python_udpla(o)}; }), py::arg("udpla"));

    // Exposing the Lambert problem class
    py::class_<kep3::lambert_problem> lambert_problem(m, "lambert_problem", pykep::lambert_problem_docstring().c_str());
    lambert_problem
        .def(py::init<const std::array<double, 3> &, const std::array<double, 3> &, double, double, bool, unsigned>(),
             py::arg("r0") = std::array<double, 3>{{1., 0., 0}}, py::arg("r1") = std::array<double, 3>{{0., 1., 0}},
             py::arg("tof") = kep3::pi / 2, py::arg("mu") = 1., py::arg("cw") = false, py::arg("multi_revs") = 1)
        // repr().
        .def("__repr__", &pykep::ostream_repr<kep3::lambert_problem>)
        // Copy and deepcopy.
        .def("__copy__", &pykep::generic_copy_wrapper<kep3::lambert_problem>)
        .def("__deepcopy__", &pykep::generic_deepcopy_wrapper<kep3::lambert_problem>)
        // Pickle support.
        .def(py::pickle(&pykep::pickle_getstate_wrapper<kep3::lambert_problem>,
                        &pykep::pickle_setstate_wrapper<kep3::lambert_problem>))
        .def_property_readonly("v0", &kep3::lambert_problem::get_v0, "The velocity at the first point.")
        .def_property_readonly("v1", &kep3::lambert_problem::get_v1, "The velocity at the second point.")
        .def_property_readonly("r0", &kep3::lambert_problem::get_r0, "The first point.")
        .def_property_readonly("r1", &kep3::lambert_problem::get_r1, "The second point.")
        .def_property_readonly("tof", &kep3::lambert_problem::get_tof, "The time of flight between the two points.")
        .def_property_readonly("mu", &kep3::lambert_problem::get_mu,
                               "The gravitational parameter of the attracting body.")
        .def_property_readonly("x", &kep3::lambert_problem::get_x,
                               "The Battin variable x along the time of flight curves.")
        .def_property_readonly("iters", &kep3::lambert_problem::get_iters, "The number of iterations made.")
        .def_property_readonly("Nmax", &kep3::lambert_problem::get_Nmax, "The maximum number of iterations allowed.");

    // Exposing propagators
    m.def(
        "propagate_lagrangian",
        [](const std::array<std::array<double, 3>, 2> &pos_vel, double dt, double mu, bool request_stm) {
            auto pl_retval = kep3::propagate_lagrangian(pos_vel, dt, mu, request_stm);
            if (pl_retval.second) {
                // The stm was requested lets transfer ownership to python
                std::array<double, 36> &stm = pl_retval.second.value();

                // We create a capsule for the py::array_t to manage ownership change.
                auto vec_ptr = std::make_unique<std::array<double, 36>>(stm);

                py::capsule vec_caps(vec_ptr.get(), [](void *ptr) {
                    std::unique_ptr<std::array<double, 36>> vptr(static_cast<std::array<double, 36> *>(ptr));
                });

                // NOTE: at this point, the capsule has been created successfully (including
                // the registration of the destructor). We can thus release ownership from vec_ptr,
                // as now the capsule is responsible for destroying its contents. If the capsule constructor
                // throws, the destructor function is not registered/invoked, and the destructor
                // of vec_ptr will take care of cleaning up.
                auto *ptr = vec_ptr.release();

                auto computed_stm = py::array_t<double>(
                    py::array::ShapeContainer{static_cast<py::ssize_t>(6), static_cast<py::ssize_t>(6)}, // shape
                    ptr->data(), std::move(vec_caps));
                return py::make_tuple(pl_retval.first, computed_stm);
            } else {
                return py::make_tuple(pl_retval.first[0], pl_retval.first[1]);
            }
        },
        py::arg("rv") = std::array<std::array<double, 3>, 2>{{{1, 0, 0}, {0, 1, 0}}}, py::arg("tof") = kep3::pi / 2,
        py::arg("mu") = 1, py::arg("stm") = false, pykep::propagate_lagrangian_docstring().c_str());
    m.def("propagate_lagrangian_v", &kep3::propagate_lagrangian_v, py::arg("rv") = std::array<std::array<double, 3>, 2>{{{1, 0, 0}, {0, 1, 0}}}, py::arg("tofs") = std::vector<double>{kep3::pi / 2,},
        py::arg("mu") = 1, py::arg("stm") = false, pykep::propagate_lagrangian_v_docstring().c_str());
}
