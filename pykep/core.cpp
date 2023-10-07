// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <kep3/epoch.hpp>
#include <string>

#include <boost/optional.hpp>
#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/planet.hpp>
#include <kep3/planets/keplerian.hpp>
#include <pybind11/numpy.h>
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
    py::class_<kep3::epoch>(m, "epoch").def(py::init<double>());

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
