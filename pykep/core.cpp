// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <boost/optional.hpp>
#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "docstrings.hpp"

namespace py = pybind11;
namespace pk = pykep;

PYBIND11_MODULE(core, m) {
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
      .value("KEP_M", kep3::KEP_M,
             "Keplerian Elements a,e,i,W,w,M (Mean anomaly)")
      .value("KEP_F", kep3::KEP_F,
             "Keplerian Elements a,e,i,W,w,f (True anomaly)")
      .value("MEQ", kep3::MEQ,
             "Modified Equinoctial Elements p,f,g,h,k,L (Mean Longitude)")
      .value("MEQ_R", kep3::MEQ_R,
             "Modified Equinoctial Elements (retrograde) p,f,g,h,k,L (Mean "
             "Longitude)")
      .value("POSVEL", kep3::POSVEL, "Position and Velocity")
      .export_values();

  // We expose the various anomaly conversions
  m.def("m2e", &kep3::m2e, pk::m2e_doc().c_str());
}