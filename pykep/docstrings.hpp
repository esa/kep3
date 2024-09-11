// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef PYKEP_DOCSTRINGS_HPP
#define PYKEP_DOCSTRINGS_HPP

#include "kep3/lambert_problem.hpp"
#include <string>

namespace pykep
{
// Modules
std::string core_module_doc();

// Anomaly conversions
std::string m2e_doc();
std::string e2m_doc();
std::string m2f_doc();
std::string f2m_doc();
std::string e2f_doc();
std::string f2e_doc();

std::string n2h_doc();
std::string h2n_doc();
std::string n2f_doc();
std::string f2n_doc();
std::string h2f_doc();
std::string f2h_doc();

std::string zeta2f_doc();
std::string f2zeta_doc();

// Anomaly conversions (vectorized)
std::string m2e_v_doc();
std::string e2m_v_doc();
std::string m2f_v_doc();
std::string f2m_v_doc();
std::string e2f_v_doc();
std::string f2e_v_doc();

std::string n2h_v_doc();
std::string h2n_v_doc();
std::string n2f_v_doc();
std::string f2n_v_doc();
std::string h2f_v_doc();
std::string f2h_v_doc();

std::string zeta2f_v_doc();
std::string f2zeta_v_doc();

// Epoch
std::string epoch_from_float_doc();
std::string epoch_from_datetime_doc();
std::string epoch_from_string_doc();

// Planet
std::string planet_docstring();
std::string planet_get_name_docstring();
std::string planet_get_extra_info_docstring();
std::string planet_get_mu_central_body_docstring();
std::string planet_get_mu_self_docstring();
std::string planet_get_radius_docstring();
std::string planet_get_safe_radius_docstring();
std::string planet_eph_docstring();
std::string planet_eph_v_docstring();
std::string planet_period_docstring();
std::string planet_elements_docstring();

// UDPLAS
std::string udpla_keplerian_from_elem_docstring();
std::string udpla_keplerian_from_posvel_docstring();
std::string udpla_jpl_lp_docstring();
std::string udpla_vsop2013_docstring();

// Taylor Adaptive propagators
std::string get_stark_docstring();
std::string get_stark_var_docstring();
std::string stark_dyn_docstring();


// Lambert Problem
std::string lambert_problem_docstring();

// Stark problem
std::string stark_problem_docstring();
std::string stark_problem_propagate_docstring();
std::string stark_problem_propagate_var_docstring();


// Propagators
std::string propagate_lagrangian_docstring();
std::string propagate_lagrangian_v_docstring();

// LEG
// Sims Flanagan
std::string leg_sf_docstring();
std::string leg_sf_rvs_docstring();
std::string leg_sf_ms_docstring();
std::string leg_sf_throttles_docstring();
std::string leg_sf_rvf_docstring();
std::string leg_sf_mf_docstring();
std::string leg_sf_tof_docstring();
std::string leg_sf_max_thrust_docstring();
std::string leg_sf_isp_docstring();
std::string leg_sf_mu_docstring();
std::string leg_sf_cut_docstring();
std::string leg_sf_mc_docstring();
std::string leg_sf_tc_docstring();
std::string leg_sf_mc_grad_docstring();
std::string leg_sf_tc_grad_docstring();
std::string leg_sf_nseg_docstring();
std::string leg_sf_nseg_fwd_docstring();
std::string leg_sf_nseg_bck_docstring();

} // namespace pykep

#endif