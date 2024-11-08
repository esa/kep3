// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_SF_CHECKS_H
#define kep3_SF_CHECKS_H

#include <vector>

// These checks are used for the low- and high-fidelity legs (in sims_flanagan.cpp and sims_flanagan_hf.cpp)

void _check_tof(double tof);
void _check_throttles(const std::vector<double> &throttles);
void _check_max_thrust(double max_thrust);
void _check_isp(double isp);
void _check_mu(double mu);
void _check_cut(double cut);
void _check_tol(double tol);
void _sanity_checks(const std::vector<double> &throttles, double tof, double max_thrust, double isp, double mu,
                    double cut, unsigned nseg, unsigned nseg_fwd, unsigned nseg_bck);
void _sanity_checks(const std::vector<double> &throttles, double tof, double max_thrust, double isp, double mu,
                    double cut, double tol, unsigned nseg, unsigned nseg_fwd, unsigned nseg_bck);

#endif