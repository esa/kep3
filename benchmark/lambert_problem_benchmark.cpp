/*****************************************************************************
 *   Copyright (C) 2004-2018 The pykep development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://gitter.im/esa/pykep                                             *
 *   https://github.com/esa/pykep                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#include <iomanip>
#include <iostream>
#include <random>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/lambert_problem.hpp>
#include <stdexcept>

int main() {
  // Preamble
  std::array<double, 3> r1{{0, 0, 0}}, r2{{0, 0, 0}};
  double tof = 0.;
  // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
  std::mt19937 rng_engine(122012203u);
  std::uniform_int_distribution<unsigned> dist(0, 1);
  std::uniform_real_distribution<double> dist1(-2, 2);

  double acc = 0, err_max = 0;
  unsigned count = 0u;

  // Experiment Settings
  unsigned int Ntrials = 120000;

  // Start Experiment
  for (unsigned int i = 0; i < Ntrials; ++i) {
    // 1 - generate a random problem geometry
    r1[0] = dist1(rng_engine);
    r1[1] = dist1(rng_engine);
    r1[2] = dist1(rng_engine);
    r2[0] = dist1(rng_engine);
    r2[1] = dist1(rng_engine);
    r2[2] = dist1(rng_engine);
    tof = (dist1(rng_engine) + 2) / 4 * 100 + 0.1;

    // 2 - Solve the lambert problem
    double mu = 1.0;
    unsigned revs_max = 20;
    bool cw = static_cast<bool>(dist(rng_engine));
    kep3::lambert_problem lp(r1, r2, tof, mu, cw, revs_max);

    // 3 - Check its precision using propagate_lagrangian
    for (const auto &v1 : lp.get_v1()) {
      std::array<std::array<double, 3>, 2> pos_vel{{r1, v1}};
      kep3::propagate_lagrangian(pos_vel, tof, mu);
      double err = std::sqrt((pos_vel[0][0] - r2[0]) * (pos_vel[0][0] - r2[0]) +
                             (pos_vel[0][1] - r2[1]) * (pos_vel[0][1] - r2[1]) +
                             (pos_vel[0][2] - r2[2]) * (pos_vel[0][2] - r2[2]));
      err_max = std::max(err_max, err);
      acc += err;
    }
    count += (lp.get_Nmax() * 2 + 1);
  }
  std::cout << "Max error: " << err_max << std::endl;
  std::cout << "Average Error: " << acc / count << std::endl;
  std::cout << "Number of Problems Solved: " << count << std::endl;
}