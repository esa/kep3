// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <stdexcept>
#include <vector>
#include <kep3/leg/sf_checks.hpp>

namespace kep3::leg {

void _check_tof(double tof)
{
    // SC: One should be able to give this as a negative number to run the system backwards, no?
    if (tof < 0.) {
        throw std::domain_error("The time of flight of a sims_flanagan leg needs to be larger or equal to zero.");
    }
}
void _check_throttles(const std::vector<double> &throttles)
{
    if ((throttles.size() % 3) != 0u) {
        throw std::logic_error("The throttles of a sims_flanagan leg are detected to be not a multiple of 3 in size "
                               "[u0x, u0y, u0z, .....].");
    }
    if (throttles.empty()) {
        throw std::logic_error(
            "The throttles of a sims_flanagan leg are detected to be empty! At least one segment is necessary.");
    }
}
void _check_max_thrust(double max_thrust)
{
    if (max_thrust < 0.) {
        throw std::domain_error(
            "The maximum allowed thrust of a sims_flanagan leg is detected to be smaller than zero.");
    }
}
void _check_isp(double isp)
{
    if (isp < 0.) {
        throw std::domain_error("The specific impulse of a sims_flanagan leg is detected to be smaller than zero.");
    }
}
void _check_mu(double mu)
{
    if (mu < 0.) {
        throw std::domain_error(
            "The gravitational parameter of a sims_flanagan leg is detected to be smaller than zero.");
    }
}
void _check_cut(double cut)
{
    if (cut < 0. || cut > 1.) {
        throw std::domain_error("The parameter cut of a sims_flanagan leg must be in [0, 1].");
    }
}
void _check_tol(double tol)
{
    if (tol <= 0. || tol > 1.) {
        throw std::domain_error("The parameter tol of a high-fidelity sims-flanagan leg leg must be in <0, 1].");
    }
}
void _check_nseg(unsigned nseg, unsigned nseg_fwd, unsigned nseg_bck)
{
    if (nseg_fwd + nseg_bck != nseg)
    {
        throw std::logic_error("The number of segments provided does not add up.");
    }
}
void _sanity_checks(const std::vector<double> &throttles, double tof, double max_thrust, double isp, double mu,
                    double cut, unsigned nseg, unsigned nseg_fwd, unsigned nseg_bck)
{
    _check_throttles(throttles);
    _check_tof(tof);
    _check_max_thrust(max_thrust);
    _check_isp(isp);
    _check_mu(mu);
    _check_cut(cut);
    _check_nseg(nseg, nseg_fwd, nseg_bck);
}
void _sanity_checks(const std::vector<double> &throttles, double tof, double max_thrust, double isp, double mu,
                    double cut, double tol, unsigned nseg, unsigned nseg_fwd, unsigned nseg_bck)
{
    _check_throttles(throttles);
    _check_tof(tof);
    _check_max_thrust(max_thrust);
    _check_isp(isp);
    _check_mu(mu);
    _check_cut(cut);
    _check_tol(tol);
    _check_nseg(nseg, nseg_fwd, nseg_bck);
}

} // namespace kep3::leg