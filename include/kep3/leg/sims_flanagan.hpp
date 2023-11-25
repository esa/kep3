// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_LEG_SIMS_SLANAGAN_H
#define kep3_LEG_SIMS_SLANAGAN_H

#include "kep3/core_astro/constants.hpp"
#include <array>
#include <vector>

#include <fmt/ostream.h>

#include <kep3/detail/visibility.hpp>
#include <kep3/epoch.hpp>

namespace kep3::leg
{
/// The Sims-Flanagan leg model
/**
 * This class represents, generically, a low-thrust leg as a sequence of successive
 * impulses of magnitude compatible with the low-thrust propulsion system of a spacecraft.
 * The leg achieves to transfer a given spacecraft from an initial to a final state in the
 * time given (and can be considered as feasible) whenever the method evaluate_mismatch
 * returns all zeros and the method get_throttles_con returns all values less than zero.
 */
class kep3_DLL_PUBLIC sims_flanagan
{
public:
    // Default Constructor.
    sims_flanagan() = default;
    // Constructors
    sims_flanagan(const std::array<double, 7> &xs, std::vector<double> throttles, const std::array<double, 7> &xf,
                  double tof, double max_thrust, double isp, double mu, double cut = 0.5);

    // Setters
    void set_tof(double tof);
    void set_xs(std::array<double, 7> xs);
    void set_throttles(std::vector<double> throttles);
    void set_xf(std::array<double, 7> xf);
    void set_max_thrust(double max_thrust);
    void set_isp(double isp);
    void set_mu(double mu);
    void set_cut(double cut);
    void set(const std::array<double, 7> &xs, const std::vector<double>& throttles, const std::array<double, 7> &xf,
             double tof, double max_thrust, double isp, double mu, double cut = 0.5);

    // Getters
    [[nodiscard]] double get_tof() const;
    [[nodiscard]] const std::array<double, 7> &get_xs() const;
    [[nodiscard]] const std::vector<double> &get_throttles() const;
    [[nodiscard]] const std::array<double, 7> &get_xf() const;
    [[nodiscard]] double get_max_thrust() const;
    [[nodiscard]] double get_isp() const;
    [[nodiscard]] double get_mu() const;
    [[nodiscard]] double get_cut() const;

    // Compute constraints
    [[nodiscard]] std::array<double, 7> compute_mismatch_constraints() const;
    [[nodiscard]] std::vector<double> compute_throttle_constraints() const;

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar &m_xs;
        ar &m_throttles;
        ar &m_tof;
        ar &m_xf;
        ar &m_max_thrust;
        ar &m_isp;
        ar &m_mu;
        ar &m_cut;
    }

    // Initial spacecraft state.
    std::array<double, 7> m_xs{1., 0., 0., 0, 1., 0., 1.};
    // Sequence of throttles.
    std::vector<double> m_throttles{0., .0, 0., 0., 0., 0.};
    // Final spacecraft state.
    std::array<double, 7> m_xf{0., 1., 0., -1., 0., 0., 1.};
    // Time of flight (defaults to 1/4 of the period)
    double m_tof = kep3::pi / 2;
    // Spacecraft propulsion system maximum thrust.
    double m_max_thrust{1.};
    // Spacecraft propulsion system specific impulse.
    double m_isp{1.};
    // Spacecraft gravitational parameter.
    double m_mu{1.};
    // The cut parameter
    double m_cut = 0.5;
};

// Streaming operator for the class kep3::leg::sims_flanagan.
kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &, const sims_flanagan &);

} // namespace kep3::leg

template <>
struct fmt::formatter<kep3::leg::sims_flanagan> : fmt::ostream_formatter {
};

#endif