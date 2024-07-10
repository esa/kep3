// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_LEG_SIMS_FLANAGAN_HF_H
#define kep3_LEG_SIMS_FLANAGAN_HF_H

#include <array>
#include <tuple>
#include <vector>

#include <fmt/ostream.h>

#include <heyoka/taylor.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/detail/visibility.hpp>
#include <kep3/epoch.hpp>

namespace kep3::leg
{
/// The High-fidelity leg model expanding upon Sims-Flanagan
/**
 * This class represents, generically, a low-thrust leg as a sequence of successive
 * fixed low-thrust segments.
 * The leg achieves to transfer a given spacecraft from an initial to a final state in the
 * time given (and can be considered as feasible) whenever the method evaluate_mismatch
 * returns all zeros and the method get_throttles_con returns all values less than zero.
 *
 * The dynamics is by default Keplerian, but the user can construct his own and pass it to the
 * leg. The user constructed dynamics will need to meet the requirements of having the variables
 * x,y,z,vx,vy,vz,m and the first three parameters heyoka::par[0], heyoka::par[1], heyoka::par[2] that represent
 * the thrust direction.
 */
class kep3_DLL_PUBLIC sims_flanagan_hf
{
public:
    // Default Constructor.
    sims_flanagan_hf() = default;
    // Constructors
    sims_flanagan_hf(
        const std::array<std::array<double, 3>, 2> &rvs, double ms, std::vector<double> throttles,
        const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust, double isp,
        double mu, double cut, double ts,
        const std::optional<const std::pair<heyoka::taylor_adaptive<double> &, heyoka::taylor_adaptive<double> &>>
            &tas);

    // Setters
    void set_tof(double tof);
    void set_rvs(std::array<std::array<double, 3>, 2> rv);
    void set_ms(double mass);
    void set_throttles(std::vector<double> throttles);
    void set_throttles(std::vector<double>::const_iterator it1, std::vector<double>::const_iterator it2);
    void set_rvf(std::array<std::array<double, 3>, 2> rv);
    void set_mf(double mass);
    void set_max_thrust(double max_thrust);
    void set_isp(double isp);
    void set_mu(double mu);
    void set_cut(double cut);
    void set_ts(double ts);
    void set(const std::array<std::array<double, 3>, 2> &rvs, double ms, const std::vector<double> &throttles,
             const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust, double isp,
             double mu, double cut = 0.5);

    // Getters
    [[nodiscard]] double get_tof() const;
    [[nodiscard]] const std::array<std::array<double, 3>, 2> &get_rvs() const;
    [[nodiscard]] double get_ms() const;
    [[nodiscard]] const std::vector<double> &get_throttles() const;
    [[nodiscard]] const std::array<std::array<double, 3>, 2> &get_rvf() const;
    [[nodiscard]] double get_mf() const;
    [[nodiscard]] double get_max_thrust() const;
    [[nodiscard]] double get_isp() const;
    [[nodiscard]] double get_mu() const;
    [[nodiscard]] double get_cut() const;
    [[nodiscard]] double get_ts() const;
    [[nodiscard]] unsigned get_nseg() const;
    [[nodiscard]] unsigned get_nseg_fwd() const;
    [[nodiscard]] unsigned get_nseg_bck() const;

    // Compute constraints
    [[nodiscard]] std::array<double, 7> compute_mismatch_constraints() const;
    [[nodiscard]] std::vector<double> compute_throttle_constraints() const;

    // Compute mismatch constraint gradients (w.r.t. rvm state and w.r.t. throttles, tof)
    [[nodiscard]] std::tuple<std::array<double, 49>, std::array<double, 49>, std::vector<double>>
    compute_mc_grad() const;

    // Compute throttle constraint gradients
    [[nodiscard]] std::vector<double> compute_tc_grad() const;

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar &m_rvs;
        ar &m_ms;
        ar &m_throttles;
        ar &m_tof;
        ar &m_rvf;
        ar &m_mf;
        ar &m_max_thrust;
        ar &m_isp;
        ar &m_mu;
        ar &m_cut;
        ar &m_ts;
        ar &m_nseg;
        ar &m_nseg_fwd;
        ar &m_nseg_bck;
        ar &m_tas;
    }

    // Initial spacecraft state.
    std::array<std::array<double, 3>, 2> m_rvs{{{1., 0., 0.}, {0, 1., 0.}}};
    double m_ms = 1.;
    // Sequence of throttles.
    std::vector<double> m_throttles{0., .0, 0., 0., 0., 0.};
    // Final spacecraft state.
    std::array<std::array<double, 3>, 2> m_rvf{{{0., 1., 0.}, {-1., 0., 0.}}};
    double m_mf = 1.;
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
    // The reference epoch
    double m_ts = 0.;
    // The adaptive Taylor integrators
    std::optional<std::pair<heyoka::taylor_adaptive<double>, heyoka::taylor_adaptive<double>>> m_tas = std::nullopt;
    // Segment sizes
    unsigned m_nseg = 2u;
    unsigned m_nseg_fwd = 1u;
    unsigned m_nseg_bck = 1u;
};

// Streaming operator for the class kep3::leg::sims_flanagan.
kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &, const sims_flanagan_hf &);

} // namespace kep3::leg

template <>
struct fmt::formatter<kep3::leg::sims_flanagan_hf> : fmt::ostream_formatter {
};

#endif