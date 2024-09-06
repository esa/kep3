// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_STARK_PROBLEM_H
#define kep3_STARK_PROBLEM_H

#include <array>
#include <utility>
#include <vector>

#include <fmt/ostream.h>

#include <heyoka/taylor.hpp>

#include <kep3/detail/s11n.hpp>
#include <kep3/detail/visibility.hpp>

namespace kep3
{

/// Stark Problem
/**
 * This class represent the Stark's problem. When instantiated it gets a
 * Taylor adaptive integrator from the kep3 ta cache and copies it as to be able to use it
 * later for numerical propagation. The information on the solution sensitivities can be retreived too
 * in which case a variational integrator is employed.
 */
class kep3_DLL_PUBLIC stark_problem
{
public:
    explicit stark_problem(double mu=1., double veff=1., double tol=1e-16);
    std::array<double, 7> propagate(const std::array<double, 7> &x0, std::array<double, 3> thrust, double tof);
    std::tuple<std::array<double, 7>, std::array<double, 49>, std::array<double, 21>>
    propagate_var(const std::array<double, 7> &x0, std::array<double, 3> thrust, double tof);
    // Getters
    [[nodiscard]] double get_mu() const;
    [[nodiscard]] double get_veff() const;
    [[nodiscard]] double get_tol() const;
    [[nodiscard]] heyoka::taylor_adaptive<double> get_ta() const;
    [[nodiscard]] heyoka::taylor_adaptive<double> get_ta_var() const;
    // Setters
    void set_mu(double);
    void set_veff(double);

private:
    // Serialization code
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar &m_mu;
        ar &m_veff;
        ar &m_tol;
        ar &m_ta;
        ar &m_ta_var;
        ar &m_var_ic;
    }
    double m_mu;
    double m_veff;
    double m_tol;
    heyoka::taylor_adaptive<double> m_ta;
    heyoka::taylor_adaptive<double> m_ta_var;
    std::vector<double> m_var_ic;
};
kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &, const stark_problem &);
} // namespace kep3

template <>
struct fmt::formatter<kep3::stark_problem> : fmt::ostream_formatter {
};

#endif // kep3_STARK_PROBLEM_H