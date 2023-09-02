
// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_UDPLA_KEPLERIAN_H
#define kep3_UDPLA_KEPLERIAN_H

#include <array>
#include <utility>

#include <fmt/ostream.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/detail/visibility.hpp>
#include <kep3/epoch.hpp>
#include <kep3/planet.hpp>

namespace kep3::udpla {

class kep3_DLL_PUBLIC keplerian {

  kep3::epoch m_ref_epoch;
  std::array<std::array<double, 3>, 2> m_pos_vel_0;
  std::string m_name;
  double m_mu_central_body;
  double m_mu_self;
  double m_radius;
  double m_safe_radius;
  bool m_ellipse;

  friend class boost::serialization::access;
  template <typename Archive> void serialize(Archive &ar, unsigned) {
    ar &m_ref_epoch;
    ar &m_pos_vel_0;
    ar &m_name;
    ar &m_mu_central_body;
    ar &m_mu_self;
    ar &m_radius;
    ar &m_safe_radius;
    ar &m_ellipse;
  }

public:
  // NOTE: in here elem is a,e,i,W,w,M (Mean anomaly, not true anomaly)
  // NOTE: added_param contains mu_self, radius and safe_radius
  explicit keplerian(const epoch &ref_epoch, const std::array<double, 6> &par,
                     double mu_central_body = 1., std::string name = "Unknown",
                     std::array<double, 3> added_params = {-1., -1., -1.});
  explicit keplerian(const epoch &ref_epoch = kep3::epoch(0),
                     const std::array<std::array<double, 3>, 2> &pos_vel =
                         {{{1.0, 0.0, 0.0}, {0., 1.0, 0.0}}},
                     double mu_central_body = 1., std::string name = "Unknown",
                     std::array<double, 3> added_params = {-1., -1., -1.});
  // Mandatory UDPLA methods
  [[nodiscard]] std::array<std::array<double, 3>, 2> eph(const epoch &) const;

  // Optional UDPLA methods
  [[nodiscard]] std::string get_name() const;
  [[nodiscard]] double get_mu_central_body() const;
  [[nodiscard]] double get_mu_self() const;
  [[nodiscard]] double get_radius() const;
  [[nodiscard]] double get_safe_radius() const;
  [[nodiscard]] std::string get_extra_info() const;

  // Other methods
  [[nodiscard]] kep3::epoch get_ref_epoch() const;
  [[nodiscard]] std::array<double, 6>
      elements(kep3::elements_type = kep3::elements_type::KEP_F) const;
};
kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &,
                                         const kep3::udpla::keplerian &);
} // namespace kep3::udpla

template <>
struct fmt::formatter<kep3::udpla::keplerian> : ostream_formatter {};
kep3_S11N_PLANET_EXPORT_KEY(kep3::udpla::keplerian);

#endif // kep3_EPOCH_H