// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <kep3/planet.hpp>

namespace kep3::detail {

planet_inner_base::planet_inner_base() = default;
planet_inner_base::~planet_inner_base() = default;
null_udpla::null_udpla() = default;

std::array<std::array<double, 3>, 2> null_udpla::eph(const epoch &) {
  std::array<double, 3> pos = {1., 0., 0.};
  std::array<double, 3> vel = {0., 1., 0.};
  return {pos, vel};
};
} // namespace kep3::detail

namespace kep3 {
planet::planet() : planet(detail::null_udpla{}){};
planet::planet(const planet &other) : m_ptr(other.m_ptr->clone()){};
planet::planet(planet &&other) noexcept : m_ptr(std::move(other.m_ptr)){};

planet &planet::operator=(planet &&other) noexcept {
  if (this != &other) {
    m_ptr = std::move(other.m_ptr);
  }
  return *this;
}

planet &planet::operator=(const planet &other) { return *this = planet(other); }

bool planet::is_valid() const { return static_cast<bool>(m_ptr); }

std::type_index planet::get_type_index() const {
  return ptr()->get_type_index();
}

const void *planet::get_ptr() const { return ptr()->get_ptr(); }

void *planet::get_ptr() { return ptr()->get_ptr(); }

detail::planet_inner_base const *planet::ptr() const {
  assert(m_ptr.get() != nullptr);
  return m_ptr.get();
}

detail::planet_inner_base *planet::ptr() {
  assert(m_ptr.get() != nullptr);
  return m_ptr.get();
}

std::array<std::array<double, 3>, 2> planet::eph(const epoch &ep) {
  return ptr()->eph(ep);
}

std::string planet::get_name() const { return ptr()->get_name(); }

std::string planet::get_extra_info() const { return ptr()->get_extra_info(); }

} // namespace kep3