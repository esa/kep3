// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef KEP3_DETAIL_S11N_HPP
#define KEP3_DETAIL_S11N_HPP

#include <cstddef>
#include <locale>
#include <random>
#include <sstream>
#include <string>
#include <tuple>

#include <boost/serialization/access.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/binary_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

// Include the archives.
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

namespace kep3::detail
{

// Implementation of tuple serialization.
template <std::size_t N>
struct tuple_s11n {
    template <class Archive, typename... Args>
    static void serialize(Archive &ar, std::tuple<Args...> &t, unsigned version)
    {
        ar &std::get<N - 1u>(t);
        tuple_s11n<N - 1u>::serialize(ar, t, version);
    }
};

template <>
struct tuple_s11n<0> {
    template <class Archive, typename... Args>
    static void serialize(Archive &, std::tuple<Args...> &, unsigned)
    {
    }
};

} // namespace kep3::detail

namespace boost::serialization
{

// Implement serialization for std::tuple.
template <class Archive, typename... Args>
inline void serialize(Archive &ar, std::tuple<Args...> &t, unsigned version)
{
    kep3::detail::tuple_s11n<sizeof...(Args)>::serialize(ar, t, version);
}

// Implement serialization for the Mersenne twister engine.
template <class Archive, class UIntType, std::size_t w, std::size_t n, std::size_t m, std::size_t r, UIntType a,
          std::size_t u, UIntType d, std::size_t s, UIntType b, std::size_t t, UIntType c, std::size_t l, UIntType f>
inline void save(Archive &ar, std::mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f> const &e,
                 unsigned)
{
    std::ostringstream oss;
    // Use the "C" locale.
    oss.imbue(std::locale::classic());
    oss << e;
    ar << oss.str();
}

template <class Archive, class UIntType, std::size_t w, std::size_t n, std::size_t m, std::size_t r, UIntType a,
          std::size_t u, UIntType d, std::size_t s, UIntType b, std::size_t t, UIntType c, std::size_t l, UIntType f>
inline void load(Archive &ar, std::mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f> &e,
                 unsigned)
{
    std::istringstream iss;
    // Use the "C" locale.
    iss.imbue(std::locale::classic());
    std::string tmp;
    ar >> tmp;
    iss.str(tmp);
    iss >> e;
}

template <class Archive, class UIntType, std::size_t w, std::size_t n, std::size_t m, std::size_t r, UIntType a,
          std::size_t u, UIntType d, std::size_t s, UIntType b, std::size_t t, UIntType c, std::size_t l, UIntType f>
inline void serialize(Archive &ar, std::mersenne_twister_engine<UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f> &e,
                      unsigned version)
{
    split_free(ar, e, version);
}

} // namespace boost::serialization

#endif