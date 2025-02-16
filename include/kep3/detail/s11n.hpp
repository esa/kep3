// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

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

#endif