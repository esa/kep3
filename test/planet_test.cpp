// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.


#include <string>

#include <fmt/core.h>

#include <kep3/planet.hpp>

#include "catch.hpp"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>


// IMPORTANT: put the IMPLEMENT here (in the executable), not in the dylib
TANUKI_S11N_WRAP_EXPORT_IMPLEMENT(kep3::detail::null_udpla, kep3::detail::planet_iface)


using kep3::planet;


TEST_CASE("serialization_test_2")
{
    // Instantiate a planet
    using udp_type = kep3::detail::null_udpla;
    planet pla{kep3::detail::null_udpla()};
    // Store the string representation.
    std::stringstream ss;
    auto before = boost::lexical_cast<std::string>(pla);
    // Now serialize, deserialize and compare the result.
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << pla;
    }
}