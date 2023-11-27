## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)## 
## This file is part of the kep3 library.## 
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
User defined planets that can construct a pykep.planet
"""

from ._tle import tle
from ._spice import spice, de440s
from .. import core

# Renaming cpp udplas (we need to create an alias first and then 
# to fool sphinx into thinking these are not aliases, else the sphinx built docs
# would report them as aliases and fail to document these classes)
keplerian = core._keplerian
keplerian.__name__ = "keplerian"
keplerian.__module__ = "udpla"
null_udpla = core._null_udpla
null_udpla.__name__ = "null_udpla"
null_udpla.__module__ = "udpla"
jpl_lp = core._jpl_lp
jpl_lp.__name__ = "jpl_lp"
jpl_lp.__module__ = "udpla"
vsop2013 = core._vsop2013
vsop2013.__name__ = "vsop2013"
vsop2013.__module__ = "udpla"

# Removing core from the list of imported symbols.
del core
