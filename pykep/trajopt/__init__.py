## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)## 
## This file is part of the kep3 library.## 
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
User defined problems (compatible to pagmo) that represent interplanetary optimization problems 
"""

# Direct methods for low-thrust problems
from ._direct_point2point import direct_point2point
from ._direct_pl2pl import direct_pl2pl

# Evolutionary encodings for high energy transfers (chemical propulsion)
from ._mga import mga

# The interplanetary trajectory gym
from . import gym 