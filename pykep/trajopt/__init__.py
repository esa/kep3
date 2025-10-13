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
from ._direct_cr3bp import direct_cr3bp
from ._direct_pl2pl_alpha import direct_pl2pl_alpha
from ._direct_cr3bp_alpha import direct_cr3bp_alpha

# Evolutionary encodings for high energy transfers (chemical propulsion)
from ._mga import mga
from ._mga_1dsm import mga_1dsm
from ._pl2pl_N_impulses import pl2pl_N_impulses

# Indirect methods for low-thrust problems
from ._pontryagin_cartesian import pontryagin_cartesian_mass, pontryagin_cartesian_time
from ._pontryagin_equinoctial import pontryagin_equinoctial_mass, pontryagin_equinoctial_time
from ._mim import mim_from_hop

# MIT (multiple Impulse Trajectories)
from ._primer_vector import primer_vector, primer_vector_surrogate
from ._min_Bu_bu import minBu_bu_p, minBu_bu

# The launchers models
from ._launchers import _launchers
launchers = _launchers()

# The interplanetary trajectory gym
from . import gym