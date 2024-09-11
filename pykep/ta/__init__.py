## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)## 
## This file is part of the kep3 library.## 
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
Various types of interplanetary trajectory legs
"""
from .. import core as _core

# Importing here heyoka we make sure to register the various conversions from/to heyoka objects.
# We then delete the symbol at the end, but the import did its work by then.
import heyoka as _hy 

# Renaming cpp taylor adaptive integrators (we need to create an alias first and then 
# to fool sphinx into thinking these are not aliases, else the sphinx built docs
# would report them as aliases and fail to document these classes)
get_stark = _core._get_stark
get_stark.__module__ = "ta"
get_stark_var = _core._get_stark_var
get_stark_var.__module__ = "ta"
stark_dyn = _core._stark_dyn
stark_dyn.__module__ = "ta"
# Removing core from the list of imported symbols.
del _core
del _hy