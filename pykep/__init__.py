"""Pykep is a coolbox for interplanetary trajectory design developed by ESA's Advanced Concepts Team. Its main
 purpose is fast prototyping of reseacrh ideas, and is not intended for operational usage.

 Some important conventions followed:
 1 - All units expected are S.I. (m,sec,kg,N) unless explicitly stated.
 2 - The default set of osculating orbital parameters is, in this order: [sma, ecc, incl, W, w, f], where f is the true anomaly
 3 - The default option to represent epochs as floats is the modified julian date 2000 (MJD2000). By default, time durations are in days."""

# Version setup.
from ._version import __version__
del _version


# Importing cpp functionalities
from .core import *

# Importing python udplas
from . import udpla
# Importing cpp udplas (we need to create an alias first and then 
# to fool sphinx into thinking these are not aliases, else the sphinx built docs
# would report them as aliases and fail to document these classes)
udpla.keplerian = core._keplerian
udpla.keplerian.__name__ = "keplerian"
udpla.keplerian.__module__ = "udpla"
udpla.null_udpla = core._null_udpla
udpla.null_udpla.__name__ = "null_udpla"
udpla.null_udpla.__module__ = "udpla"
udpla.jpl_lp = core._jpl_lp
udpla.jpl_lp.__name__ = "jpl_lp"
udpla.jpl_lp.__module__ = "udpla"

# Importing the python utils
from .utils import *

# Patch the problem class.
from . import _patch_planet

# Import the plot module
from . import plot

# We import the unit test submodule
from . import test