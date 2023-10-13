"""Pykep is a coolbox for interplanetary trajectory design developed by ESA's Advanced Concepts Team. Its main
 purpose is fast prototyping of reseacrh ideas, and is not intended for operational usage.

 Some important conventions followed:
 1 - All units expected are S.I. (m,sec,kg,N) unless explicitly stated.
 2 - The default set of osculating orbital parameters is, in this order: [sma, ecc, incl, W, w, f], where f is the true anomaly
 3 - The default option to represent epochs as floats is the modified julian date 2000 (MJD2000)."""

# Version setup.
from ._version import __version__
del _version


# Importing cpp functionalities
from .core import *

from .udpla import tle_satellite

# Patch the problem class.
from . import _patch_planet

# We import the unit test submodule
from . import test