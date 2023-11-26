"""Pykep is a coolbox for interplanetary trajectory design developed by ESA's Advanced Concepts Team. Its main
 purpose is fast prototyping of reseacrh ideas, and is not intended for operational usage.

 Some important conventions followed:
 1 - All units expected are S.I. (m,sec,kg,N) unless explicitly stated.
 2 - The default set of osculating orbital parameters is, in this order: [sma, ecc, incl, W, w, f], where f is the true anomaly
 3 - The default option to represent epochs as floats is the modified julian date 2000 (MJD2000). By default, time durations are in days."""

# Version setup.
from ._version import __version__
del _version

import os as _os

if _os.name == "posix":
    # NOTE: on some platforms Python by default opens extensions
    # with the RTLD_LOCAL flag, which creates problems because
    # public symbols used by heyoka (e.g., sleef functions, quad
    # precision math) are then not found by the LLVM jit machinery.
    # Thus, before importing core, we temporarily flip on the
    # RTLD_GLOBAL flag, which makes the symbols visible and
    # solves these issues. Another possible approach suggested
    # in the llvm discord is to manually and explicitly add
    # libheyoka.so to the DL search path:
    # DynamicLibrarySearchGenerator::Load(“/path/to/libheyoka.so”)
    # See:
    # https://docs.python.org/3/library/ctypes.html
    import ctypes as _ctypes
    import sys as _sys

    _orig_dlopen_flags = _sys.getdlopenflags()
    _sys.setdlopenflags(_orig_dlopen_flags | _ctypes.RTLD_GLOBAL)

    try:
        # Importing cpp functionalities
        from .core import *
    finally:
        # Restore the original dlopen flags whatever
        # happens.
        _sys.setdlopenflags(_orig_dlopen_flags)

        del _ctypes
        del _sys
        del _orig_dlopen_flags
else:
    # Importing cpp functionalities
    from .core import *

del _os

# Importing user defined planets
from . import udpla
# Renaming cpp udplas (we need to create an alias first and then 
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
udpla.vsop2013 = core._vsop2013
udpla.vsop2013.__name__ = "vsop2013"
udpla.vsop2013.__module__ = "udpla"

# Importing trajectory legs udplas
from . import leg
# Renaming cpp legs (we need to create an alias first and then 
# to fool sphinx into thinking these are not aliases, else the sphinx built docs
# would report them as aliases and fail to document these classes)
leg.sims_flanagan = core._sims_flanagan
udpla.sims_flanagan.__name__ = "sims_flanagan"
udpla.sims_flanagan.__module__ = "leg"

# Importing the python utils
from .utils import *

# Patch the problem class.
from . import _patch_planet

# Import the plot module
from . import plot

# We import the unit test submodule
from . import test