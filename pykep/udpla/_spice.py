## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)## 
## This file is part of the kep3 library.## 
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import spiceypy as pyspice
import pykep as pk
import numpy as np


class spice:
    """__init__(body, ref_frame, obs)

    This User Defined Planet (UDPLA) represents a planet/object/spacecraft as defined by a pre-loaded
    SPICE kernel. The interface to the NAIF SPICE code is provided via the third-party python package called
    spiceypy: https://spiceypy.readthedocs.io/en/stable/

    The resulting ephemerides will be returned in SI units and in the selected reference frame (J2000 ECLIPTIC by default).

    .. note::
       SPICE kernels containing interesting objects can be, for example, found at:
       https://naif.jpl.nasa.gov/naif/data_archived.html -> Spacecarft from NAIF
       https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/ -> Solar System planets
       https://www.cosmos.esa.int/web/spice -> ESA missions

    .. warning::
       The utility function :func:`~pykep.utils.load_spice_kernels` must be called to pre-load into
       memory the necessary SPICE kernels. Else, the call to the :func:`~pykep.udpla.spice.eph` method will
       fail rising an exception.

    """
    def __init__(self, body, ref_frame, obs):
        # We check if the various names are supported by NAIF.
        # If not the conversions will fail and rise and exception.
        # This may be too restrictive and maybe lifted in the future.
        self.naifid = pk.utils.name2naifid(body)
        self.frame_naifid = pk.utils.framename2naifid(ref_frame)
        self.obs_naifid = pk.utils.name2naifid(obs)
        # Store the strings
        self.body = body
        self.ref_frame = ref_frame
        self.obs = obs
        # And we also store the default values for physical parameters 
        # to allow the user to overwrite them
        self.mu_central_body = -1
        self.mu_self = -1
        self.radius = -1
        self.safe_radius = -1

    def eph(self, mjd2000):
        """Mandatory method of the :class:`~pykep.planet` interface.

        Args:
            *mjd2000* (:class:`float`): Modified Julian Date 2000

        Returns:
            :class:`list` [:class:`list`, :class:`list`]: the planet ephemerides.
        """
        spice_epoch = (mjd2000-0.5)*pk.DAY2SEC
        rv, _ = pyspice.spkezr(self.body, spice_epoch, self.ref_frame, "NONE", self.obs)
        return [rv[:3] * 1000, rv[3:] * 1000]
    
    def eph_v(self, mjd2000s):
        spice_epochs = (np.array(mjd2000s) - 0.5)*pk.DAY2SEC
        rvs, _ = pyspice.spkezr(self.body, spice_epochs, self.ref_frame, "NONE", self.obs)
        return np.array(rvs) * 1000

    def get_name(self):
        """Optional method of the :class:`~pykep.planet` interface.

        Returns:
            :class:`str`: The body name
        """
        return self.body + " - SPICE"
    
    def get_extra_info(self):
        """Optional method of the :class:`~pykep.planet` interface.

        Returns:
            :class:`str`: Extra info on the udpla
        """
        return "Body: " + self.body + "\nObserver: " + self.obs + "\nReference Frame: " + self.ref_frame 
    
    def get_mu_central_body(self):
        return self.mu_central_body
    
    def get_mu_self(self):
        return self.mu_self
    
    def get_radius(self):
        return self.radius
    
    def get_safe_radius(self):
        return self.safe_radius
    

class de440s(spice):
    """__init__(body = "EARTH BARYCENTER", ref_frame = "ECLIPJ2000", obs = "SSB")

    This User Defined Planet (UDPLA) represents a solar system planet as described by JPL 440s ephemerides.
    The correct spice kernel is preloaded upon instantiation of the class as its shipped with the pykep module.

    Args:
        *body* (:class:`str`): NAIF name of the solar system body. Defaults to: "EARTH BARYCENTER".

        *ref_frame* (:class:`str`): NAIF name of the reference frame. For example: "ECLIPJ2000".

        *obs* (:class:`str`): NAIF name of the observer. For example: "SSB" (solar system barycenter).
    """
    def __init__(self, body = "EARTH BARYCENTER", ref_frame = "ECLIPJ2000", obs = "SSB"):
        path = pk.__path__[0]
        self.kernel = path+"/data/de440s.bsp"
        pk.utils.load_spice_kernels(self.kernel)
        spice.__init__(self, body, ref_frame = ref_frame, obs = obs)

    def __del__(self):
        """Destructor.
       Unloads the kernel from memory as not needed anymore.
        """
        pk.utils.unload_spice_kernels(self.kernel)

    def get_name(self):
        """Optional method of the :class:`~pykep.planet` interface.

        Returns:
            :class:`str`: The body name
        """
        return self.body + " - de440s"
    
    def kernel_file():
        """The full path for the kernel file "de440s.bsp" shipped in `pykep`

        Returns:
            :class:`str`: The full file name.
        """
        path = pk.__path__[0] + "/data/de440s.bsp"
        return path
    
    def body_list():
        """The list of bodies contained in the de440s ephemerides. 
        This is a static method and can be quieried before constructing the object.

        Returns:
            :class:`list`: The list of possible body names.
        """
        kernel_file = pk.udpla.de440s.kernel_file()
        list = pk.utils.inspect_spice_kernel(kernel_file)
        return [pk.utils.naifid2name(item) for item in list]
    
