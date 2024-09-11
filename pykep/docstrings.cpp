// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <string>

#include "docstrings.hpp"

namespace pykep
{

std::string core_module_doc()
{
    return R"(core is the Pykep module that contains most of its core routines efficiently coded in c++
)";
}

std::string m2e_doc()
{
    return R"(m2e(M, ecc)
    
    Converts from Mean to Eccentric anomaly. Requires ecc < 1.

    Args:
      *M* (:class:`float`): the Mean anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Eccentric anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> M = 1.2
      >>> ecc = 0.1
      >>> pk.m2e(M, ecc)
      1.296254963787226
)";
}

std::string e2m_doc()
{
    return R"(e2m(E, ecc)
    
    Converts from Eccentric to Mean anomaly. Requires ecc < 1.

    Args:
      *E* (:class:`float`): the Eccentric anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Mean anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> E = 0.5
      >>> ecc = 0.1
      >>> pk.e2m(E, ecc)
      0.4520574461395797
)";
}

std::string e2f_doc()
{
    return R"(e2f(E, ecc)
    
    Converts from eccentric to true anomaly. Requires ecc < 1.

    Args:
      *E* (:class:`float`): the Eccentric anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> E = 0.5
      >>> ecc = 0.1
      >>> pk.e2f(E, ecc)
      0.5502639747136633
)";
}

std::string f2e_doc()
{
    return R"(f2e(f, ecc)
    
    Converts from True to Eccentric anomaly. Requires ecc < 1.

    Args:
      *f* (:class:`float`): the True anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Eccentric anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> f = 1.2
      >>> ecc = 0.1
      >>> pk.f2e(f, ecc)
      1.1082931139529482
)";
}

std::string f2m_doc()
{
    return R"(f2m(f, ecc)
    
    Converts from True to Mean anomaly. Requires ecc < 1.

    Args:
      *f* (:class:`float`): the True anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Mean anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> f = -0.34
      >>> ecc = 0.67
      >>> pk.f2m(f, ecc)
      -0.05065883735669101
)";
}

std::string m2f_doc()
{
    return R"(m2f(M, ecc)
    
    Converts from Mean to True anomaly. Requires ecc < 1.

    Args:
      *M* (:class:`float`): the Mean anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> M = 0.32
      >>> ecc = 0.65
      >>> pk.m2f(M, ecc)
      1.4497431281728277
)";
}

std::string h2n_doc()
{
    return R"(h2n(H, ecc)
    
    Converts from Hyperbolic to Hyperbolic Mean anomaly. Requires ecc > 1.

    Args:
      *H* (:class:`float`): the Hyperbolic anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Hyperbolic Mean anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> H = 1.2
      >>> ecc = 10.32
      >>> pk.h2n(H, ecc)
      14.377641187853621
)";
}

std::string n2h_doc()
{
    return R"(n2h(N, ecc)
    
    Converts from Hyperbolic Mean to Hyperbolic anomaly. Requires ecc > 1.

    Args:
      *N* (:class:`float`): the Hyperbolic Mean anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Hyperbolic anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> N = 1.2
      >>> ecc = 10.32
      >>> pk.n2h(N, ecc)
      0.12836469743916526
)";
}

std::string h2f_doc()
{
    return R"(h2f(H, ecc)
    
    Converts from Hyperbolic to True anomaly. Requires ecc > 1.

    Args:
      *H* (:class:`float`): the Hyperbolic anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> H = 10.32
      >>> ecc = 4.5
      >>> pk.h2f(H, ecc)
      1.7948251330114304
)";
}

std::string f2h_doc()
{
    return R"(f2h(f, ecc)
    
    Converts from True to Hyperbolic anomaly. Requires ecc > 1.

    Args:
      *f* (:class:`float`): the True anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Hyperbolic anomaly

    Examples:
      >>> import pykep as pk
      >>> f = 1.2
      >>> ecc = 1.1
      >>> pk.f2h(f, ecc)
      0.30083016696826936
)";
}

std::string f2n_doc()
{
    return R"(f2n(f, ecc)
    
    Converts from True to Hyperbolic Mean anomaly. Requires ecc > 1.

    Args:
      *f* (:class:`float`): the True anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Hyperbolic Mean anomaly

    Examples:
      >>> import pykep as pk
      >>> f = 1.2
      >>> ecc = 5.7
      >>> pk.f2n(f, ecc)
      8.421335633880908
)";
}

std::string n2f_doc()
{
    return R"(n2f(N, ecc)
    
    Converts from Hyperbolic Mean to True anomaly. Requires ecc > 1.

    Args:
      *N* (:class:`float`): the Hyperbolic Mean anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the True anomaly

    Examples:
      >>> import pykep as pk
      >>> N = 10.32
      >>> ecc = 13.45
      >>> pk.n2f(N, ecc)
      0.7373697968359353
)";
}

std::string zeta2f_doc()
{
    return R"(zeta2f(zeta, ecc)
    
    Converts from Gudermannian to True anomaly. Requires ecc > 1.

    See Battin: "An Introduction to the Mathematics and Methods of Astrodynamics" for a 
    definition of zeta and the treatment of the resulting equations.

    Args:
      *zeta* (:class:`float`): the Gudermannian (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the True anomaly

    Examples:
      >>> import pykep as pk
      >>> zeta = 8.2
      >>> ecc = 2.2
      >>> pk.zeta2f(zeta, ecc)
      2.3290929552114266
)";
}

std::string f2zeta_doc()
{
    return R"(f2zeta(f, ecc)
    
    Converts from True anomaly to Gudermannian. Requires ecc > 1.

    Args:
      *f* (:class:`float`): the True anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Gudermannian 

    Examples:
      >>> import pykep as pk
      >>> f = 0.5
      >>> ecc = 3.3
      >>> pk.f2zeta(f, ecc)
      0.36923933496389816
)";
}

std::string m2e_v_doc()
{
    return R"(m2e_v(Ms, eccs)
    
    Converts from Mean to Eccentric anomaly (vectorized version). Requires ecc < 1.

    Args:
      *Ms* (:class:`numpy.ndarray` or :class:`float`): the Mean anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Eccentric anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Ms = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 0.375
      >>> Es = pk.m2e_v(Ms, ecc)
      >>> np.shape(Es)
      (100,)
)";
}

std::string e2m_v_doc()
{
    return R"(e2m_v(Es, eccs)
    
    Converts from Eccentric to Mean anomaly (vectorized version). Requires ecc < 1.

    Args:
      *Es* (:class:`numpy.ndarray` or :class:`float`): the Eccentric anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Mean anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Es = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 0.86345
      >>> Ms = pk.e2m_v(Es, ecc)
      >>> np.shape(Ms)
      (100,)
)";
}

std::string e2f_v_doc()
{
    return R"(e2f_v(Es, eccs)
    
    Converts from eccentric to true anomaly (vectorized version). Requires ecc < 1.

    Args:
      *Es* (:class:`numpy.ndarray` or :class:`float`): the Eccentric anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Es = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 0.0256
      >>> fs = pk.e2f_v(Es, ecc)
      >>> np.shape(fs)
      (100,)
)";
}

std::string f2e_v_doc()
{
    return R"(f2e_v(fs, eccs)
    
    Converts from True to Eccentric anomaly (vectorized version). Requires ecc < 1.

    Args:
      *fs* (:class:`numpy.ndarray` or :class:`float`): the True anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Eccentric anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> fs = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 0.23
      >>> Es = pk.f2e_v(fs, ecc)
      >>> np.shape(Es)
      (100,)
)";
}

std::string f2m_v_doc()
{
    return R"(f2m_v(fs, eccs)
    
    Converts from True to Mean anomaly (vectorized version). Requires ecc < 1.

    Args:
      *fs* (:class:`numpy.ndarray` or :class:`float`): the True anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Mean anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> fs = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 0.4
      >>> Ms = pk.f2m_v(fs, ecc)
      >>> np.shape(Ms)
      (100,)
)";
}

std::string m2f_v_doc()
{
    return R"(m2f_v(Ms, eccs)
    
    Converts from Mean to True anomaly (vectorized version). Requires ecc < 1.

    Args:
      *Ms* (:class:`numpy.ndarray` or :class:`float`): the Mean anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Ms = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 0.4
      >>> fs = pk.m2f_v(Ms, ecc)
      >>> np.shape(fs)
      (100,)
)";
}

std::string h2n_v_doc()
{
    return R"(h2n_v(Hs, eccs)
    
    Converts from Hyperbolic to Hyperbolic Mean anomaly (vectorized version). Requires ecc > 1.

    Args:
      *Hs* (:class:`numpy.ndarray` or :class:`float`): the Hyperbolic anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Hyperbolic Mean anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Hs = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 4.5
      >>> Ns = pk.h2n_v(Hs, ecc)
      >>> np.shape(Ns)
      (100,)
)";
}

std::string n2h_v_doc()
{
    return R"(n2h_v(Ns, eccs)
    
    Converts from Hyperbolic Mean to Hyperbolic anomaly (vectorized version). Requires ecc > 1.

    Args:
      *Ns* (:class:`numpy.ndarray` or :class:`float`): the Hyperbolic Mean anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Hyperbolic anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Ns = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 4.5
      >>> Hs = pk.n2h_v(Ns, ecc)
      >>> np.shape(Hs)
      (100,)
)";
}

std::string h2f_v_doc()
{
    return R"(h2f_v(Hs, eccs)
    
    Converts from Hyperbolic to True anomaly (vectorized version). Requires ecc > 1.

    Args:
      *Hs* (:class:`numpy.ndarray` or :class:`float`): the Hyperbolic anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Hs = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 4.5
      >>> fs = pk.h2f_v(Hs, ecc)
      >>> np.shape(fs)
      (100,)
)";
}

std::string f2h_v_doc()
{
    return R"(f2h_v(fs, eccs)
    
    Converts from True to Hyperbolic anomaly (vectorized version). Requires ecc > 1.

    Args:
      *fs* (:class:`numpy.ndarray` or :class:`float`): the True anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Hyperbolic anomaly

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> fs = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 5.7
      >>> Hs = pk.n2h_v(fs, ecc)
      >>> np.shape(Hs)
      (100,)
)";
}

std::string f2n_v_doc()
{
    return R"(f2n_v(fs, eccs)
    
    Converts from True to Hyperbolic Mean anomaly (vectorized version). Requires ecc > 1.

    Args:
      *fs* (:class:`numpy.ndarray` or :class:`float`): the True anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Hyperbolic Mean anomaly

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> fs = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 5.7
      >>> Ns = pk.n2f_v(fs, ecc)
      >>> np.shape(Ns)
      (100,)
)";
}

std::string n2f_v_doc()
{
    return R"(n2f_v(Ns, eccs)
    
    Converts from Hyperbolic Mean to True anomaly (vectorized version). Requires ecc > 1.

    Args:
      *Ns* (:class:`numpy.ndarray` or :class:`float`): the Hyperbolic Mean anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the True anomaly

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Ns = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 13.45
      >>> fs = pk.n2f_v(Ns, ecc)
      >>> np.shape(fs)
      (100,)
)";
}

std::string zeta2f_v_doc()
{
    return R"(zeta2f_v(zetas, eccs)
    
    Converts from Gudermannian to True anomaly (vectorized version). Requires ecc > 1.

    See Battin: "An Introduction to the Mathematics and Methods of Astrodynamics" for a 
    definition of zeta and the treatment of the resulting equations.

    Args:
      *zetas* (:class:`numpy.ndarray` or :class:`float`): the Gudermannian (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the True anomaly

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> zetas = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 2.2
      >>> fs = pk.zeta2f_v(zetas, ecc)
      >>> np.shape(fs)
      (100,)
)";
}

std::string f2zeta_v_doc()
{
    return R"(f2zeta_v(fs, eccs)
    
    Converts from True anomaly to Gudermannian (vectorized version). Requires ecc > 1.

    Args:
      *fs* (:class:`numpy.ndarray` or :class:`float`): the True anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Gudermannian 

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> fs = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 10.2
      >>> zetas = pk.f2zeta_v(fs, ecc)
      >>> np.shape(zetas)
      (100,)
)";
}

std::string epoch_from_float_doc()
{
    return R"(__init__(when: float, julian_type = MJD2000)
    
    Constructs an epoch from a Julian Date.

    Args:
      *when* (:class:`float`): the Julian Date (days since reference)

      *julian_type* (:class:`~pk.epoch.julian_type`): one of MJD2000, JD or MJD

    Examples:
      >>> import pykep as pk
      >>> pk.epoch(12.3, pk.epoch.julian_type.MJD2000)
      2000-01-13T07:12:00.000000
)";
}

std::string epoch_from_datetime_doc()
{
    return R"(**Alternative Constructor:**
    **__init__(** *when: datetime.datetime* **)**
    
    Constructs an epoch from a datetime object.

    Args:
      *when* (:class:`datetime.datetime`): a date

    Examples:
      >>> import pykep as pk
      >>> from datetime import datetime
      >>> pk.epoch(datetime(year=2000, month=1, day=13))
      2000-01-13T00:00:00.000000
)";
}

std::string epoch_from_string_doc()
{
    return R"(**Alternative Constructor:**
    **__init__(** *when: str*, *string_format = pk.epoch.string_format.ISO* **)**
    
    Constructs an epoch from a string.

    Args:
      *when* (:class:`str`): a date

      *string_format* (:class`~pykep.epoch.string_format`): string format.

    Examples:
      >>> import pykep as pk
      >>> pk.epoch("2000-01-14T00:00:00.000001")
      2000-01-14T00:00:00.000001
)";
}

std::string planet_docstring()
{
    return R"(__init__(udpla)

Planet class.

This type-erasing class represents a generic object moving in space. 
Basically anything which can be defined by its position and velocity in some reference frame.

In order to define a planet in pykep, the user must first define a class
whose methods describe the properties of the planet and allow to compute
its ephemerides (position and velocity), possibly its osculating elements, etc.. 
In pykep, we refer to such a class as a **user-defined planet**, or UDPLA for short. 
Once defined and instantiated, a UDPLA can then be used to construct an instance
of this class, the :class:`~pykep.planet`.

Every UDPLA must implement at least the following method:

.. code-block::

   def eph(self, epoch):
     ...

The ``eph()`` method is expected to return the Cartesian position and velocity at epoch
in some chosen reference frame.

The ``eph()`` method of the UDPLA will then be accessible from the corresponding
:func:`pykep.planet.eph()` method (see its documentation for information on how the method should be implemented
in the UDPLA and other details).

The mandatory method above allow to define a simple planet, which, in a minimal case,
could actually be also just a fixed point in space, for example if its ``eph()`` method
returns a constant position and zero velocity. 

In order to consider more complex cases, the UDPLA may implement one or more of the following methods:

.. code-block::

   def eph_v(self, mjd2000s):
     ...
   def get_mu_central_body(self):
     ...
   def get_mu_self(self):
     ...
   def get_radius(self):
     ...
   def get_safe_radius(self):
     ...
   def period(self, mjd2000):
     ...
   def elements(self, mjd2000, elements_type):
     ...
   def get_name(self):
     ...
   def get_extra_info(self):
     ...

See the documentation of the corresponding methods in this class for details on how the optional
methods in the UDPLA should be implemented and on how they are used by :class:`~pykep.planet`.


Args:
    *udpla*: a user-defined planet, either C++ or Python

Raises:
    *NotImplementedError*: if *udpla* does not implement the mandatory methods detailed above
    unspecified: any exception thrown by methods of the UDP invoked during construction,
    the deep copy of the UDP, the constructor of the underlying C++ class,
    failures at the intersection between C++ and Python (e.g., type conversion errors, mismatched function
    signatures, etc.)

)";
}

std::string planet_get_extra_info_docstring()
{
    return R"(get_extra_info()

PLanet's extra info.

If the UDPLA provides a ``get_extra_info()`` method, then this method will return the output of its ``get_extra_info()``
method. Otherwise, an empty string will be returned. 

The string representation of a :class:`~pykep.planet` contains the output of a call to this method.

Returns:
  :class:`str`: extra info about the UDPLA

Raises:
  unspecified: any exception thrown by the ``get_extra_info()`` method of the UDPLA

)";
}

std::string planet_get_name_docstring()
{
    return R"(get_name()

Planet's name.

If the UDPLA provides a ``get_name()`` method, then this method will return the output of its ``get_name()`` method.
Otherwise, an implementation-defined name based on the type of the UDPLA will be returned.

The string representation of a :class:`~pykep.planet` contains the output of a call to this method.

Returns:
    :class:`str`: the problem's name

)";
}

std::string planet_get_mu_central_body_docstring()
{
    return R"(get_mu_central_body()

The gravitational parameter in SI units (m^3/sec^2) of a main body of attraction.

If the UDPLA provides a ``get_mu_central_body()`` method, then this method will return the output of its ``get_mu_central_body()`` method.
Otherwise, -1 will be returned.


Returns:
    :class:`float`: the central body gravitational parameter.

)";
}

std::string planet_get_mu_self_docstring()
{
    return R"(get_mu_self()

The gravitational parameter in SI units (m^3/sec^2) of the planet.

If the UDPLA provides a ``get_mu_self()`` method, then this method will return the output of its ``get_mu_self()`` method.
Otherwise, -1 will be returned.


Returns:
    :class:`float`: the planet's gravitational parameter.

)";
}

std::string planet_get_radius_docstring()
{
    return R"(get_radius()

An average radius in SI units (m^3/sec^2) of the planet.

If the UDPLA provides a ``get_radius()`` method, then this method will return the output of its ``get_radius()`` method.
Otherwise, -1 will be returned.


Returns:
    :class:`float`: the planet's average radius.

)";
}

std::string planet_get_safe_radius_docstring()
{
    return R"(get_safe_radius()

The safe radius in SI units (m^3/sec^2) of the planet. This is mainly for use in planetary fly-manouvres as to avoid
the planet atmosphere or circumvent its radiation environment.

If the UDPLA provides a ``get_safe_radius()`` method, then this method will return the output of its ``get_safe_radius()`` method.
Otherwise, -1 will be returned.


Returns:
    :class:`float`: the planet's safe radius.

)";
}

std::string planet_eph_docstring()
{
    return R"(eph(when = 0.)

The planet ephemerides, i.e. its position and velocity.

In order to be able to construct a :class:`~pykep.planet` object, the user must provide his own UDPLA (User-Defined-Planet).
This is a class that must implement the method: 

.. code-block::

   def eph(self, mjd2000: float):
      ...
      return [[float, float, float], [float, float, float]]

.. note::
   In the udpla, the signature for eph demands a float as epoch (mjd2000). The planet, instead, constructed from the same udpla, will also allow :class:`~pykep.epoch`.
   

Args:
    *when* (:class:`float` or :class:`~pykep.epoch`): the epoch at which compute the period. When a :class:`float` is passed mjd2000 is assumed.

Returns:
    :class:`list` [:class:`list`, :class:`list`]: r and v, that is the final position and velocity after the propagation.

)";
}

std::string planet_eph_v_docstring()
{
    return R"(eph_v(mjd2000s)

The planet ephemerides, i.e. position and velocity (vectorized version over many epochs).

This method is the vectorized version of its companion :func:`~pykep.planet.eph` and, in its default implementation, it just
calls it in a loop. This behaviour can be changed by the user (for efficiency purposes) who can provide a more efficient version
in his UDPLA by coding a method having the signature: 


.. code-block::

   def eph_v(self, mjd2000s):
      ...
      return np.array((len(mjd2000s), 6))

see, for example, the python implementation of the UDPLAS :class:`~pykep.udpla.tle` and :class:`~pykep.udpla.spice`.

Args:
    *mjd2000s* (:class:`numpy.ndarray` or :class:`list`): the Modified Julian Dates at which to compute the ephemerides.

Returns:
    :class:`list` [:class:`list`, :class:`list`]: r and v, that is the final position and velocity after the propagation.

)";
}

std::string planet_period_docstring()
{
    return R"(period(when = 0.)

The period of the planet in seconds.

If the UDPLA provides a ``period(float)`` method, ``planet.period`` will call it.
Otherwise, if the UDPLA provides a ``get_mu_self()`` method ``planet.period`` will return the period as computed by the
equation:

.. math::
   T = 2 \pi \sqrt{\frac{a^3}{\mu}}

Else, -1 will be returned.

Args:
    *when* (:class:`float` or :class:`~pykep.epoch`): the epoch at which compute the period. When a :class:`float` is passed mjd2000 is assumed.

Returns:
    :class:`float`: the planet's period.

)";
}

std::string planet_elements_docstring()
{
    return R"(elements(when = 0., el_type = KEP_F)

The period of the planet in seconds.

If the UDPLA provides a ``elements(float, pk.el_type)`` method, then ``planet.elements`` will call it.
Otherwise, if the UDPLA provides a ``get_mu_self()`` method ``planet.elements`` will return the elements as computed by the
:func:`pykep.ic2par`. Otherwise, -1 will be returned.

Args:
    *when* (:class:`float` or :class:`~pykep.epoch`): the epoch at which compute the period. When a :class:`float` is passed mjd2000 is assumed.

    *el_type* (:class:`~pykep.el_type`): the elements type.

Returns:
    :class:`list`: the planet's elements at epoch.

)";
}

std::string udpla_keplerian_from_posvel_docstring()
{
    return R"(**Alternative Constructor:**
    __init__(ep, posvel, mu_central_body, name = "unkown", added_params = [-1,-1,-1])

Constructs a Keplerian udpla from its position and velocity at epoch.

Args:
    *ep* (:class:`~pykep.epoch`): the epoch at which the orbital elements are provided.

    *posvel* (:class:`list` [:class:`list`, :class:`list`]): the body position and velocty.

    *mu_central_body* (:class:`float`): the gravitational parameter of the main attracting body.

    *name* (:class:`str`): the name of the orbiting body.

    *added_params* (:class:`list`): the body gravitational parameter, its radius and its safe radius. (if -1 they are assumed unkown)

Examples:
    >>> import pykep as pk
    >>> r = [1, 0, 0]
    >>> v = [0, 1, 0]
    >>> ep = pk.epoch("2025-03-22")
    >>> udpla = pk.udpla.keplerian(ep = ep, posvel = [r, v], mu_central_body =1, name = "my_pla")
    >>> pla = pk.planet(udpla)
)";
}

std::string udpla_keplerian_from_elem_docstring()
{
    return R"(__init__(ep, elem, mu_central_body, name = "unkown", added_params = [-1,-1,-1], elem_type = KEP_F)

Constructs a Keplerian udpla from its orbital elements at epoch.

Args:
    *ep* (:class:`~pykep.epoch`): the epoch at which the orbital elements are provided.

    *elem* (:class:`list` or :class:`numpy.ndarray`): the orbital elements. by default.

    *mu_central_body* (:class:`float`): the gravitational parameter of the main attracting body.

    *name* (:class:`str`): the name of the orbiting body.

    *added_params* (:class:`list`): the body gravitational parameter, its radius and its safe radius. (if -1 they are assumed unkown)

    *el_type* (:class:`~pykep.el_type`): the elements type. Deafulets to osculating Keplrian (a ,e ,i, W, w, f) with true anomaly.

Examples:
    >>> import pykep as pk
    >>> elem = [1, 0, 0, 0, 0, 0]
    >>> ep = pk.epoch("2025-03-22")
    >>> udpla = pk.udpla.keplerian(ep = ep, elem = elem, mu_central_body =1, name = "my_pla")
    >>> pla = pk.planet(udpla)
)";
}

std::string udpla_jpl_lp_docstring()
{
    return R"(__init__(name = "earth")

Constructs a solar system planet with ephemerides computed using a low-precision (non Keplerian)
model from JPL (https://ssd.jpl.nasa.gov/planets/approx_pos.html).

Args:
    *name* (:class:`str`): the name of the solar system planet.

Examples:
    >>> import pykep as pk
    >>> udpla = pk.udpla.jpl_lp(name="mercury")
    >>> pla = pk.planet(udpla)
)";
}

std::string udpla_vsop2013_docstring()
{
    return R"(__init__(body = "mercury", thresh = 1e-5)

Constructs a solar system planet with ephemerides computed using the VSOP2013 analytical
theory (https://en.wikipedia.org/wiki/VSOP_model).

Args:
    *body* (:class:`str`): the name of the solar system planet.

    *thresh* (:class:`float`): the truncation threshold for the theory's coefficients.

Examples:
    >>> import pykep as pk
    >>> udpla = pk.udpla.vsop2013(body="venus")
    >>> pla = pk.planet(udpla)
)";
}

std::string lambert_problem_docstring()
{
    return R"(__init__(r0 = [1,0,0], r1 = [0,1,0], tof = pi/2, mu = 1., cw = False, max_revs = 0)

      Args:
          *r0* (1D array-like): Cartesian components of the first position vector [xs, ys, zs]. Defaults to [1,0,0].

          *r1* (1D array-like): Cartesian components of the second position vector [xf, yf, zf]. Defaults tot [0,1,0].

          *tof* (:class:`float`): time of flight. Defaults to :math:`\frac{\pi}{2}`.

          *mu* (:class:`float`): gravitational parameter. Defaults to 1.

          *cw* (:class:`bool`): True for retrograde motion (clockwise). Defaults to False.

          *max_revs* (:class:`float`): Maximum number of multiple revolutions to be computed. Defaults to 0.

      .. note::

        Units need to be consistent. The multirev Lambert's problem will be solved upon construction
        and its solution stored in data members.

      Examples:
        >>> import pykep as pk
        >>> import numpy as np
        >>> r0 = [1,0,0]
        >>> r1 = [0,1,0]
        >>> tof = np.pi/2
        >>> mu = 1.
        >>> lp = pk.lambert_problem(r0, r1, tof, mu)
        >>> lp.v0[0]
        [-4.1028493158958256e-16, 1.0000000000000002, 0.0]
)";
}

std::string stark_problem_docstring()
{
    return R"(__init__(mu = 1., veff = 1., tol = 1e-16)

Class representing the Stark problem. That is the initial value problem
of a fixed inertial thrust dynamics described by the equations:

.. math::
   \left\{
   \begin{array}{l}
       \dot{\mathbf r} = \mathbf v \\
       \dot{\mathbf v} = -\frac{\mu}{r^3} \mathbf r + \frac{\mathbf T}{m} \\
       \dot m = - \frac{|\mathbf T|}{I_{sp} g_0}
   \end{array}\right.

Args:
    *mu* (:class:`float`): central body gravitational parameter. Defaults to 1.

    *veff* (:class:`float`): propulsion system effective velocity (Isp g0). Defaults to 1.

    *tol* (:class:`float`): tolerance of the Taylor adaptive integration. Defaults to 1e-16.

.. note::

  Units need to be consistent upon construction as well as later on when calling the propagate methods.

Examples:
  >>> import pykep as pk
  >>> import numpy as np
  >>> mu = pk.MU_SUN
  >>> veff = 3000. * pk.G0
  >>> tol = 1e-14
  >>> sp = pk.stark_problem(mu, veff, tol)
  >>> sp.propagate(rvm_state = [1., 0., 0., 0., 1., 0., 1], thrust = [0., 0., 1e-8], tof = 7.32)
  [0.5089647068650076, 0.8607873878989034, 0.0, -0.8607873878989032, 0.5089647068650074, 0.0, 1.0]
)";
}

std::string stark_problem_propagate_docstring()
{
    return R"(propagate(rvm_state, thrust, tof)

Stark problem numerical propagation. 
The propagation will be singular for vanishing masses (infinite acceleration) and raise an exception.

Args:
    *rvm_state* (:class:`list` (7,)): position, velocity and mass flattened into a 7D list. 

    *thrust* (:class:`list` (3,)): thrust flattened into a 3D list. 

    *tof* (:class:`float`): time of flight.

Returns:
    :class:`list` (7,) : position, velocity and mass after the numerical propagation, flattened into a 7D list. 

.. note::

  Units need to be consistent with the ones used upon constructing the instance.

Examples:
  >>> import pykep as pk
  >>> import numpy as np
  >>> mu = pk.MU_SUN
  >>> veff = 3000. * pk.G0
  >>> tol = 1e-14
  >>> sp = pk.stark_problem(mu, veff, tol)
  >>> sp.propagate(rvm_state = [1., 0., 0., 0., 1., 0., 1], thrust = [0., 0., 1e-8], 7.32)
  [0.5089647068650076, 0.8607873878989034, 0.0, -0.8607873878989032, 0.5089647068650074, 0.0, 1.0]
)";
}

std::string stark_problem_propagate_var_docstring()
{
    return R"(propagate_var(rvm_state, thrust, tof)

Stark problem numerical propagation via variational equations. 
The propagation will be singular for vanishing masses (infinite acceleration) and raise an exception.

It also computes the system State Transition Matrix:

.. math::
    \mathbf M = \frac{d\mathbf x_f}{d\mathbf x_0}

as well as the gradients of the final states with respect to the thrust direction.

.. math::
    \mathbf U = = \frac{d\mathbf x_f}{d\mathbf u}

Args:
    *rvm_state* (:class:`list` (7,)): position, velocity and mass flattened into a 7D list. 

    *thrust* (:class:`list` (3,)): thrust flattened into a 3D list. 

    *tof* (:class:`float`): time of flight.

Returns:
    :class:`tuple` (:class:`list` (7,), :class:`numpy.ndarray` (7,7), :class:`numpy.ndarray` (7,3)): position, velocity and mass after ropagation flattened into a 7D list, state transition matrix 
    and 

.. note::

  Units need to be consistent upon construction and later on when calling the propagate methods.

Examples:
  >>> import pykep as pk
  >>> import numpy as np
  >>> mu = pk.MU_SUN
  >>> veff = 3000. * pk.G0
  >>> tol = 1e-14
  >>> sp = pk.stark_problem(mu, veff, tol)
  >>> sp.propagate(state = [1., 0., 0., 0., 1., 0., 1], thrust = [0., 0., 1e-8], 7.32)
)";
}

std::string ta_stark_docstring()
{return "";}
std::string ta_stark_var_docstring()
{return "";}
std::string ta_stark_dyn_docstring()
{return "";}


std::string propagate_lagrangian_docstring()
{
    return R"(propagate_lagrangian(rv = [[1,0,0], [0,1,0]], tof = pi/2, mu = 1, stm = False)

    Propagates (Keplerian) the state for an assigned time and computes the State Transition Matrix (if requested) using the Lagrangian coefficients.

    Args:
          *rv* (2D array-like): Cartesian components of the initial position vector and velocity [[x0, y0, z0], [v0, vy0, vz0]]. Defaults to [[1,0,0], [0,1,0]].

          *tof* (:class:`float`): time of flight. Defaults to :math:`\frac{\pi}{2}`.

          *mu* (:class:`float`): gravitational parameter. Defaults to 1.

          *stm* (:class:`bool`): requests the computations of the State Transition Matrix

    Returns:
          :class:`tuple` (:class:`list`, :class:`list`): r and v, that is the final position and velocity after the propagation. (if *stm* is False)
          :class:`tuple` (:class:`list` [:class:`list`, :class:`list`], :class:`numpy.ndarray` (6,6)): [r,v] and the STM. (if *stm* is True)

    Examples:
        >>> import pykep as pk
        >>> import numpy as np
        >>> r0 = [1,0,0]
        >>> v0 = [0,1,0]
        >>> tof = pi/2
        >>> mu = 1
        >>> [r1,v1], stm = pk.propagate_lagrangian(rv=[r0,v0], tof = tof, mu = mu, stm = True)
        >>> [r1,v1] = pk.propagate_lagrangian(rv=[r0,v0], tof = tof, mu = mu, stm = False)

)";
}

std::string propagate_lagrangian_v_docstring()
{
    return R"(propagate_lagrangian_v(rv = [[1,0,0], [0,1,0]], tofs = [pi/2], mu = 1, stm = False)

    This is the vectorized version of :func:`pykep.propagate_lagrangian`. Vectorization allows to compute many
    different time of flights at once. Note that this is not necessarily more efficient than calling
    :func:`pykep.propagate_lagrangian` in a loop, since there is no parallelization nor SIMD magic implemented atm. 
    Nevertheless we offer this interface for cenvenience as it may allow more compact code. 

    Args:
          *rv* (2D array-like): Cartesian components of the initial position vector and velocity [[x0, y0, z0], [v0, vy0, vz0]]. Defaults to [[1,0,0], [0,1,0]].

          *tof* (1D array-like): time of flight. Defaults to [:math:`\frac{\pi}{2}`].

          *mu* (:class:`float`): gravitational parameter. Defaults to 1.

          *stm* (:class:`bool`): requests the computations of the State Transition Matrix

    Returns:
          :class:`list` [:class:`tuple` ( :class:`list` , :class:`list` ) ]: For each time of flight: [r,v], that is the final position
          and velocity after the propagation and the flattened stm (if requested).

    Examples:
        >>> import pykep as pk
        >>> import numpy as np
        >>> r0 = [1,0,0]
        >>> v0 = [0,1,0]
        >>> tofs = [pi/2, pi, 3/4pi]
        >>> mu = 1
        >>> res = pk.propagate_lagrangian_v(rv = [r0, v0], tofs = tofs, mu = mu, stm = True)
        >>> rs = [it[0][0] for it in res]
        >>> vs = [it[0][1] for it in res]
        >>> stms = [it[1] for it in res]
)";
}

std::string leg_sf_docstring()
{
    return R"(__init__(rvs = [[1,0,0], [0,1,0]], ms = 1., throttles = [0,0,0,0,0,0], rvf = [[0,1,0], [-1,0,0]], mf = 1., tof = pi/2, max_thrust = 1., isp = 1., mu=1., cut = 0.5)

      This class represents an interplanetary low-thrust transfer between a starting and a final point in the augmented state-space :math:`[\mathbf r, \mathbf v, m]`.
      The low-thrust transfer is described by a sequence of equally spaced impulses as described in:

      Sims, J., Finlayson, P., Rinderle, E., Vavrina, M. and Kowalkowski, T., 2006, August. Implementation of a low-thrust trajectory optimization algorithm for preliminary design. 
      In AIAA/AAS Astrodynamics specialist conference and exhibit (p. 6746).

      The low-thrust transfer will be feasible is the state mismatch equality constraints and the throttle mismatch inequality constraints are satisfied.

      Args:
          *rvs* (2D array-like): Cartesian components of the initial position vector and velocity [[xs, ys, zs], [vxs, vys, vzs]]. Defaults to [[1,0,0], [0,1,0]].

          *ms* (:class:`float`): initial mass. Defaults to 1.

          *throttles* (1D array-like): the Cartesan components of the throttle history [ux1, uy1, uz1, ux2, uy2, uz2, .....]. Defaults to a ballistic, two segments profile [0,0,0,0,0,0].

          *rvf* (2D array-like): Cartesian components of the final position vector and velocity [[xf, yf, zf], [vxf, vyf, vzf]]. Defaults to [[0,1,0], [-1,0,0]].

          *mf* (:class:`float`): final mass. Defaults to 1.

          *tof* (:class:`float`): time of flight. Defaults to :math:`\frac{\pi}{2}`.

          *max_thrust* (:class:`float`): maximum level for the spacecraft thrust. Defaults to 1.

          *isp* (:class:`float`): specific impulse of the propulasion system. Defaults to 1.

          *mu* (:class:`float`): gravitational parameter. Defaults to 1.

          *cut* (:class:`float`): the leg cut, in [0,1]. It determines the number of forward and backward segments. Defaults to 0.5.

      .. note::

        Units need to be consistent. 

      Examples:
        >>> import pykep as pk
        >>> import numpy as np
        >>> sf = pk.leg.sims_flanagan()
)";
}
std::string leg_sf_rvs_docstring()
{
    return "The initial position vector and velocity: [[xs, ys, zs], [vxs, vys, vzs]].";
};
std::string leg_sf_ms_docstring()
{
    return "Initial mass.";
};
std::string leg_sf_throttles_docstring()
{
    return "The Cartesan components of the throttle history [ux1, uy1, uz1, ux2, uy2, uz2, .....].";
};
std::string leg_sf_rvf_docstring()
{
    return "The final position vector and velocity: [[xs, ys, zs], [vxs, vys, vzs]].";
};
std::string leg_sf_mf_docstring()
{
    return "Final mass.";
};
std::string leg_sf_tof_docstring()
{
    return "Time of flight.";
};
std::string leg_sf_max_thrust_docstring()
{
    return "Maximum spacecraft thruet.";
};
std::string leg_sf_isp_docstring()
{
    return "Specific impulse of the propulasion system";
};
std::string leg_sf_mu_docstring()
{
    return "Central body gravitational parameter.";
};
std::string leg_sf_cut_docstring()
{
    return "The leg cut: it determines the number of forward and backward segments.";
};
std::string leg_sf_nseg_docstring()
{
    return "The total number of segments";
};
std::string leg_sf_nseg_bck_docstring()
{
    return "The total number of backward segments";
};
std::string leg_sf_nseg_fwd_docstring()
{
    return "The total number of forward segments";
};
std::string leg_sf_mc_docstring()
{
    return R"(compute_mismatch_constraints()

      In the Sims-Flanagan trajectory leg model, a forward propagation is performed from the starting state as well as a backward from the final state.
      The state values thus computed need to match in some middle control point. This is typically imposed as 7 independent constraints called mismatch-constraints
      computed by this method. 

      Returns:
          :class:`list` [:class:`float`]: The seven mismatch constraints in the same units used to construct the leg.
      
      Examples:
        >>> import pykep as pk
        >>> import numpy as np
        >>> sf = pk.leg.sims_flanagan()
        >>> sf.compute_mismatch_constraints()
)";
};
std::string leg_sf_tc_docstring()
{
    return R"(compute_throttle_constraints()

      In the Sims-Flanagan trajectory leg model implemented in pykep, we introduce the concept of throttles. Each throttle is defined by three numbers
      :math:`[u_x, u_y, u_z] \in [0,1]` indicating that a certain component of the thrust vector has reached a fraction of its maximum allowed value. 
      As a consequence, along the segment along which the throttle is applied, the constraint  :math:`u_x ^2 + u_y ^2 + u_z^2 = 1`, called a throttle constraint,
      has to be met. 

      Returns:
          :class:`list` [:class:`float`]: The throttle constraints.
      
      Examples:
        >>> import pykep as pk
        >>> import numpy as np
        >>> sf = pk.leg.sims_flanagan()
        >>  sf.throttles = [0.8]*3
        >>> sf.compute_throttle_constraints()
)";
};
std::string leg_sf_mc_grad_docstring()
{
    return R"(compute_mc_grad()

Computes the gradients of the mismatch constraints. Indicating the initial augmented state with :math:`\mathbf x_s = [\mathbf r_s, \mathbf v_s, m_s]`, the
final augmented state with :math:`\mathbf x_f = [\mathbf r_f, \mathbf v_f, m_f]`, the total time of flight with :math:`T` and the introducing the augmented throttle vector
:math:`\mathbf u = [u_{x0}, u_{y0}, u_{z0}, u_{x1}, u_{y1}, u_{z1} ..., T]` (note the time of flight at the end), this method computes the following gradients:

.. math::
  \frac{\partial \mathbf {mc}}{\partial \mathbf x_s}

.. math::
  \frac{\partial \mathbf {mc}}{\partial \mathbf x_f}

.. math::
  \frac{\partial \mathbf {mc}}{\partial \mathbf u}

Returns:
    :class:`tuple` [:class:`numpy.ndarray`, :class:`numpy.ndarray`, :class:`numpy.ndarray`]: The three gradients. sizes will be (7,7), (7,7) and (7,nseg*3)

Examples:
  >>> import pykep as pk
  >>> import numpy as np
  >>> sf = pk.leg.sims_flanagan()
  >>  sf.throttles = [0.8]*3
  >>> sf.compute_mc_grad()
)";
};

std::string leg_sf_tc_grad_docstring()
{
    return R"(compute_tc_grad()

Computes the gradients of the throttles constraints. Indicating the total time of flight with :math:`T` and introducing the augmented throttle vector
:math:`\mathbf u = [u_{x0}, u_{y0}, u_{z0}, u_{x1}, u_{y1}, u_{z1} ..., T]` (note the time of flight at the end), this method computes the following gradient:

.. math::
  \frac{\partial \mathbf {tc}}{\partial \mathbf u}

Returns:
    :class:`tuple` [:class:`numpy.ndarray`]: The gradient. Size will be (nseg,nseg*3).

Examples:
  >>> import pykep as pk
  >>> import numpy as np
  >>> sf = pk.leg.sims_flanagan()
  >>  sf.throttles = [0.8]*3
  >>> sf.compute_tc_grad()
)";
};

} // namespace pykep