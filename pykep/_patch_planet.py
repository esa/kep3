## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)## 
## This file is part of the kep3 library.## 
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

from .core import planet

def _planet_extract(self, t):
    """Extract the user-defined planet.

    This method allows to extract a reference to the user-defined planet (UDPLA) stored within this
    :class:`~pykep.planet` instance. The behaviour of this function depends on the value
    of *t* (which must be a :class:`type`) and on the type of the internal UDPLA:

    * if the type of the UDPLA is *t*, then a reference to the UDP will be returned
      (this mirrors the behaviour of the corresponding C++ method
      :cpp:func:`kep3::planet::extract()`),
    * if *t* is :class:`object` and the UDP is a Python object (as opposed to an
      :ref:`exposed C++ planet `), then a reference to the
      UDPLA will be returned (this allows to extract a Python UDPLA without knowing its type),
    * otherwise, :data:`None` will be returned.

    Args:
        t (:class:`type`): the type of the user-defined planet to extract

    Returns:
        a reference to the internal user-defined planet, or :data:`None` if the extraction fails

    Raises:
        TypeError: if *t* is not a :class:`type`

    Examples:
        >>> import pygmo as pg
        >>> p1 = pg.problem(pg.rosenbrock())
        >>> p1.extract(pg.rosenbrock) # doctest: +SKIP
        <pygmo.core.rosenbrock at 0x7f56b870fd50>
        >>> p1.extract(pg.ackley) is None
        True
        >>> class prob:
        ...     def fitness(self, x):
        ...         return [x[0]]
        ...     def get_bounds(self):
        ...         return ([0],[1])
        >>> p2 = pg.problem(prob())
        >>> p2.extract(object) # doctest: +SKIP
        <__main__.prob at 0x7f56a66b6588>
        >>> p2.extract(prob) # doctest: +SKIP
        <__main__.prob at 0x7f56a66b6588>
        >>> p2.extract(pg.rosenbrock) is None
        True

    """
    if not isinstance(t, type):
        raise TypeError("the 't' parameter must be a type")
    if hasattr(t, "_pykep_cpp_udpla"):
        return self._cpp_extract(t())
    return self._py_extract(t)


def _planet_is(self, t):
    """Check the type of the user-defined problem.

    This method returns :data:`False` if :func:`extract(t) <pygmo.problem.extract>` returns
    :data:`None`, and :data:`True` otherwise.

    Args:
        t (:class:`type`): the type that will be compared to the type of the UDP

    Returns:
        bool: whether the UDP is of type *t* or not

    Raises:
        unspecified: any exception thrown by :func:`~pygmo.problem.extract()`

    """
    return not self.extract(t) is None


# Do the actual patching.
setattr(planet, "extract", _planet_extract)
setattr(planet, "is_", _planet_is)