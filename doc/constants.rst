.. _constants:

Global constants
=======================

In `pykep` the access a number of common constants are provided for convenience. The user can overwrite their values if needed.
These constants are not used in the `pykep` internals, they are only provided for convenience for the user to instantiate / use
the various `pykep` objects and functionalities.

.. list-table:: Pykep global constants
   :widths: 50 25 25 50
   :header-rows: 1

   * - Constant's Name
     - Symbol in pykep
     - Units
     - Value
   * - Astronomical Unit
     - pykep.AU
     - :math:`m` 
     - 149597870700.0
   * - Cavendish constant
     - pykep.CAVENDISH
     - :math:`\frac{N m^2}{kg^2}` 
     - 7.36687e-10
   * - Sun's gravitational parameter
     - pykep.MU_SUN
     - :math:`\frac{m^3}{sec^2}` 
     - 1.32712440018e+20
   * - Earth's gravitational parameter
     - pykep.MU_EARTH
     - :math:`\frac{m^3}{sec^2}` 
     - 398600441800000.0
   * - Earth's velocity
     - pykep.EARTH_VELOCITY
     - :math:`\frac{m}{sec}` 
     - 29784.691831696804
   * - Earth's radius
     - pykep.EARTH_RADIUS
     - :math:`\frac{m}{sec}` 
     - 6378137.0
   * - Earth's :math:`J_2`
     - pykep.EARTH_J2
     - --
     - 0.00108262668
   * - Seconds in one day
     - pykep.DAY2SEC
     - --
     - 86400.0
   * - Degrees in one radians
     - pykep.RAD2DEG
     - --
     - 57.29577951308232
