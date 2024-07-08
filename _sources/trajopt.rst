.. _trajopt:

Trajectory Optimization
=======================

In `pykep` both direct and indirect optimization techniques are provided
to perform spacecraft trajectory optimization. Most direct techniques provided are
based on the "Sims" transcription :cite:p:`sims`, while indirect techniques are based on some
version of Pontryagin Maximum Principle.

Most of the classes in this Trajectory Optimization module are provided as
optimization problems (NLP) compatible to the pagmo software suite.

.. currentmodule:: pykep.trajopt

Direct
------
.. autoclass:: direct_point2point
    :members: pretty, plot

.. autoclass:: direct_pl2pl
    :members: pretty, plot


