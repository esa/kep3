.. _trajopt:

Trajectory Optimization
###########################

In `pykep` both direct and indirect optimization transcriptions as well as evolutionary encodings
are provided to perform spacecraft trajectory optimization. Most direct techniques provided are
based on the "Sims" transcription :cite:p:`sims`, while indirect techniques are based on some
version of Pontryagin Maximum Principle. The evolutionary encoding are mostly based on the work
performed at `ESA' Advanced Concepts Team <https://www.esa.int/gsp/ACT/>`_ (:cite:p:`izzo2010global`, :cite:p:`izzo2013search`).

Most of the classes in this Trajectory Optimization module are provided as
optimization problems (NLPs) in the for of User Defined Problems (UDP) compatible
with the `pygmo <https://esa.github.io/pygmo2/>`_ :cite:p:`pagmo` python package.

.. currentmodule:: pykep.trajopt

Direct
************************
When solving optimal control problems (OCPs), the unknown is a function (i.e. the control), 
hence one is dealing with an infinite dimensional search space.
A common approach to overcome this difficulty is to use direct methods. 
These are based on the idea of discretizing the problem introducing some time grid
and thus transforming the OCP into a NLP (Non Linear Programming) problem
to then solve it using numerical solvers available for the task.

-------------------------------------------------------

.. autoclass:: direct_point2point
    :members: pretty, plot

-------------------------------------------------------

.. autoclass:: direct_pl2pl
    :members: pretty, plot

-------------------------------------------------------

Indirect
******************
Indirect methods are based on the Pontryagin Maximum Principle (PMP), which provides necessary conditions
for optimality. These methods involve deriving the optimal control laws and
the corresponding state trajectories by solving a two point boundary value problem
(TPBVP) derived from the PMP.

-------------------------------------------------------

.. autoclass:: pontryagin_cartesian_mass
    :members: plot, plot_misc

.. autoclass:: pontryagin_cartesian_time
    :members: plot, plot_misc

-------------------------------------------------------


Evolutionary Encodings
************************
Some type of interplanetary trajectories can be *evolved*: shocking?. 
This AI approach to trajectory design resulted in two silver (:cite:p:`terrile2005evolutionary`, :cite:p:`petropoulos2018gtoc9`) 
and one gold (:cite:p:`izzo2013search`)
`Humies award <https://www.human-competitive.org>`_  for human-competitive
results that were produced by any form of genetic and evolutionary computation.
In `pykep` we offer some classes that help instantiating these type of transfers
amenable to evolutionary techniques.
 
-------------------------------------------------------

.. autoclass:: mga
    :members: pretty, plot, to_planet

.. autoclass:: mga_1dsm
    :members: pretty, plot, to_planet

.. autoclass:: pl2pl_N_impulses
    :members: pretty, plot, plot_primer_vector

Utilities
*********
In order to facilitate the use of the classes in this module, some utilities are provided.

-------------------------------------------------------

.. autofunction:: primer_vector

.. autofunction:: primer_vector_surrogate

.. autoclass:: _launchers
    :members: atlas501, atlas551, soyuzf, ariane5