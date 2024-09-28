.. _trajopt:

Trajectory Optimization
###########################

In `pykep` both direct and indirect optimization transcriptions as well as evolutionary encodings
are provided to perform spacecraft trajectory optimization. Most direct techniques provided are
based on the "Sims" transcription :cite:p:`sims`, while indirect techniques are based on some
version of Pontryagin Maximum Principle. The evolutionary encoding are mostly based on the work
performed at `ESA' Advanced Concepts Team <https://www.esa.int/gsp/ACT/>`_ (:cite:p:`izzo2010global`, :cite:p:`izzo2013search`).

Most of the classes in this Trajectory Optimization module are provided as
optimization problems (NLP) compatible to the pagmo software suite.

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