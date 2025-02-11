.. _taylor_adaptive:

Taylor adaptive propagators
############################


Taylor adaptive integrators are offered in `pykep` wrapping some of the functionalities of
`Heyoka <https://bluescarni.github.io/heyoka.py/index.html>`_ :cite:p:`biscaniheyoka1` python package.
Their variational version is also offered (at order one) as to be able to produce stms and, where needed,
more. Higher order variational equations can also be obtained directly using the available dynamics and 
using `Heyoka <https://bluescarni.github.io/heyoka.py/index.html>`_ :cite:p:`biscaniheyoka1` syntax.

Some of the Taylor adaptive integrators are associated to OCPs (Optimal Control Problems) of
relevance to interplanetary flight. In particular they are birn when applying Pontryagin 
principle to dynamics of interests. In these cases, a number of auxiliary functiond describing
the control, the Hamiltnian, the switching function, etc., are also provided.

--------------------------------------------------------

Stark
*************************
 
.. currentmodule:: pykep.ta

.. autofunction:: get_stark

.. autofunction:: get_stark_var

.. autofunction:: stark_dyn 


Circular Restricted Three Body Problem
*****************************************
 
.. currentmodule:: pykep.ta

.. autofunction:: get_cr3bp

.. autofunction:: get_cr3bp_var

.. autofunction:: cr3bp_dyn 

Low-thrust Pontryagin Cartesian Problem
*****************************************
 
.. currentmodule:: pykep.ta

.. autofunction:: get_pc

.. autofunction:: get_pc_var

.. autofunction:: pc_dyn
