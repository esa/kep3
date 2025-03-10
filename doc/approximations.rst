.. _approximations:

Various approximations for orbital transfers
############################################

Computing the exact optimal transfer between orbits is often more expensive than allowed in
preliminary phases of the mission design. For this reason, `pykep` offers a number of approximations
that can be used as surrogates of the actual complete computations. Many of these approximations have
been used and developed during the various edition of the [GTOC competition](https://sophia.estec.esa.int/gtoc_portal/)
and can be found in papers :cite:p:`approximations`, :cite:p:`gtoc12`.

.. currentmodule:: pykep

.. autofunction:: mima

.. autofunction:: mima2