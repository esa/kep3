.. _gym:

Trajectory Optimization Gym
###########################

.. currentmodule:: pykep.trajopt.gym

MGA problems
************************

.. autoattribute:: pykep.trajopt.gym.cassini1

This is an MGA problem inspired to the Cassini spacecraft interplanetary transfer to Saturn. 
The objective of this mission is to reach Saturn and to be captured by its gravity into an orbit having pericenter radius :math:`r_p=108950` km, 
and eccentricity :math:`e=0.98`. The planetary fly-by sequence considered is E-VVEJ-S (as the one used by Cassini spacecraft). 
As objective function we use the total :math:`\Delta V` accumulated during the mission, including the launch :math:`\Delta V` and the various :math:`\Delta V` one 
needs to give at the planets and upon arrival to perform the final orbit injection. 