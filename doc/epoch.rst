.. _epoch:

Relatated to epochs
====================

The representation of an epoch, that is of a specific point in time, be it in the future or in the past, can be rather confusing. 
In `pykep` we opted to offer a dedicated class called epoch to offer a simple interface and, under the hoods, interfacing
seamlessly both to the c++ `std::chrono <https://en.cppreference.com/w/cpp/header/chrono>`_ 
library and to the python `datetime <https://docs.python.org/3/library/datetime.html>`_ module. 

.. note::
    In pykep the default Julian Date is the Modified Julian Date, defined as a `float` representing the number
    of days since the start of 2000-1-1. 


.. currentmodule:: pykep

-----------------------------------

.. autoclass:: pykep.epoch
   :members:
   :special-members: __init__

----------------------------------

.. autofunction:: pykep.utc_now


