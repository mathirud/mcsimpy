Control
=======

The ``mcsimpy.control`` package contains a set of DP and Maneuvering controllers.

Basic Controllers
-----------------

The module contains simple PD and PID controllers.

PD Controller
^^^^^^^^^^^^^

.. autoclass:: mcsimpy.control.basic.PD
    :members:


PID Controller
^^^^^^^^^^^^^^

.. autoclass:: mcsimpy.control.basic.PID
    :members:

PI Controller
^^^^^^^^^^^^^

.. autoclass:: mcsimpy.control.basic.PI
    :members:


Direct Bias Compensation Controller
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: mcsimpy.control.basic.DirectBiasCompensationController
    :members:


Adaptive Controller
-------------------

.. autoclass:: mcsimpy.control.adaptiveFS.AdaptiveFSController
    :members:

Backstepping Controllers
------------------------

In the ``mcsimpy.control.backstepping`` module there are controllers based on backstepping
theory. Currently, only a backstepping controller for maneuvering purposes has been implemented.

.. autoclass:: mcsimpy.control.backstepping.BacksteppingController
    :members:
    :no-undoc-members: