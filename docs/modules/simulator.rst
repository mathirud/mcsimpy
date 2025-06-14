Simulator
=========

Base Vessel
-----------

The `Vessel` class works as an abstract base class for all vessel models. All subclasses of  `Vessel` inherits
the mehtods of `Vessel`.

.. autoclass:: mcsimpy.simulator.vessel.Vessel
    :members:

Gunnerus
--------

Simulation models for R/V Gunnerus.


RVG 3 DOF Manuevering Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: mcsimpy.simulator.gunnerus.GunnerusManeuvering3DoF
    :members:
    :inherited-members:


RVG 6 DOF DP Model
^^^^^^^^^^^^^^^^^^

.. autoclass:: mcsimpy.simulator.gunnerus.RVG_DP_6DOF
    :members:
    :inherited-members:

CSAD
----

Simulation models for C/S Arctic Drillship.

CSAD 3DOF Maneuvering Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: mcsimpy.simulator.csad.CSADMan3DOF
    :members:
    :inherited-members:


CSAD 6DOF DP Model
^^^^^^^^^^^^^^^^^^

.. autoclass:: mcsimpy.simulator.csad.CSAD_DP_6DOF
    :members:
    :inherited-members: