"""
Observer (:mod:`mcsimpy.observer`)
=======================================

The package consist of different observers
for state estimation.

Contents
--------

 - ekf.py: Extended Kalman Filter
 - nonlinobs.py: Nonlinear observers

Description
-----------

The different observers are implemented as
classes, and are structured in different modules as
listed in `contents`.

Examples
--------

>>> from mcsimpy.observer import NonlinObs3dof

Alternative import

>>> from mcsimpy.observer.nonlinobs import NonlinObs3dof

Third option

>>> import mcsimpy.observer as obs
>>> estimator = obs.NonlinObs3dof(*args, **kwargs)
"""

from mcsimpy.observer.nonlinobs import NonlinObs3dof
from mcsimpy.observer.ekf import EKF