"""
Control (:mod:`mcsimpy.control`)
====================================

Controller package for DP and maneuvering
control.

Contents
--------

 - backstepping.py: Controllers based on backstepping theory.
 - basic.py: Basic PD and PID controller.

Examples
--------

>>> from mcsimpy.control import PD, PID
"""

from mcsimpy.control.backstepping import BacksteppingController
from mcsimpy.control.basic import PD, PID