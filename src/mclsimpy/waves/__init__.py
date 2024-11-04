"""
Waves (:mod:`mclsimpy.waves`)
================================

Waves package for simulations.

Contents
--------

 - wave_loads.py: Wave load calculations
 - wave_spectra.py: Wave spectra module

Examples
--------

>>> from mclsimpy.waves import JONSWAP, WaveLoad
"""

from mclsimpy.waves.wave_loads import WaveLoad, FluidMemory
from mclsimpy.waves.wave_spectra import JONSWAP, ModifiedPiersonMoskowitz