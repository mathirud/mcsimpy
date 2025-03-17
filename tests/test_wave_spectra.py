#!/usr/bin/env python3

from mclsimpy.waves.wave_spectra import ModifiedPiersonMoskowitz, JONSWAP
from mclsimpy.waves.wave_spreading import MultiDirectional, MultiDirectionalAlt


import numpy as np

import sys
import os

cwd = os.getcwd()
sys.path.insert(0, cwd)


class TestSpectra:
    def test_realization_pm(self):
        # Test statistical output from time realization of wave MPM spectra.
        hs = 3
        tp = 12
        wp = 2 * np.pi / tp

        freqs = np.arange(wp - wp / 2, 3 * wp, 0.01)
        time = np.arange(0, 1000, 0.1)

        spectrum = ModifiedPiersonMoskowitz(freqs, False)
        timeseries = spectrum.realization(time, hs=hs, tp=tp)

        hs_estimated = 4 * np.sqrt(np.average(timeseries**2))

        assert np.isclose(hs, hs_estimated, rtol=1e-1)

    def test_realization_jonswap(self):
        # Test that estimated hs of time realization is close to that of the specified hs
        hs = 3
        tp = 12
        wp = 2 * np.pi / tp
        gamma = 3.3

        freqs = np.arange(wp - wp / 2, 3 * wp, 0.01)
        time = np.arange(0, 1000, 0.1)

        spectrum = JONSWAP(freqs)
        timeseries = spectrum.realization(time, hs=hs, tp=tp, gamma=gamma)

        hs_estimated = 4 * np.sqrt(np.average(timeseries**2))

        assert np.isclose(hs, hs_estimated, rtol=1e-1)

    def test_equal_spectra(self):
        hs = 3
        tp = 12
        wp = 2 * np.pi / tp
        gamma = 1

        freqs = np.arange(wp - wp / 2, 3 * wp, 0.01)

        jonswap_spectrum = JONSWAP(freqs)
        pm_spectrum = ModifiedPiersonMoskowitz(freqs)

        assert np.all(
            np.isclose(jonswap_spectrum(hs, tp, gamma), pm_spectrum(hs, tp), rtol=1e-1)
        )

    def test_area_spreading_function(self):
        # Verify that the integral of the spreading function equals to 1
        N = 100
        theta_min, theta_max = -np.pi, np.pi
        theta = np.linspace(theta_min, theta_max, N)

        spreading_spectra = MultiDirectional(theta)
        angles, spreading = spreading_spectra()
        try:
            integrate_func = np.trapezoid  # New function in latest NumPy versions
        except AttributeError:
            integrate_func = np.trapz  # Fallback for older versions
        integral = integrate_func(spreading, angles)
        print(integral)
        assert np.isclose(integral, 1, rtol=1e-2)

    def test_area_spreading_function_alt(self):
        N = 100
        theta_min, theta_max = -np.pi, np.pi
        theta = np.linspace(theta_min, theta_max, N)
        theta0 = np.pi / 4
        s = 2

        spreading_spectra = MultiDirectionalAlt(theta)
        angles, spreading = spreading_spectra(theta0, s)

        integral = np.trapz(spreading, angles)
        print(integral)
        assert np.isclose(integral, 1, rtol=1e-2)
