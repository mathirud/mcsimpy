.. mclsimpy documentation master file, created by
   sphinx-quickstart on Fri Jan 20 09:45:08 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to mclsimpy's documentation!
=========================================

|License: GPL v3| |Python: v3.8| |Python: v3.9| |Python: v3.10| |Python: v3.11| |Python: v3.12| |Python: v3.13|

.. |License: GPL v3| image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0

.. |Python: v3.8| image:: https://shields.io/badge/Python-v3.8-green.svg
   :target: https://www.python.org/downloads/release/python-380/

.. |Python: v3.9| image:: https://shields.io/badge/Python-v3.9-green.svg
   :target: https://www.python.org/downloads/release/python-390

.. |Python: v3.10| image:: https://shields.io/badge/Python-v3.10-green.svg
   :target: https://www.python.org/downloads/release/python-3100

.. |Python: v3.11| image:: https://shields.io/badge/Python-v3.11-green.svg
   :target: https://www.python.org/downloads/release/python-3110

.. |Python: v3.12| image:: https://shields.io/badge/Python-v3.12-green.svg
   :target: https://www.python.org/downloads/release/python-3120

.. |Python: v3.13| image:: https://shields.io/badge/Python-v3.13-green.svg
   :target: https://www.python.org/downloads/release/python-3130


The ``mclsimpy`` package is developed as a part of the master thesis of Jan-Erik Hygen, Marie Kongshaug, and Harald Mo, at the
`Department of Marine Technology (NTNU, Trondheim) <https://www.ntnu.edu/imt>`_. The package is used for simulation of DP stationkeeping models and
manuevering models subject to first- and second order wave loads.

The package have some of the same features as the `MSS Toolbox <https://se.mathworks.com/matlabcentral/fileexchange/86393-marine-systems-simulator-mss?requestedDomain=>`_ in Matlab/Simulink.

.. toctree::
   :maxdepth: 1
   :caption: Contents

   getting_started
   modules/simulator
   modules/control
   modules/guidance
   modules/thrust_allocation
   modules/observer
   modules/waves
   modules/utils

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/simulation
   tutorials/waveload

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
