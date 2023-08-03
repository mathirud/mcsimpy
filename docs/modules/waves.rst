=====
Waves
=====

Wave spectra and wave loads functionality can be found in the ``mcsimpy.waves`` subpackage.

Wave Spectra
------------

The ``mcsimpy.waves.wave_spectra`` module containts basic wave spectra functionality.

.. autoclass:: mcsimpy.waves.wave_spectra.BaseSpectrum
    :members:
    :special-members: __init__


Modified Pierson Moskowitz
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: mcsimpy.waves.wave_spectra.ModifiedPiersonMoskowitz
    :members:
    :special-members: __call__

JONSWAP
^^^^^^^
.. autoclass:: mcsimpy.waves.wave_spectra.JONSWAP
    :special-members: __call__



Wave Loads
----------

.. autoclass:: mcsimpy.waves.wave_loads.WaveLoad
    :members: first_order_loads, second_order_loads, QTF_METHODS
    :private-members: _set_force_raos, _full_qtf_6dof
    :special-members: __call__, __init__
