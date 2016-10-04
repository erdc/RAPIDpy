Generating Stream Network
=========================

Using ArcHydro to Generate Stream Network
-----------------------------------------

See: 

- https://github.com/Esri/python-toolbox-for-rapid
- https://github.com/erdc-cm/python-toolbox-for-rapid

Using TauDEM to Generate Stream Network
---------------------------------------

For more information about taudem, see:
http://hydrology.usu.edu/taudem/taudem5/index.html 

Installation
------------

Step 1: Install TauDEM
~~~~~~~~~~~~~~~~~~~~~~

::

    $ cd /path/to/scripts
    $ git clone https://github.com/dtarb/TauDEM.git
    $ cd TauDEM/src
    $ make

Step 2: Install RAPIDpy with GIS Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See: :doc:`installation`

How To Use
----------

Initialize TauDEM Manager
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: RAPIDpy.gis.taudem.TauDEM

Generate network from DEM
~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: RAPIDpy.gis.taudem.TauDEM.demToStreamNetwork(output_directory,raw_elevation_dem="",pit_filled_elevation_grid="",flow_dir_grid_d8="",contributing_area_grid_d8="",flow_dir_grid_dinf="",contributing_area_grid_dinf="",use_dinf=False,threshold=1000,delineate=False)

Add Length in meters attribute
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: RAPIDpy.gis.taudem.TauDEM.addLengthMeters(stream_network)

Extract Sub Network
~~~~~~~~~~~~~~~~~~~

STEP 1: Extract sub network from stream network
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
There are two options to do this.

1. Choose your own outlet point: :func:`~RAPIDpy.gis.taudem.TauDEM.extractSubNetwork()`
2. Or let the code find the larges network: :func:`~RAPIDpy.gis.taudem.TauDEM.extractLargestSubNetwork()`.

.. automethod:: RAPIDpy.gis.taudem.TauDEM.extractSubNetwork(network_file,out_subset_network_file,outlet_ids,river_id_field,next_down_id_field,river_magnitude_field,safe_mode=True)

.. automethod:: RAPIDpy.gis.taudem.TauDEM.extractLargestSubNetwork(network_file,out_subset_network_file,river_id_field,next_down_id_field,river_magnitude_field,safe_mode=True)

STEP 2: Extract sub network from catchments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: RAPIDpy.gis.taudem.TauDEM.extractSubsetFromWatershed(subset_network_file,subset_network_river_id_field,watershed_file,watershed_network_river_id_field,out_watershed_subset_file)
