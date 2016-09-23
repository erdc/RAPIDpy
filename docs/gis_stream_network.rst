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

.. code:: python

    from RAPIDpy.gis.taudem import TauDEM
    td = TauDEM("/path/to/scripts/TauDEM",
                num_processors=1, #default is one, can change this for processing
                use_all_processors=False, #default is False, if True it will use all cores
                )

Generate network from DEM

.. code:: python

    elevation_dem = '/path/to/dem.tiff'
    output_directory = '/path/to/output/files'
    td.demToStreamNetwork(elevation_dem, 
                          output_directory, 
                          threshold=1000, #See: http://hydrology.usu.edu/taudem/taudem5/help53/StreamDefinitionByThreshold.html
                          )

Add Length in meters attribute

.. code:: python

    td.addLengthMeters(os.path.join(output_directory,"stream_reach_file.shp"))

Extract Sub Network

.. code:: python

    #to let the code find the largest network
    td.extractLargestSubNetwork(network_file=os.path.join(output_directory,"stream_reach_file.shp"),                                         
                                out_subset_network_file=os.path.join(output_directory,"stream_reach_file_subset.shp"),
                                river_id_field="LINKNO",
                                next_down_id_field="DSLINKNO",
                                river_magnitude_field="Magnitude",
                                safe_mode=True, #this will prevent your machine from crashing. If you are sure it will work, set to False.
                                )
    #to extract a specific network
    td.extractLargestSubNetwork(network_file=os.path.join(output_directory,"stream_reach_file.shp"),                                         
                                out_subset_network_file=os.path.join(output_directory,"stream_reach_file_subset.shp"),
                                outlet_ids=[60830], #list of outlet ids
                                river_id_field="LINKNO",
                                next_down_id_field="DSLINKNO",
                                river_magnitude_field="Magnitude",
                                safe_mode=True, #this will prevent your machine from crashing. If you are sure it will work, set to False.
                                )
    #to extract the subset watersheds using subset river network
    td.extractSubsetFromWatershed(subset_network_file=os.path.join(output_directory,"stream_reach_file_subset.shp"),
                                  subset_network_river_id_field="LINKNO",
                                  watershed_file=os.path.join(output_directory,"watershed_shapefile.shp"),
                                  watershed_network_river_id_field="LINKNO",
                                  out_watershed_subset_file=os.path.join(output_directory,"watershed_shapefile_subset.shp"))
