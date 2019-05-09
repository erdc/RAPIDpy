RAPID GIS Tools
===============

These tools generate the RAPID input files and weight table files from the GIS stream networks.


.. note:: To generate your own network from a DEM see :doc:`gis_stream_network`

.. note:: For these tools to work, you need GIS dependencies installed (See :doc:`installation`).


Workflows
---------

Static RAPID Files
~~~~~~~~~~~~~~~~~~

.. autofunction:: RAPIDpy.gis.workflow.CreateAllStaticRAPIDFiles


Weight Table Files
~~~~~~~~~~~~~~~~~~

.. autofunction:: RAPIDpy.gis.workflow.CreateAllStaticECMWFFiles

Static RAPID Files and Weight Table Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: RAPIDpy.gis.workflow.CreateAllStaticECMWFRAPIDFiles


Individual Tools
----------------

Static RAPID Files
~~~~~~~~~~~~~~~~~~
.. autofunction:: RAPIDpy.gis.network.CreateNetworkConnectivity

.. autofunction:: RAPIDpy.gis.network.CreateNetworkConnectivityNHDPlus

.. autofunction:: RAPIDpy.gis.network.CreateSubsetFile

.. autofunction:: RAPIDpy.gis.muskingum.CreateMuskingumKfacFile

.. autofunction:: RAPIDpy.gis.muskingum.CreateMuskingumKFile

.. autofunction:: RAPIDpy.gis.muskingum.CreateMuskingumXFileFromDranageLine

.. autofunction:: RAPIDpy.gis.muskingum.CreateConstMuskingumXFile

Weight Tables
~~~~~~~~~~~~~
.. autofunction:: RAPIDpy.gis.weight.CreateWeightTableECMWF

.. autofunction:: RAPIDpy.gis.weight.CreateWeightTableLDAS


Utilities
---------

.. autofunction:: RAPIDpy.gis.centroid.FlowlineToPoint
