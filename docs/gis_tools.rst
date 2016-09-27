RAPID GIS Tools
===============

These tools generate the RAPID input files and weight table files from the GIS stream networks.


.. note:: To generate your own network from a DEM see :doc:`gis_stream_network`

.. note:: For these tools to work, you need GIS dependencies installed (See :doc:`installation`).

There are also tools by Esri for ArcMap located here:

- https://github.com/Esri/python-toolbox-for-rapid
- https://github.com/erdc-cm/python-toolbox-for-rapid

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


How it works:
-------------

Snow, Alan D., Scott D. Christensen, Nathan R. Swain, E. James Nelson,
Daniel P. Ames, Norman L. Jones, Deng Ding, Nawajish S. Noman, Cedric H.
David, Florian Pappenberger, and Ervin Zsoter, 2016. A High-Resolution
National-Scale Hydrologic Forecast System from a Global Ensemble Land
Surface Model. *Journal of the American Water Resources Association
(JAWRA)* 1-15, DOI: 10.1111/1752-1688.12434

Snow, Alan Dee, "A New Global Forecasting Model to Produce
High-Resolution Stream Forecasts" (2015). All Theses and Dissertations.
Paper 5272. http://scholarsarchive.byu.edu/etd/5272
