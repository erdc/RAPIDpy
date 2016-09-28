RAPID to GSSHA
==============

It is possible to use RAPID streamflow as a boundary condition to the 
Gridded Surface Subsurface Hydrologic Analysis (GSSHA) model.

What is GSSHA?
--------------

GSSHA is a physically-based, distributed hydrologic model. GSSHA is developed 
and maintained by Coastal and Hydraulics Laboratory (CHL) which is
a member of the Engineer Research & Development Center of the United
States Army Corps of Engineers (USACE).

.. note::
	
	For more information about GSSHA please visit the the gsshawiki_ .

.. _gsshawiki: http://www.gsshawiki.com/Main_Page

Tutorial
--------

Step 1: Generate XYS File
~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: RAPIDpy.dataset.RAPIDDataset.write_flows_to_gssha_time_series_xys(path_to_output_file,series_name,series_id,reach_index=None,reach_id=None,date_search_start=None,date_search_end=None,daily=False,mode="mean")
    :noindex:

Step 2: Determine Link & Node Connection to RAPID river ID
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step 2.1: Add XYS File in WMS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
See Aquaveo's WMS Tutorial: http://wmstutorials-10.1.aquaveo.com/55%20Gssha-Applications-OverlandBoundaryConditions.pdf

Step 2.2: Look at Generated IHG File to find Link & Node
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Step 3: Generate other IHG Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: RAPIDpy.dataset.RAPIDDataset.write_flows_to_gssha_time_series_ihg(path_to_output_file,point_list,date_search_start=None,date_search_end=None,daily=False,mode="mean")
    :noindex:

