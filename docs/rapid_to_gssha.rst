RAPID to GSSHA
==============

It is possible to use RAPID streamflow as an overland flow boundary condition  
to the Gridded Surface Subsurface Hydrologic Analysis (GSSHA) model.

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
There are two ways to input RAPID as a boundary condition for GSSHA. 
One is to connect the GSSHA stream network link and node to the RAPID
river ID and generate the IHG file. The other is to generate an XYS 
timeseries file and add it to the netork using WMS.  


Method 1: Generate IHG File
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step 1.1: Look at Stream Network in WMS to find Link & Node
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1. Open GSSHA project in WMS

    - Switch to 2-D Grid Module
    - In the top menu: *GSSHA -> Open Project File*

2. Turn on *Stream Link Numbers* Display

    - In the top menu: *Display -> Display Options*
    - Select *Map Data* in top left box
    - In the center box, make sure *Stream Link Numbers* is checked under the Arcs subsection.

3. Determin the Link ID by looking on model.

Step 1.2: Connect RAPID river ID to GSSHA Link & Node and Generate IHG
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: RAPIDpy.dataset.RAPIDDataset.write_flows_to_gssha_time_series_ihg
    :noindex:


Method 2: Generate XYS File
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step 2.1: Generate XYS File
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: RAPIDpy.dataset.RAPIDDataset.write_flows_to_gssha_time_series_xys
    :noindex:

Step 2.2: Add XYS File in WMS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the Distributed Hydrology section, go to *Overland Flow Boundary Conditions in GSSHA (Tutorial 55)* at http://www.aquaveo.com/software/wms-learning-tutorials

Here is a direct link to the document: http://wmstutorials-10.1.aquaveo.com/55%20Gssha-Applications-OverlandBoundaryConditions.pdf