Inflow from Land Surface Models
===============================

Code to use to prepare input data for RAPID from Land Surface Models
(LSM) such as:

- ECMWF's ERA Interim Data
- NASA's GLDAS/NLDAS/LIS Data

Step 1: Retrieve Land Surface Model Runoff Output
-------------------------------------------------

Download the data into a local directory.

- http://apps.ecmwf.int/datasets
- http://ldas.gsfc.nasa.gov/index.php


Step 2: Create folders for RAPID input and output
-------------------------------------------------

In this instance:

::

    $ cd $HOME
    $ mkdir -p rapid-io/input rapid-io/output

Step 3: Create script using LSM process
---------------------------------------

Here is the API for the function to use. Follow the example to create a script to use the code such as in ~/run\_lsm.py.

.. autofunction:: RAPIDpy.inflow.lsm_rapid_process.run_lsm_rapid_process

Step 4: Add RAPID files to the rapid-io/input directory
---------------------------------------------------------

Make sure the directory is in the format [watershed name]-[subbasin
name] with lowercase letters, numbers, and underscores only. No spaces!

Example:

::

    $ ls /rapid-io/input
    nfie_texas_gulf_region-huc_2_12

    $ ls /rapid-io/input/nfie_texas_gulf_region-huc_2_12
    comid_lat_lon_z.csv
    k.csv
    rapid_connect.csv
    riv_bas_id.csv
    weight_era_t511.csv
    weight_nldas.csv
    weight_gldas.csv
    weight_lis.csv
    weight_wrf.csv
    x.csv

If you have not generated these files yet, see :doc:`gis_tools` 

Step 5: Run the code
--------------------

::

    $ python ~/run_lsm.py
