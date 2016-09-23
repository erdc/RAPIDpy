Inflow from Land Surface Models
===============================

Code to use to prepare input data for RAPID from Land Surface Models
(LSM) such as ECMWF ERA Interim Data or NASA GLDAS/NLDAS/LIS data.

Step 1: Create folders for RAPID input and output
-------------------------------------------------

In this instance:

::

    $ cd $HOME
    $ mkdir -p rapid-io/input rapid-io/output

Step 2: Create script using LSM process
---------------------------------------

Here is an example script to use the code (ex. ~/run\_lsm.py).

.. code:: python

    from datetime import datetime
    from RAPIDpy.inflow import run_lsm_rapid_process
    #------------------------------------------------------------------------------
    #main process
    #------------------------------------------------------------------------------
    if __name__ == "__main__":
        run_lsm_rapid_process(
            rapid_executable_location='/home/alan/rapid/src/rapid',
            rapid_io_files_location='/home/alan/rapid-io',
            lsm_data_location='/home/alan/era_data', #path to folder with LSM data
            simulation_start_datetime=datetime(1980, 1, 1),
            simulation_end_datetime=datetime(2014, 12, 31),
            generate_rapid_namelist_file=True, #if you want to run RAPID manually later
            run_rapid_simulation=True, #if you want to run RAPID after generating inflow file
            generate_return_periods_file=False, #if you want to get return period file from RAPID simulation
            generate_seasonal_initialization_file=False, #if you want to get seasonal init file from RAPID simulation
            generate_initialization_file=False, #if you want to generate qinit file from end of RAPID simulation
            use_all_processors=True, #defaults to use all processors available
            num_processors=1, #you can change this number if use_all_processors=False
            cygwin_bin_location="" #if you are using Windows with Cygwin
        )

Step 3: Add RAPID files to the work/rapid/input directory
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
    x.csv

Step 4: Run the code
--------------------

::

    $ python ~/run_lsm.py
