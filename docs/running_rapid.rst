Running RAPID
=============

Step 1: Initialize the RAPID manager class.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  First, add the path to you rapid executable location.
-  Next, you need to either tell it how many processors to use using the
   *num\_processors* input variable or to use all available processors
   set *use\_all\_processors* to true.
-  After that, add any other parameters you would like to use that would
   normally be in the rapid namelist file (this is case sensitive).

Example:

.. code:: python

    from RAPIDpy import RAPID
    rapid_manager = RAPID(rapid_executable_location='~/work/rapid/run/rapid'
                          use_all_processors=True, #optional, default is False
                          #num_processors=1, #optional, default is 1, overridden if use_all_processors is True                         
                          #ksp_type="richardson", #optional, default is richardson
                          ZS_TauR=24*3600, #duration of routing procedure (time step of runoff data)
                          ZS_dtR=15*60, #internal routing time step
                          ZS_TauM=365*24*3600, #total simulation time 
                          ZS_dtM=24*3600 #input time step 
                         )

If you are using Cygwin on Windows:

.. code:: python

    from RAPIDpy import RAPID
    rapid_manager = RAPID(rapid_executable_location='C:\\cygwin64\\home\\username\\work\\rapid\\run\\rapid',
                          cygwin_bin_location='C:\\cygwin64\\bin',
                          use_all_processors=True, #optional, default is False
                          #num_processors=1, #optional, default is 1, overridden if use_all_processors is True                         
                          #ksp_type="richardson", #optional, default is richardson
                          ZS_TauR=24*3600, #duration of routing procedure (time step of runoff data)
                          ZS_dtR=15*60, #internal routing time step
                          ZS_TauM=365*24*3600, #total simulation time 
                          ZS_dtM=24*3600 #input time step 
                         )

Step 2 (optional): Add/update additional namelist parameters later
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can add or update rapid namelist parameters using the
*update\_parameters* function by using the name of the variable in the
rapid namelist file (this is case sensitive).

Example:

.. code:: python

    rapid_manager.update_parameters(rapid_connect_file='../rapid_input_directory/rapid_connect.csv',
                                    Vlat_file='../rapid_input_directory/m3_riv.nc',
                                    riv_bas_id_file='../rapid_input_directory/riv_bas_id.csv,
                                    k_file='../rapid_input_directory/k.csv',
                                    x_file='../rapid_input_directory/x.csv',
                                    Qout_file='../OUTPUT/Qout.nc',
                                    )

Step 3 (optional): Update reach number data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you don't want to manually count the numbers for the rapid\_connect
or riv\_bas\_id files, use the *update\_reach\_number\_data* function.

Example:

.. code:: python

    rapid_manager.update_reach_number_data()

Step 4 (optional): Update simulation runtime data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you don't want to manually determine the total simulation duration,
use the *update\_simulation\_runtime* function.

Example:

.. code:: python

    rapid_manager.update_simulation_runtime()

Step 5: Run RAPID
~~~~~~~~~~~~~~~~~

This will generate your rapid\_namelist file and run RAPID from wherever
you call this script (your working directory).

Example:

.. code:: python

    rapid_manager.run()

Step 6 (optional): Convert RAPID output to be CF Compliant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This will convert the RAPID output to be CF compliant. This will require
a comid\_lat\_lon\_z file. Additionally, it prepends time zero to you
simulation. If no qinit file is given, a value of zero is added.

Example:

.. code:: python

    rapid_manager.make_output_CF_compliant(simulation_start_datetime=datetime.datetime(1980, 1, 1),
                                           comid_lat_lon_z_file='/rapid_input_directory/comid_lat_lon_z.csv', #optional
                                           project_name="ERA Interim Historical flows by US Army ERDC") 
