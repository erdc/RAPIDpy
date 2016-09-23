Postprocessing
==============

Merge RAPID Output
------------------

You can use this to combine consecutive RAPID output files into one
file.

WARNING: This code replaces the first file with the combined output and
deletes the second file. BACK UP YOUR FILES!!!!

Example:

.. code:: python

    import datetime
    from RAPIDpy.postprocess import ConvertRAPIDOutputToCF
    file1 = ""
    file2 = ""
    cv = ConvertRAPIDOutputToCF(rapid_output_file=[file1, file2],
                                start_datetime=datetime.datetime(2005,1,1),
                                time_step=[3*3600, 3*3600],
                                qinit_file="", #optional
                                comid_lat_lon_z_file="", #optional
                                rapid_connect_file="", #optional
                                project_name="NLDAS(VIC)-RAPID historical flows by US Army ERDC",
                                output_id_dim_name='COMID',
                                output_flow_var_name='Qout',
                                print_debug=False)
    cv.convert()

Generate qinit from past qout
-----------------------------

RAPIDpy also creates a qinit file from a RAPID qout file. This example
shows how.

.. code:: python

    from RAPIDpy import RAPID
    rapid_manager = RAPID(Qout_file='/output_mississippi-nfie/Qout_k2v1_2005to2009.nc', 
                          rapid_connect_file='/input_mississippi_nfie/rapid_connect_ECMWF.csv'
                         )
    rapid_manager.generate_qinit_from_past_qout(qinit_file='/input_mississippi_nfie/Qinit_2008_flood.csv',
                            time_index=10162) #time_index is optional, if not included it will be last time step

Generate seasonal qinit from past qout
--------------------------------------

RAPIDpy also creates a seasonal qinit file from a RAPID qout file. This
requires a simulation Qout file with a longer time period of record and
to be CF compliant. It takes the average of the current date +- 3 days
and goes back as far as possible. This example shows how to use it.

.. code:: python

    from RAPIDpy.rapid import RAPID
    rapid_manager = RAPID(Qout_file='/output_mississippi-nfie/Qout_2000to2015.nc', 
                          rapid_connect_file='/input_mississippi_nfie/rapid_connect_ECMWF.csv'
                          )

    rapid_manager.generate_seasonal_intitialization(qinit_file='/input_mississippi_nfie/Qinit_seasonal_avg_jan_1.csv',
                                                    #datetime_start_initialization, #this is optional, default is utcnow
                                                    )
                                                    

Goodness of Fit
---------------

To check how well your simulation performed versus observations, this
function can help you.

.. code:: python

    from RAPIDpy.postprocess import find_goodness_of_fit_csv, find_goodness_of_fit

    #if you have observations and model results next to each other in columns, this works for you
    find_goodness_of_fit_csv('/united_kingdom-thames/flows_kingston_gage_noah.csv')

    #if you have gage files in the format RAPID needs for calibration, this will generate a file for you with results
    reach_id_file = os.path.join(INPUT_DATA_PATH, 'obs_reach_id.csv') 
    observed_file = os.path.join(INPUT_DATA_PATH, 'obs_flow.csv') 

    #using CF-compliant file
    cf_input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_CF.nc')
    cf_out_analysis_file = os.path.join(OUTPUT_DATA_PATH, 'cf_goodness_of_fit_results-daily.csv') 
    find_goodness_of_fit(cf_input_qout_file, reach_id_file, observed_file,
                         cf_out_analysis_file, daily=True)

    #using original RAPID file
    original_input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_original.nc')
    original_out_analysis_file = os.path.join(OUTPUT_DATA_PATH, 'original_goodness_of_fit_results-daily.csv') 
    find_goodness_of_fit(rapid_qout_file=original_input_qout_file, 
                         reach_id_file=reach_id_file, 
                         observed_file=observed_file,
                         out_analysis_file=original_out_analysis_file, 
                         steps_per_group=8) #for raw rapid output (8 is produces daily flows for 3-hr timesteps)
