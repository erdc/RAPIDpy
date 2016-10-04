Postprocessing
==============

Merge RAPID Output
------------------

.. autoclass:: RAPIDpy.postprocess.ConvertRAPIDOutputToCF

Generate qinit from past qout
-----------------------------

RAPIDpy also creates a qinit file from a RAPID qout file. This example
shows how.

.. automethod:: RAPIDpy.rapid.RAPID.generate_qinit_from_past_qout(qinit_file, time_index=-1)
    :noindex:                                                

Generate seasonal qinit from past qout
--------------------------------------

.. automethod:: RAPIDpy.rapid.RAPID.generate_seasonal_intitialization(qinit_file,datetime_start_initialization=datetime.datetime.utcnow())
    :noindex:                                                

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
