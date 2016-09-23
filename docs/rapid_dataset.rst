RAPIDDataset
============

This is a wrapper for the RAPID Qout netCDF file. Here are some basic
examples for useage.

.. code:: python

    from datetime import datetime
    from RAPIDpy.dataset import RAPIDDataset

    with RAPIDDataset('/path/to/Qout.nc') as qout_nc:

This example demonstrates how to retrieve or generate a time array to go
along with your RAPID streamflow series

.. code:: python

        #CF-Compliant
        time_array = qout_nc.get_time_array()
        #or, to get datetime array
        time_datetime = qout_nc.get_time_array(return_datetime=True)
        
        #Original Qout
        time_array = qout_nc.get_time_array(datetime_simulation_start=datetime(1980, 1, 1),
                                            simulation_time_step_seconds=3*3600)
        #or, to get datetime array
        time_datetime = qout_nc.get_time_array(datetime_simulation_start=datetime(1980, 1, 1),
                                               simulation_time_step_seconds=3*3600,
                                               return_datetime=True)

This example demonstrates how to get the river ids in the RAPID Qout
dataset.

.. code:: python

        river_ids = qout_nc.get_river_id_array()

This example demonstrates how to retrieve the streamflow associated with
the reach you are interested in

.. code:: python

        river_id = 500
        streamflow_array = qout_nc.get_qout(river_id)

This example demonstrates how to retrieve the streamflow within a date
range associated with the reach you are interested in

.. code:: python

        river_id = 500
        #CF-Compliant
        streamflow_array = qout_nc.get_qout(river_id,
                                            date_search_start=datetime(1985,1,1),
                                            date_search_end=datetime(1985,2,4))

        #Original RAPID Qout
        streamflow_array = qout_nc.get_qout(river_id,
                                            time_index_start=20,
                                            time_index_end=25)

This example demonstrates how to get daily streamflow averages as an
array

.. code:: python

        river_id = 500
        #CF-Compliant
        streamflow_array = qout_nc.get_daily_qout(river_id,
                                                  mode="mean", #(optional) default is mean, but max is another option
                                                  )
        #Original RAPID Qout
        streamflow_array = qout_nc.get_daily_qout(river_id,
                                                  mode="mean", #(optional) default is mean, but max is another option
                                                  steps_per_group=8, #average 8 timesteps together for 1 day
                                                  )
        #OR if you know the index
        river_index = qout_nc.get_river_index(river_id)
        #CF-Compliant
        streamflow_array = qout_nc.get_daily_qout_index(river_index)
        #Original RAPID Qout
        streamflow_array = qout_nc.get_daily_qout_index(river_index,
                                                        steps_per_group=8, #average 8 timesteps together for 1 day
                                                        )

Write Qout timeseries to csv file This function simplifies writing time
series for each stream reach to a csv file. NOTE: Need either
*reach\_id* or *reach\_index* parameter, but either can be used.

.. code:: python

        #for writing entire time series to file
        qout_nc.write_flows_to_csv('/timeseries/Qout_3624735.csv', 
                                   reach_id=3624735, #COMID or rivid
                                   ) 
        #if file is CF compliant, you can write out daily average
        qout_nc.write_flows_to_csv('/timeseries/Qout_daily.csv',
                                   reach_index=20, #index of COMID or rivid (example if you already know index instead if reach_id),
                                   daily=True, #if file is CF compliant, write out daily flows
                                   mode="mean" #(optional) default is "mean", however "max" is available
                                   )
        #if file is CF compliant, you can filter by date
        qout_nc.write_flows_to_csv('/timeseries/Qout_daily_date_filter.csv',
                                   reach_index=20, #index of COMID or rivid (example if you already know index instead if reach_id),
                                   daily=True, #if file is CF compliant, write out daily flows
                                   date_search_start=datetime(2002, 8, 31), #optional start date filter
                                   date_search_end=datetime(2002, 9, 15), #optional end date filter
                                   mode="max" #(optional) default is "mean", however "max" is available
                                   )
