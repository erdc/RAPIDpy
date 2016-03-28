# -*- coding: utf-8 -*-
##
##  generate_return_periods.py
##  spt_lsm_autorapid_process
##
##  Created by Alan D. Snow and Scott D. Christensen.
##  Copyright © 2015-2016 Alan D Snow and Scott D. Christensen. All rights reserved.
##  License: BSD-3 Clause

from datetime import datetime
import netCDF4 as nc
import numpy as np
from RAPIDpy import RAPIDDataset

def generate_return_periods(qout_file, return_period_file, storm_duration_days=7):
    """
    Generate return period from RAPID Qout file
    """

    #get ERA Interim Data Analyzed
    with RAPIDDataset(qout_file) as qout_nc_file:
        print "Setting up Return Periods File ..."
        return_period_nc = nc.Dataset(return_period_file, 'w')
        
        return_period_nc.createDimension('rivid', qout_nc_file.size_river_id)

        timeSeries_var = return_period_nc.createVariable('rivid', 'i4', ('rivid',))
        timeSeries_var.long_name = (
            'Unique NHDPlus COMID identifier for each river reach feature')

        max_flow_var = return_period_nc.createVariable('max_flow', 'f8', ('rivid',))
        return_period_20_var = return_period_nc.createVariable('return_period_20', 'f8', ('rivid',))
        return_period_10_var = return_period_nc.createVariable('return_period_10', 'f8', ('rivid',))
        return_period_2_var = return_period_nc.createVariable('return_period_2', 'f8', ('rivid',))

        lat_var = return_period_nc.createVariable('lat', 'f8', ('rivid',),
                                                  fill_value=-9999.0)
        lat_var.long_name = 'latitude'
        lat_var.standard_name = 'latitude'
        lat_var.units = 'degrees_north'
        lat_var.axis = 'Y'

        lon_var = return_period_nc.createVariable('lon', 'f8', ('rivid',),
                                                  fill_value=-9999.0)
        lon_var.long_name = 'longitude'
        lon_var.standard_name = 'longitude'
        lon_var.units = 'degrees_east'
        lon_var.axis = 'X'

        return_period_nc.variables['lat'][:] = qout_nc_file.qout_nc.variables['lat'][:]
        return_period_nc.variables['lon'][:] = qout_nc_file.qout_nc.variables['lon'][:]

        river_id_list = qout_nc_file.get_river_id_array()
        return_period_nc.variables['rivid'][:] = river_id_list

        print "Extracting Data and Generating Return Periods ..."
        time_array = qout_nc_file.get_time_array()
        num_years = int((datetime.utcfromtimestamp(time_array[-1])-datetime.utcfromtimestamp(time_array[0])).days/365.2425)
        time_steps_per_day = (24*3600)/float((datetime.utcfromtimestamp(time_array[1])-datetime.utcfromtimestamp(time_array[0])).seconds)
        step = max(1,int(time_steps_per_day * storm_duration_days))

        for comid_index, comid in enumerate(river_id_list):

            filtered_flow_data = qout_nc_file.get_daily_qout_index(comid_index, 
                                                                   steps_per_group=step,
				                                   mode="max")
            sorted_flow_data = np.sort(filtered_flow_data)[:num_years:-1]

            rp_index_20 = int((num_years + 1)/20.0)
            rp_index_10 = int((num_years + 1)/10.0)
            rp_index_2 = int((num_years + 1)/2.0)
            
            max_flow_var[comid_index] = sorted_flow_data[0]
            return_period_20_var[comid_index] = sorted_flow_data[rp_index_20]
            return_period_10_var[comid_index] = sorted_flow_data[rp_index_10]
            return_period_2_var[comid_index] = sorted_flow_data[rp_index_2]

        return_period_nc.close()
