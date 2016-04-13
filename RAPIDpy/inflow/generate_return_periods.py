# -*- coding: utf-8 -*-
##
##  generate_return_periods.py
##  spt_lsm_autorapid_process
##
##  Created by Alan D. Snow and Scott D. Christensen.
##  Copyright © 2015-2016 Alan D Snow and Scott D. Christensen. All rights reserved.
##  License: BSD-3 Clause

import netCDF4 as nc
import numpy as np
from HydroStats.VDF import VDF

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
        return_period_100_var = return_period_nc.createVariable('return_period_100', 'f8', ('rivid',))
        return_period_50_var = return_period_nc.createVariable('return_period_50', 'f8', ('rivid',))
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
        time_array = qout_nc_file.get_time_array(return_datetime=True)
        
        for rivid_index, rivid in enumerate(river_id_list):
            streamflow = np.nan_to_num(qout_nc_file.get_qout_index(rivid_index))
            max_flow = np.amax(streamflow)
            if max_flow > 0.1:
                vdf_calc = VDF(time_array, 
                               streamflow,
                               np.array([storm_duration_days*24*60]),
                               ['kde', 'silverman'],
                               True
                               )
                               
##                vdf_calc = VDF(time_array, 
##                               streamflow,
##                               np.array([storm_duration_days*24*60]),
##                               ['gev'],
##                               True
##                               )
                max_flow_var[rivid_index] = np.amax(streamflow)
                return_period_100_var[rivid_index] = vdf_calc.calculate_return_period_curve(100)[0]
                return_period_50_var[rivid_index] = vdf_calc.calculate_return_period_curve(50)[0]
                return_period_20_var[rivid_index] = vdf_calc.calculate_return_period_curve(20)[0]
                return_period_10_var[rivid_index] = vdf_calc.calculate_return_period_curve(10)[0]
                return_period_2_var[rivid_index] = vdf_calc.calculate_return_period_curve(2)[0]
        else:
                return_period_100_var[rivid_index] = max_flow
                return_period_50_var[rivid_index] = max_flow
                return_period_20_var[rivid_index] = max_flow
                return_period_10_var[rivid_index] = max_flow
                return_period_2_var[rivid_index] = max_flow
            

        return_period_nc.close()
