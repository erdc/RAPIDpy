# -*- coding: utf-8 -*-
##
##  generate_return_periods.py
##  RAPIDpy
##
##  Created by Alan D. Snow and Scott D. Christensen.
##  Copyright © 2015-2016 Alan D Snow and Scott D. Christensen. All rights reserved.
##  License: BSD-3 Clause

from datetime import datetime
import multiprocessing
import netCDF4 as nc
import numpy as np

#local
from ..dataset import RAPIDDataset
from ..utilities import partition

def generate_single_return_period(args):
    """
    This function calculates a single return period for a single reach
    """
    qout_file=args[0]
    return_period_file=args[1]
    rivid_index_list=args[2]
    step=args[3]
    num_years=args[4]
    mp_lock=args[5]
    
    with RAPIDDataset(qout_file) as qout_nc_file: 
        #get index of return period data
        rp_index_20 = int((num_years + 1)/20.0)
        rp_index_10 = int((num_years + 1)/10.0)
        rp_index_2 = int((num_years + 1)/2.0)
        
        #iterate through rivids to generate return periods
        max_flow_array = np.zeros(len(rivid_index_list))
        return_20_array = np.zeros(len(rivid_index_list))
        return_10_array = np.zeros(len(rivid_index_list))
        return_2_array = np.zeros(len(rivid_index_list))
        
        for iter_idx, rivid_index in enumerate(rivid_index_list):
            filtered_flow_data = qout_nc_file.get_qout_index(rivid_index,
                                                             pd_filter="{0}D".format(step),
                                                             filter_mode="max")
                                                             
            sorted_flow_data = np.sort(filtered_flow_data)[:num_years:-1]
            max_flow = sorted_flow_data[0]
            if max_flow < 0.01:
                print("WARNING: Return period data < 0.01 generated for rivid {0}" \
                      .format(qout_nc_file.qout_nc.variables[qout_nc_file.river_id_dimension][rivid_index]))
            max_flow_array[iter_idx] = max_flow
            return_20_array[iter_idx] = sorted_flow_data[rp_index_20]
            return_10_array[iter_idx] = sorted_flow_data[rp_index_10]
            return_2_array[iter_idx] = sorted_flow_data[rp_index_2]

        mp_lock.acquire()
        return_period_nc = nc.Dataset(return_period_file, 'a')
        return_period_nc.variables['max_flow'][rivid_index_list] = max_flow_array
        return_period_nc.variables['return_period_20'][rivid_index_list] = return_20_array
        return_period_nc.variables['return_period_10'][rivid_index_list] = return_10_array
        return_period_nc.variables['return_period_2'][rivid_index_list] = return_2_array
        return_period_nc.close()
        mp_lock.release()

def generate_return_periods(qout_file, return_period_file, num_cpus=multiprocessing.cpu_count(), storm_duration_days=7):
    """
    Generate return period from RAPID Qout file
    """

    #get ERA Interim Data Analyzed
    with RAPIDDataset(qout_file) as qout_nc_file:
        print("Setting up Return Periods File ...")
        return_period_nc = nc.Dataset(return_period_file, 'w')
        
        return_period_nc.createDimension('rivid', qout_nc_file.size_river_id)

        timeSeries_var = return_period_nc.createVariable('rivid', 'i4', ('rivid',))
        timeSeries_var.long_name = (
            'unique identifier for each river reach')

        max_flow_var = return_period_nc.createVariable('max_flow', 'f8', ('rivid',))
        max_flow_var.long_name = 'maxumum streamflow'
        max_flow_var.units = 'm3/s'
        
        return_period_20_var = return_period_nc.createVariable('return_period_20', 'f8', ('rivid',))
        return_period_20_var.long_name = '20 year return period flow'
        return_period_20_var.units = 'm3/s'
        
        return_period_10_var = return_period_nc.createVariable('return_period_10', 'f8', ('rivid',))
        return_period_10_var.long_name = '10 year return period flow'
        return_period_10_var.units = 'm3/s'
        
        return_period_2_var = return_period_nc.createVariable('return_period_2', 'f8', ('rivid',))
        return_period_2_var.long_name = '2 year return period flow'
        return_period_2_var.units = 'm3/s'

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
        return_period_nc.close()

        time_array = qout_nc_file.get_time_array()
        
    print("Extracting Data and Generating Return Periods ...")
    num_years = int((datetime.utcfromtimestamp(time_array[-1])-datetime.utcfromtimestamp(time_array[0])).days/365.2425)
    time_steps_per_day = (24*3600)/float((datetime.utcfromtimestamp(time_array[1])-datetime.utcfromtimestamp(time_array[0])).total_seconds())
    step = max(1,int(time_steps_per_day * storm_duration_days))
    
    #generate multiprocessing jobs
    mp_lock = multiprocessing.Manager().Lock()
    job_combinations = []
    partition_list, partition_index_list = partition(river_id_list, num_cpus*2)
    for sub_partition_index_list in partition_index_list:
        if len(sub_partition_index_list) > 0:
            job_combinations.append((qout_file,
                                     return_period_file,
                                     sub_partition_index_list, 
                                     step,
                                     num_years,
                                     mp_lock
                                     ))

    pool = multiprocessing.Pool(num_cpus)
    pool.map(generate_single_return_period,
             job_combinations)
    pool.close()
    pool.join()