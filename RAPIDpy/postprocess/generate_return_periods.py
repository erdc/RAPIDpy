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
    method=args[5]
    mp_lock=args[6]
    
    with RAPIDDataset(qout_file) as qout_nc_file: 
        #get index of return period data
        rp_index_20 = int((num_years + 1)/20.0)
        rp_index_10 = int((num_years + 1)/10.0)
        rp_index_2 = int((num_years + 1)/2.0)
        
        #iterate through rivids to generate return periods
        max_flow_array = np.zeros(len(rivid_index_list))
        if method in ('gumble', 'log_pearson', 'gev'):
            return_100_array = np.zeros(len(rivid_index_list))
            return_50_array = np.zeros(len(rivid_index_list))
        return_20_array = np.zeros(len(rivid_index_list))
        return_10_array = np.zeros(len(rivid_index_list))
        return_2_array = np.zeros(len(rivid_index_list))
        
        for iter_idx, rivid_index in enumerate(rivid_index_list):
            filtered_flow_data = qout_nc_file.get_qout_index(rivid_index,
                                                             pd_filter="{0}D".format(step),
                                                             filter_mode="max")

            if method == 'weibull':

                sorted_flow_data = np.sort(filtered_flow_data)[:num_years:-1]
                max_flow = sorted_flow_data[0]
                if max_flow < 0.01:
                    print("WARNING: Return period data < 0.01 generated for rivid {0}" \
                          .format(qout_nc_file.qout_nc.variables[qout_nc_file.river_id_dimension][rivid_index]))
                max_flow_array[iter_idx] = max_flow
                return_20_array[iter_idx] = sorted_flow_data[rp_index_20]
                return_10_array[iter_idx] = sorted_flow_data[rp_index_10]
                return_2_array[iter_idx] = sorted_flow_data[rp_index_2]

            elif method == 'gumble'

                mean_flow = np.mean(filtered_flow_data)
                stddev = np.std(filtered_flow_data)
                if mean_flow < 0.01:
                    print("WARNING: Return period data < 0.01 generated for rivid {0}" \
                          .format(qout_nc_file.qout_nc.variables[qout_nc_file.river_id_dimension][rivid_index]))
                return_100_array[iter_idx] = mean_flow + 3.14*stddev
                return_50_array[iter_idx] = mean_flow + 2.59*stddev
                return_20_array[iter_idx] = mean_flow + 1.87*stddev
                return_10_array[iter_idx] = mean_flow + 1.3*stddev
                return_2_array[iter_idx] = mean_flow - .164*stddev



            elif method == 'log_pearson'

                log_flow = numpy.log10(filtered_flow_data)
                mean_log_flow = np.mean(log_flow)
                std_log_flow = np.std(log_flow)
                log_flow_array = np.array(log_flow)
                skew = (num_years*(np.sum(np.power((log_flow_array - mean_log_flow),3))))/((num_years-1)*(num_years-2)*(std_log_flow)**3)
                K_table_2 = [-.396, -.384, -.368, -.351, -.330, -.307, -.282, -.254, -.225, -.195, -.164, -.132, -.099, -.066, -.033, 0, .033, .066, .099, .132, .164, .195, .225, .254, .282, .307, .330, .351, .368, .384, .396]
                K_table_10 = [1.18, 1.21, 1.238, 1.262, 1.284, 1.302, 1.318, 1.329, 1.337, 1.340, 1.340, 1.336, 1.328, 1.317, 1.301, 1.282, 1.258, 1.231, 1.2, 1.166, 1.128, 1.086, 1.041, .994, .945, .895, .844, .795, .747, .702, .660]
                K_table_50 = [3.152, 3.114, 3.071, 3.023, 2.97, 2.912, 2.848, 2.78, 2.706, 2.626, 2.542, 2.453, 2.359, 2.261, 2.159, 2.054, 1.945, 1.834, 1.72, 1.606, 1.492, 1.379, 1.270, 1.166, 1.069, .98, .9, .83, .768, .714, .666]
                K_table_100 = [4.051, 3.973, 3.889, 3.8, 3.705, 3.605, 3.499, 3.388, 3.271, 3.149, 3.022, 2.891, 2.755, 2.615, 2.472, 2.326, 2.178, 2.029, 1.88, 1.733, 1.588, 1.499, 1.318, 1.197, 1.087, .99, .905, .832, .769, .714, .667]
                return_100_array[iter_idx] =
                return_50_array[iter_idx] =
                return_20_array[iter_idx] =
                return_10_array[iter_idx] =
                return_2_array[iter_idx] =

        mp_lock.acquire()
        return_period_nc = nc.Dataset(return_period_file, 'a')
        return_period_nc.variables['max_flow'][rivid_index_list] = max_flow_array
        if method in ('gumble', 'log_pearson', 'gev'):
            return_period_nc.variables['return_period_100'][rivid_index_list] = return_100_array
            return_period_nc.variables['return_period_50'][rivid_index_list] = return_50_array
        return_period_nc.variables['return_period_20'][rivid_index_list] = return_20_array
        return_period_nc.variables['return_period_10'][rivid_index_list] = return_10_array
        return_period_nc.variables['return_period_2'][rivid_index_list] = return_2_array
        return_period_nc.close()
        mp_lock.release()

def generate_return_periods(qout_file, return_period_file, num_cpus=multiprocessing.cpu_count(), storm_duration_days=7, method='weibull'):
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

        if method in ('gumble', 'log_pearson', 'gev'):

            return_period_100_var = return_period_nc.createVariable('return_period_100', 'f8', ('rivid',))
            return_period_100_var.long_name = '100 year return period flow'
            return_period_100_var.units = 'm3/s'

            return_period_50_var = return_period_nc.createVariable('return_period_50', 'f8', ('rivid',))
            return_period_50_var.long_name = '50 year return period flow'
            return_period_50_var.units = 'm3/s'

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
                                     method,
                                     mp_lock
                                     ))

    pool = multiprocessing.Pool(num_cpus)
    pool.map(generate_single_return_period,
             job_combinations)
    pool.close()
    pool.join()