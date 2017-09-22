# -*- coding: utf-8 -*-
"""
    generate_return_periods.py
    RAPIDpy

    Created by: Alan D. Snow and Scott D. Christensen, 2015-2016.
    License: BSD 3-Clause
"""
from datetime import datetime
import multiprocessing
from netCDF4 import Dataset
import numpy as np

# local
from ..dataset import RAPIDDataset
from ..helper_functions import add_latlon_metadata, log
from ..utilities import partition


def generate_single_return_period(args):
    """
    This function calculates a single return period for a single reach
    """
    qout_file, return_period_file, rivid_index_list, step, num_years, \
        method, mp_lock = args

    skewvals = [-3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2,
                -1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0,
                1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
    kfac2 = [0.396, 0.384, 0.368, 0.351, 0.33, 0.307, 0.282, 0.254, 0.225,
             0.195, 0.164, 0.132, 0.099, 0.066, 0.033, 0, -0.033, -0.066,
             -0.099, -0.132, -0.164, -0.195, -0.225, -0.254, -0.282, -0.307,
             -0.33, -0.351, -0.368, -0.384, -0.396]
    kfac10 = [0.66, 0.702, 0.747, 0.795, 0.844, 0.895, 0.945, 0.994, 1.041,
              1.086, 1.128, 1.166, 1.2, 1.231, 1.258, 1.282, 1.301, 1.317,
              1.328, 1.336, 1.34, 1.34, 1.337, 1.329, 1.318, 1.302, 1.284,
              1.262, 1.238, 1.21, 1.18]
    kfac25 = [.666, .712, .764, .823, .888, .959, 1.035, 1.116, 1.198, 1.282,
              1.366, 1.448, 1.528, 1.606, 1.680, 1.751, 1.818, 1.880, 1.939,
              1.993, 2.043, 2.087, 2.128, 2.163, 2.193, 2.219, 2.240, 2.256,
              2.267, 2.275, 2.278]
    kfac50 = [0.666, 0.714, 0.768, 0.83, 0.9, 0.98, 1.069, 1.166, 1.27, 1.379,
              1.492, 1.606, 1.72, 1.834, 1.945, 2.054, 2.159, 2.261, 2.359,
              2.453, 2.542, 2.626, 2.706, 2.78, 2.848, 2.912, 2.97, 3.023,
              3.071, 3.114, 3.152]
    kfac100 = [0.667, 0.714, 0.769, 0.832, 0.905, 0.99, 1.087, 1.197, 1.318,
               1.499, 1.588, 1.733, 1.88, 2.029, 2.178, 2.326, 2.472, 2.615,
               2.755, 2.891, 3.022, 3.149, 3.271, 3.388, 3.499, 3.605, 3.705,
               3.8, 3.889, 3.973, 4.051]

    with RAPIDDataset(qout_file) as qout_nc_file:
        # get index of return period data
        if method == 'weibull':
            rp_index_20 = int((num_years + 1)/20.0)
            rp_index_10 = int((num_years + 1)/10.0)
            rp_index_2 = int((num_years + 1)/2.0)

        if method == 'weibull':
            return_20_array = np.zeros(len(rivid_index_list))
        elif method == 'gumble':
            return_100_array = np.zeros(len(rivid_index_list))
            return_50_array = np.zeros(len(rivid_index_list))
            return_20_array = np.zeros(len(rivid_index_list))
        elif method == 'log_pearson':
            return_100_array = np.zeros(len(rivid_index_list))
            return_50_array = np.zeros(len(rivid_index_list))
            return_25_array = np.zeros(len(rivid_index_list))
        return_10_array = np.zeros(len(rivid_index_list))
        return_2_array = np.zeros(len(rivid_index_list))
        max_flow_array = np.zeros(len(rivid_index_list))

        # iterate through rivids to generate return periods
        for iter_idx, rivid_index in enumerate(rivid_index_list):
            filtered_flow_data = qout_nc_file.get_qout_index(
                rivid_index,
                pd_filter="{0}D".format(step),
                filter_mode="max")
            sorted_flow_data = np.sort(filtered_flow_data)[:num_years:-1]
            max_flow = sorted_flow_data[0]
            if max_flow < 0.01:
                log("Return period data < 0.01 generated for rivid {0}"
                    .format(qout_nc_file.qout_nc.variables[
                        qout_nc_file.river_id_dimension][rivid_index]),
                    "WARNING")
            max_flow_array[iter_idx] = max_flow

            if method == 'weibull':
                return_20_array[iter_idx] = sorted_flow_data[rp_index_20]
                return_10_array[iter_idx] = sorted_flow_data[rp_index_10]
                return_2_array[iter_idx] = sorted_flow_data[rp_index_2]

            elif method == 'gumble':
                mean_flow = np.mean(filtered_flow_data)
                stddev = np.std(filtered_flow_data)
                return_100_array[iter_idx] = mean_flow + 3.14*stddev
                return_50_array[iter_idx] = mean_flow + 2.59*stddev
                return_20_array[iter_idx] = mean_flow + 1.87*stddev
                return_10_array[iter_idx] = mean_flow + 1.3*stddev
                return_2_array[iter_idx] = mean_flow - .164*stddev

            elif method == 'log_pearson':
                log_flow = np.log10(filtered_flow_data[filtered_flow_data > 0])
                if len(log_flow) <= 0:
                    continue
                mean_log_flow = np.mean(log_flow)
                std_log_flow = np.std(log_flow)
                log_flow_array = np.array(log_flow)
                skew = (num_years * (np.sum(
                    np.power((log_flow_array - mean_log_flow), 3)))) / \
                    ((num_years - 1) * (num_years - 2) * std_log_flow ** 3)
                k2 = np.interp(skew, skewvals, kfac2)
                k10 = np.interp(skew, skewvals, kfac10)
                k25 = np.interp(skew, skewvals, kfac25)
                k50 = np.interp(skew, skewvals, kfac50)
                k100 = np.interp(skew, skewvals, kfac100)
                return_100_array[iter_idx] = \
                    np.power(10, (mean_log_flow + k100*std_log_flow))
                return_50_array[iter_idx] = \
                    np.power(10, (mean_log_flow + k50*std_log_flow))
                return_25_array[iter_idx] = \
                    np.power(10, (mean_log_flow + k25*std_log_flow))
                return_10_array[iter_idx] = \
                    np.power(10, (mean_log_flow + k10*std_log_flow))
                return_2_array[iter_idx] = \
                    np.power(10, (mean_log_flow + k2*std_log_flow))

        mp_lock.acquire()
        return_period_nc = Dataset(return_period_file, 'a')
        return_period_nc.variables['max_flow'][rivid_index_list] = \
            max_flow_array
        if method == 'weibull':
            return_period_nc.variables['return_period_20'][
                rivid_index_list] = return_20_array
        elif method in 'gumble':
            return_period_nc.variables['return_period_100'][
                rivid_index_list] = return_100_array
            return_period_nc.variables['return_period_50'][
                rivid_index_list] = return_50_array
            return_period_nc.variables['return_period_20'][
                rivid_index_list] = return_20_array
        elif method == 'log_pearson':
            return_period_nc.variables['return_period_100'][
                rivid_index_list] = return_100_array
            return_period_nc.variables['return_period_50'][
                rivid_index_list] = return_50_array
            return_period_nc.variables['return_period_25'][
                rivid_index_list] = return_25_array
        return_period_nc.variables['return_period_10'][
            rivid_index_list] = return_10_array
        return_period_nc.variables['return_period_2'][
            rivid_index_list] = return_2_array
        return_period_nc.close()
        mp_lock.release()


def generate_return_periods(qout_file,
                            return_period_file,
                            num_cpus=multiprocessing.cpu_count(),
                            storm_duration_days=7,
                            method='weibull'):
    """
    Generate return period from RAPID Qout file
    """
    # get ERA Interim Data Analyzed
    with RAPIDDataset(qout_file) as qout_nc_file:
        print("Setting up Return Periods File ...")
        return_period_nc = Dataset(return_period_file, 'w')

        return_period_nc.createDimension('rivid', qout_nc_file.size_river_id)

        timeSeries_var = \
            return_period_nc.createVariable('rivid', 'i4', ('rivid',))
        timeSeries_var.long_name = (
            'unique identifier for each river reach')

        max_flow_var = \
            return_period_nc.createVariable('max_flow', 'f8', ('rivid',))
        max_flow_var.long_name = 'maximum streamflow'
        max_flow_var.units = 'm3/s'

        if method == 'weibull':

            return_period_20_var = \
                return_period_nc.createVariable('return_period_20',
                                                'f8', ('rivid',))
            return_period_20_var.long_name = '20 year return period flow'
            return_period_20_var.units = 'm3/s'

        if method == 'gumble':

            return_period_100_var = \
                return_period_nc.createVariable('return_period_100',
                                                'f8', ('rivid',))
            return_period_100_var.long_name = '100 year return period flow'
            return_period_100_var.units = 'm3/s'

            return_period_50_var = \
                return_period_nc.createVariable('return_period_50',
                                                'f8', ('rivid',))
            return_period_50_var.long_name = '50 year return period flow'
            return_period_50_var.units = 'm3/s'

            return_period_20_var = \
                return_period_nc.createVariable('return_period_20',
                                                'f8', ('rivid',))
            return_period_20_var.long_name = '20 year return period flow'
            return_period_20_var.units = 'm3/s'

        if method == 'log_pearson':

            return_period_100_var = \
                return_period_nc.createVariable('return_period_100',
                                                'f8', ('rivid',))
            return_period_100_var.long_name = '100 year return period flow'
            return_period_100_var.units = 'm3/s'

            return_period_50_var = \
                return_period_nc.createVariable('return_period_50',
                                                'f8', ('rivid',))
            return_period_50_var.long_name = '50 year return period flow'
            return_period_50_var.units = 'm3/s'

            return_period_25_var = \
                return_period_nc.createVariable('return_period_25',
                                                'f8', ('rivid',))
            return_period_25_var.long_name = '25 year return period flow'
            return_period_25_var.units = 'm3/s'

        return_period_10_var = \
            return_period_nc.createVariable('return_period_10',
                                            'f8', ('rivid',))
        return_period_10_var.long_name = '10 year return period flow'
        return_period_10_var.units = 'm3/s'

        return_period_2_var = \
            return_period_nc.createVariable('return_period_2',
                                            'f8', ('rivid',))
        return_period_2_var.long_name = '2 year return period flow'
        return_period_2_var.units = 'm3/s'

        lat_var = return_period_nc.createVariable('lat', 'f8', ('rivid',),
                                                  fill_value=-9999.0)

        lon_var = return_period_nc.createVariable('lon', 'f8', ('rivid',),
                                                  fill_value=-9999.0)

        add_latlon_metadata(lat_var, lon_var)

        return_period_nc.variables['lat'][:] = \
            qout_nc_file.qout_nc.variables['lat'][:]
        return_period_nc.variables['lon'][:] = \
            qout_nc_file.qout_nc.variables['lon'][:]

        river_id_list = qout_nc_file.get_river_id_array()
        return_period_nc.variables['rivid'][:] = river_id_list

        return_period_nc.return_period_method = method

        return_period_nc.close()

        time_array = qout_nc_file.get_time_array()

    log("Extracting Data and Generating Return Periods ...")
    num_years = int((datetime.utcfromtimestamp(time_array[-1]) -
                     datetime.utcfromtimestamp(time_array[0])).days/365.2425)
    time_steps_per_day = (24 * 3600) / float(
        (datetime.utcfromtimestamp(time_array[1]) -
         datetime.utcfromtimestamp(time_array[0])).total_seconds())
    step = max(1, int(time_steps_per_day * storm_duration_days))

    # generate multiprocessing jobs
    # pylint: disable=no-member
    mp_lock = multiprocessing.Manager().Lock()
    job_combinations = []
    partition_index_list = partition(river_id_list, num_cpus*2)[1]
    for sub_partition_index_list in partition_index_list:
        # pylint: disable=len-as-condition
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
