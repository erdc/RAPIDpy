# -*- coding: utf-8 -*-
"""
    generate_seasonal_averages.py
    RAPIDpy

    Created by: Alan D. Snow, 2016.
    License: BSD 3-Clause
"""
from calendar import isleap
import multiprocessing
from time import gmtime

from netCDF4 import Dataset
import numpy as np

from ..dataset import RAPIDDataset
from ..helper_functions import add_latlon_metadata


def generate_single_seasonal_average(args):
    """
    This function calculates the seasonal average for a single day of the year
    for all river segments
    """
    qout_file = args[0]
    seasonal_average_file = args[1]
    day_of_year = args[2]
    mp_lock = args[3]

    min_day = day_of_year - 3
    max_day = day_of_year + 3

    with RAPIDDataset(qout_file) as qout_nc_file:
        time_indices = []
        for idx, t in enumerate(qout_nc_file.get_time_array()):
            var_time = gmtime(t)
            compare_yday = var_time.tm_yday
            # move day back one past because of leap year adds
            # a day after feb 29 (day 60)
            if isleap(var_time.tm_year) and compare_yday > 60:
                compare_yday -= 1
            # check if date within range of season
            if max_day > compare_yday >= min_day:
                time_indices.append(idx)

        if not time_indices:
            raise IndexError("No time steps found within range ...")

        streamflow_array = qout_nc_file.get_qout(time_index_array=time_indices)

    avg_streamflow_array = np.mean(streamflow_array, axis=1)
    std_streamflow_array = np.std(streamflow_array, axis=1)
    max_streamflow_array = np.amax(streamflow_array, axis=1)
    min_streamflow_array = np.min(streamflow_array, axis=1)

    mp_lock.acquire()
    seasonal_avg_nc = Dataset(seasonal_average_file, 'a')
    seasonal_avg_nc.variables['average_flow'][:, day_of_year-1] = \
        avg_streamflow_array
    seasonal_avg_nc.variables['std_dev_flow'][:, day_of_year-1] = \
        std_streamflow_array
    seasonal_avg_nc.variables['max_flow'][:, day_of_year-1] = \
        max_streamflow_array
    seasonal_avg_nc.variables['min_flow'][:, day_of_year-1] = \
        min_streamflow_array
    seasonal_avg_nc.close()
    mp_lock.release()


def generate_seasonal_averages(qout_file, seasonal_average_file,
                               num_cpus=multiprocessing.cpu_count()):
    """
    This function loops through a CF compliant rapid streamflow
    file to produce a netCDF file with a seasonal average for
    365 days a year
    """
    with RAPIDDataset(qout_file) as qout_nc_file:
        print("Generating seasonal average file ...")
        seasonal_avg_nc = Dataset(seasonal_average_file, 'w')

        seasonal_avg_nc.createDimension('rivid', qout_nc_file.size_river_id)
        seasonal_avg_nc.createDimension('day_of_year', 365)

        time_series_var = seasonal_avg_nc.createVariable('rivid', 'i4',
                                                         ('rivid',))
        time_series_var.long_name = (
            'unique identifier for each river reach')

        average_flow_var = \
            seasonal_avg_nc.createVariable('average_flow', 'f8',
                                           ('rivid', 'day_of_year'))
        average_flow_var.long_name = 'seasonal average streamflow'
        average_flow_var.units = 'm3/s'

        std_dev_flow_var = \
            seasonal_avg_nc.createVariable('std_dev_flow', 'f8',
                                           ('rivid', 'day_of_year'))
        std_dev_flow_var.long_name = 'seasonal std. dev. streamflow'
        std_dev_flow_var.units = 'm3/s'

        std_dev_flow_var = \
            seasonal_avg_nc.createVariable('max_flow', 'f8',
                                           ('rivid', 'day_of_year'))
        std_dev_flow_var.long_name = 'seasonal max streamflow'
        std_dev_flow_var.units = 'm3/s'

        std_dev_flow_var = \
            seasonal_avg_nc.createVariable('min_flow', 'f8',
                                           ('rivid', 'day_of_year'))
        std_dev_flow_var.long_name = 'seasonal min streamflow'
        std_dev_flow_var.units = 'm3/s'

        lat_var = seasonal_avg_nc.createVariable('lat', 'f8', ('rivid',),
                                                 fill_value=-9999.0)

        lon_var = seasonal_avg_nc.createVariable('lon', 'f8', ('rivid',),
                                                 fill_value=-9999.0)
        add_latlon_metadata(lat_var, lon_var)

        seasonal_avg_nc.variables['lat'][:] = \
            qout_nc_file.qout_nc.variables['lat'][:]
        seasonal_avg_nc.variables['lon'][:] = \
            qout_nc_file.qout_nc.variables['lon'][:]

        river_id_list = qout_nc_file.get_river_id_array()
        seasonal_avg_nc.variables['rivid'][:] = river_id_list
        seasonal_avg_nc.close()

    # generate multiprocessing jobs
    mp_lock = multiprocessing.Manager().Lock()  # pylint: disable=no-member
    job_combinations = []
    for day_of_year in range(1, 366):
        job_combinations.append((qout_file,
                                 seasonal_average_file,
                                 day_of_year,
                                 mp_lock
                                 ))

    pool = multiprocessing.Pool(num_cpus)
    pool.map(generate_single_seasonal_average,
             job_combinations)
    pool.close()
    pool.join()
