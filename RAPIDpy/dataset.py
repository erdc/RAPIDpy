# -*- coding: utf-8 -*-
##
##  dataset.py
##  RAPIDpy
##
##  Created by Alan D Snow.
##  Copyright © 2016 Alan D Snow. All rights reserved.
##  BSD 3-Clause

from csv import writer as csv_writer
import datetime
from netCDF4 import Dataset, num2date
import numpy as np
from numpy.ma import masked
from past.builtins import xrange
from pytz import utc
import time

from .helper_functions import log, open_csv
    
#------------------------------------------------------------------------------
#Helper Function
#------------------------------------------------------------------------------
def compare_qout_files(dataset1_path, dataset2_path, Qout_var="Qout"):
    """
    This function compares the output of RAPID Qout and tells you where they are different.
    """
    qout_same = False
    
    d1 = RAPIDDataset(dataset1_path)
    d2 = RAPIDDataset(dataset2_path)

    if len(d1.get_river_id_array()) != len(d2.get_river_id_array()):
        log("Length of COMID/rivid input not the same.",
            "ERROR")

    if not (d1.get_river_id_array() == d2.get_river_id_array()).all():
        log("COMID/rivid order is different in each dataset. Reordering data for comparison.",
            "WARNING")
        
        d2_reordered_river_index_list = []
        for comid in d1.get_river_id_array():
            d2_reordered_river_index_list.append(np.where(d2.get_river_id_array()==comid)[0][0])
        d2_reordered_qout = d2.get_qout_index(d2_reordered_river_index_list)
    else:
        d2_reordered_qout = d2.get_qout()
        
    #get where the files are different
    d1_qout = d1.get_qout()
    where_diff = np.where(d1_qout != d2_reordered_qout)
    un_where_diff = np.unique(where_diff[0])
    
    #if different, check to see how different
    if un_where_diff.any():
        decimal_test = 7
        while decimal_test > 0:
            try:
                np.testing.assert_almost_equal(d1_qout,
                                               d2_reordered_qout, 
                                               decimal=decimal_test)
                log("ALMOST EQUAL to {0} decimal places.".format(decimal_test),
                    "INFO")
                qout_same = True
                decimal_test=-1
            except AssertionError as ex:
                if decimal_test <= 1:
                    print(ex)
                decimal_test-=1
                pass
        log("Number of different timeseries: {0}".format(len(un_where_diff)),
            "INFO")
        log("COMID idexes where different: {0}".format(un_where_diff),
            "INFO")
        log("COMID idexes where different: {0}".format(un_where_diff),
            "INFO")
        index = un_where_diff[0]
        log("Dataset 1 example. COMID index: {0}".format(d1.get_qout_index(index)),
            "INFO")
        log("Dataset 2 example. COMID index: {0}".format(d2_reordered_qout[index, :]),
            "INFO")
    
    else:
        qout_same = True
        log("Output Qout data is the same.",
            "INFO")

    d1.close()
    d2.close()
    return qout_same
    
#------------------------------------------------------------------------------
#Main Dataset Manager Class
#------------------------------------------------------------------------------
class RAPIDDataset(object):
    """
    This class is designed to access data from the RAPID Qout
    NetCDF file.
    
    Attributes:
        filename(str): Path to the RAPID Qout NetCDF file.
        river_id_dimension(Optional[str]): Name of the river ID dimension. Default is to search through a standard list.
        river_id_variable(Optional[str]): Name of the river ID variable. Default is to search through a standard list.
        streamflow_variable(Optional[str]): Name of the streamflow varaible. Default is to search through a standard list.

    Example::
    
        from RAPIDpy import RAPIDDataset

        path_to_rapid_qout = '/path/to/Qout.nc'
        with RAPIDDataset(path_to_rapid_qout) as qout_nc:
            #USE FUNCTIONS TO ACCESS DATA HERE
                         
    """

    def __init__(self, filename, 
                 river_id_dimension="", 
                 river_id_variable="", 
                 streamflow_variable=""):
        """
        Initialize the class with variables given by the user
        """
        self.qout_nc = Dataset(filename, mode='r')
        
        #determine river ID dimension
        self.river_id_dimension = river_id_dimension
        if not river_id_dimension:
            if 'rivid' in self.qout_nc.dimensions:
                self.river_id_dimension = 'rivid'
            elif 'COMID' in self.qout_nc.dimensions:
                self.river_id_dimension = 'COMID'
            elif 'station' in self.qout_nc.dimensions:
                self.river_id_dimension = 'station'
            elif 'DrainLnID' in self.qout_nc.dimensions:
                self.river_id_dimension = 'DrainLnID'
            elif 'FEATUREID' in self.qout_nc.dimensions:
                self.river_id_dimension = 'FEATUREID'
            else:
                raise IndexError('ERROR: Could not find river ID dimension.')
        elif river_id_dimension not in self.qout_nc.dimensions:
            raise IndexError('ERROR: Could not find river ID dimension: {0}.'.format(river_id_dimension))
            
        self.size_river_id = len(self.qout_nc.dimensions[self.river_id_dimension])
        
        variable_keys = self.qout_nc.variables.keys()

        #determin streamflow variable
        self.q_var_name = streamflow_variable
        if not streamflow_variable:
            if 'Qout' in variable_keys:
                self.q_var_name = 'Qout'
            elif 'streamflow' in variable_keys:
                self.q_var_name = 'streamflow'
            elif 'm3_riv' in variable_keys:
                self.q_var_name = 'm3_riv'
            else:
                raise IndexError('ERROR: Could not find flow variable. Looked for Qout, streamflow, and m3_riv.')
        elif not streamflow_variable in variable_keys:
            raise IndexError('ERROR: Could not find flow variable. Looked for {0}.'.format(streamflow_variable))

        self.size_q_var = len(self.qout_nc.variables[self.q_var_name])

        #determine time dimension
        if 'time' in self.qout_nc.dimensions:
            self.size_time = len(self.qout_nc.dimensions['time'])
        elif 'Time' in self.qout_nc.dimensions:
            self.size_time = len(self.qout_nc.dimensions['Time'])
        else:
            raise IndexError('ERROR: Could not find time dimension.')

        #determine river ID variable
        self.river_id_variable = river_id_variable
        if not river_id_variable:
            if 'rivid' in variable_keys:
                self.river_id_variable = 'rivid'
            elif 'COMID' in variable_keys:
                self.river_id_variable = 'COMID'
            elif 'station_id' in variable_keys:
                self.river_id_variable = 'station_id'
            elif 'DrainLnID' in variable_keys:
                self.river_id_variable = 'DrainLnID'
            elif 'FEATUREID' in variable_keys:
                self.river_id_variable = 'FEATUREID'
            else:
                print('WARNING: Could not find river ID variable in {0}.'.format(variable_keys))
        elif river_id_variable not in variable_keys:
            print('WARNING: Could not find river ID variable: {0}.'.format(river_id_variable))

    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        self.qout_nc.close()

    def is_time_variable_valid(self):
        """
        This function returns whether or not the time variable
        is valid.
        
        Returns:
            boolean: True if the time variable is valid, otherwise false.
        
        Example::
        
            from RAPIDpy import RAPIDDataset
    
            path_to_rapid_qout = '/path/to/Qout.nc'
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                if qout_nc.is_time_variable_valid():
                    #DO WORK HERE
        """
        time_var_valid = False
        if 'time' in self.qout_nc.variables.keys():
            if len(self.qout_nc.dimensions['time'])>0:
                if not (self.qout_nc.variables['time'][:] == masked).any():
                    try:
                        timestep = (datetime.datetime.utcfromtimestamp(self.qout_nc.variables['time'][1])-
                                    datetime.datetime.utcfromtimestamp(self.qout_nc.variables['time'][0])).total_seconds()
                        if timestep > 0:
                            time_var_valid = True
                    except ValueError:
                        pass
    
        return time_var_valid
    
    def get_time_array(self, 
                       datetime_simulation_start=None,
                       simulation_time_step_seconds=None,
                       return_datetime=False,
                       time_index_array=None):
        """
        This method extracts or generates an array of time. The new version of RAPID output has the time array stored. 
        However, the old version requires the user to know when the simulation began and the time step of the output.
        
        Parameters:
            datetime_simulation_start(Optional[datetime]): This is a datetime object with the date of the simulation start time.
            simulation_time_step_seconds(Optional[integer]): This is the time step of the simulation output in seconds.
            return_datetime(Optional[boolean]): If true, it converts the data to a list of datetime objects. Default is False.
            time_index_array(Optional[list or np.array]): This is used to extract the datetime vales. This can be from the *get_time_index_range* function.
            
        Returns:
            list: An array of integers representing seconds since Jan 1, 1970 UTC or datetime objects if return_datetime is set to True.
        
        This example demonstrates how to retrieve or generate a time array to go
        along with your RAPID streamflow series::
        
            from RAPIDpy import RAPIDDataset
    
            path_to_rapid_qout = '/path/to/Qout.nc'
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                #CF-Compliant version
                if qout_nc.is_time_variable_valid():
                    time_array = qout_nc.get_time_array()
                    
                    #or, to get datetime array
                    time_datetime = qout_nc.get_time_array(return_datetime=True)
                #Original version
                else:
                    time_array = qout_nc.get_time_array(datetime_simulation_start=datetime(1980, 1, 1),
                                                        simulation_time_step_seconds=3*3600)
                                                        
                    #or, to get datetime array
                    time_datetime = qout_nc.get_time_array(datetime_simulation_start=datetime(1980, 1, 1),
                                                           simulation_time_step_seconds=3*3600,
                                                           return_datetime=True)
                    
        """
        time_array = []

        epoch = datetime.datetime(1970,1,1, tzinfo=utc)
        time_units = "seconds since {0}".format(epoch)
        
        #CF-1.6 compliant file
        if self.is_time_variable_valid():
            time_array = self.qout_nc.variables['time'][:]
            if self.qout_nc.variables['time'].units:
                time_units = self.qout_nc.variables['time'].units
            
        #Original Qout file
        elif datetime_simulation_start is not None and simulation_time_step_seconds is not None:
            initial_time_seconds = (datetime_simulation_start.replace(tzinfo=utc)-
                                    epoch).total_seconds()+simulation_time_step_seconds
            final_time_seconds = initial_time_seconds + self.size_time*simulation_time_step_seconds
            time_array = np.arange(initial_time_seconds, final_time_seconds, simulation_time_step_seconds)
        else:
            raise Exception("ERROR: This file does not contain the time variable."
                            " To get time array, add datetime_simulation_start"
                            " and simulation_time_step_seconds")
        
        if time_index_array is not None:
            time_array = time_array[time_index_array]

        if return_datetime:
            return num2date(time_array, time_units)
        
        return time_array

    def get_time_index_range(self, date_search_start=None,
                             date_search_end=None,
                             time_index_start=None,
                             time_index_end=None,
                             time_index=None):
        """
        Generates a time index range based on time bounds given. This is useful for subset data extraction.
        
        Parameters:
            date_search_start(Optional[datetime]): This is a datetime object with the date of the minimum date for starting.
            date_search_end(Optional[datetime]): This is a datetime object with the date of the maximum date for ending.
            time_index_start(Optional[int]): This is the index of the start of the time array subset. Useful for the old file version.
            time_index_end(Optional[int]): This is the index of the end of the time array subset. Useful for the old file version.
            time_index(Optional[int]): This is the index of time to return in the case that your code only wants one index. Used internally.
            
        Returns:
            index_array: This is an array used to extract a subset of data.
        
        Example::
        
            from RAPIDpy import RAPIDDataset
    
            path_to_rapid_qout = '/path/to/Qout.nc'
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                #CF-Compliant version
                if qout_nc.is_time_variable_valid():
                    time_index_range = qout_nc.get_time_index_range(date_search_start=datetime(1980, 1, 1),
                                                                    date_search_end=datetime(1980, 12, 11))
                #Original version
                else:
                    time_index_range = qout_nc.get_time_index_range(time_index_start=22,
                                                                    time_index_end=40)
                    
        """
        #get the range of time based on datetime range
        time_range = None
        if self.is_time_variable_valid() and (date_search_start is not None or date_search_end is not None):
            print("Determining time range ({0} to {1})...".format(date_search_start, date_search_end))
            time_array = self.qout_nc.variables['time'][:]
            if date_search_start is not None:
                seconds_start = (date_search_start-datetime.datetime(1970,1,1)).total_seconds()
                time_range = np.where(time_array>=seconds_start)[0]
            
            if date_search_end is not None:
                seconds_end = (date_search_end-datetime.datetime(1970,1,1)).total_seconds()
                if time_range is not None:
                    time_range = np.intersect1d(time_range, np.where(time_array<=seconds_end)[0])
                else:
                    time_range = np.where(time_array<=seconds_end)[0]
    
        #get the range of time based on time index range
        elif time_index_start is not None or time_index_end is not None:
            if time_index_start is None:
                time_index_start = 0
            if time_index_end is None:
                time_index_end = self.size_time
            time_range = range(time_index_start,time_index_end)

        #get only one time step
        elif time_index is not None:
            time_range = time_index
        #return all
        else:
            time_range = range(self.size_time)
        
        return time_range

    def get_daily_time_index_array(self, time_index_range=None):
        """
        Returns an array of the first index of each day in the time array
        """
        idx = 0
        if time_index_range is None:
            datetime_array = self.get_time_array(return_datetime=True)
        else:
            idx = time_index_range[0]
            datetime_array = self.get_time_array(time_index_array=time_index_range,
                                                 return_datetime=True)
       
        current_day = datetime_array[0]
        daily_time_index_array = [idx]
        for var_time in datetime_array:
            if current_day.day != var_time.day:
                 daily_time_index_array.append(idx)
            current_day = var_time
            idx += 1
        return daily_time_index_array

    def get_river_id_array(self):
        """
        This method returns the river ID array for this file.
        
        Returns:
            numpy.array: An array of the river ID's
        
        Example::
        
            from RAPIDpy import RAPIDDataset
    
            path_to_rapid_qout = '/path/to/Qout.nc'
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                river_ids = qout_nc.get_river_id_array()
                
        """ 
        return self.qout_nc.variables[self.river_id_variable][:]
    
    def get_river_index(self, river_id):
        """
        This method retrieves the river index in the netCDF
        dataset corresponding to the river ID.

        Returns:
            int: The index of the river ID's in the file
        
        Example::
        
            from RAPIDpy import RAPIDDataset
    
            path_to_rapid_qout = '/path/to/Qout.nc'
            river_id = 53458
            
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                river_index = qout_nc.get_river_index(river_id)
        """
        try:
            return np.where(self.get_river_id_array()==river_id)[0][0]
        except IndexError:
            raise IndexError("ERROR: River ID", river_id, "not found in dataset ...")

    def get_subset_riverid_index_list(self, river_id_list):
        """
        Gets the subset riverid_list from the netcdf file
        Optional returns include the list of valid river ids in the dataset
        as well as a list of missing rive rids
        """
        netcdf_river_indices_list = []
        valid_river_ids = []
        missing_river_ids = []
        for river_id in river_id_list:
            #get where streamids are in netcdf file
            try:
                netcdf_river_indices_list.append(self.get_river_index(river_id))
                valid_river_ids.append(river_id)
            except IndexError:
                print("WARNING: ReachID {0} not found in netCDF dataset. Skipping ...".format(river_id))
                missing_river_ids.append(river_id)
                pass
        
        np_valid_river_indices_list = np.array(netcdf_river_indices_list)
        np_valid_river_ids = np.array(valid_river_ids)
        sorted_indexes = np.argsort(np_valid_river_indices_list)
        
        return np_valid_river_indices_list[sorted_indexes], \
               np_valid_river_ids[sorted_indexes], \
               np.array(missing_river_ids)

    

    def get_qout(self, river_id_array=None,
                 date_search_start=None,
                 date_search_end=None,
                 time_index_start=None,
                 time_index_end=None,
                 time_index=None,
                 time_index_array=None):
        """
        This method extracts streamflow data by a single river ID or by a river ID array. 
        It has options to extract by date or by date index.
        
        Parameters:
            river_id_array(Optional[list or int]): A single river ID or an array of river IDs.
            date_search_start(Optional[datetime]): This is a datetime object with the date of the minimum date for starting.
            date_search_end(Optional[datetime]): This is a datetime object with the date of the maximum date for ending.
            time_index_start(Optional[int]): This is the index of the start of the time array subset. Useful for the old file version.
            time_index_end(Optional[int]): This is the index of the end of the time array subset. Useful for the old file version.
            time_index(Optional[int]): This is the index of time to return in the case that your code only wants one index. Used internally.
            time_index_array(Optional[list or np.array]): This is used to extract the vales only for particular dates. This can be from the *get_time_index_range* function.
            
        Returns:
            numpy.array: This is a 1D or 2D array or a single value depending on your input search. 
        
        This example demonstrates how to retrieve the streamflow associated with
        the reach you are interested in::

            from RAPIDpy import RAPIDDataset
    
            path_to_rapid_qout = '/path/to/Qout.nc'
            river_id = 500
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                streamflow_array = qout_nc.get_qout(river_id)

        This example demonstrates how to retrieve the streamflow within a date
        range associated with the reach you are interested in::

            from RAPIDpy import RAPIDDataset
    
            path_to_rapid_qout = '/path/to/Qout.nc'
            river_id = 500
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:

                #CF-Compliant version
                if qout_nc.is_time_variable_valid():
                    streamflow_array = qout_nc.get_qout(river_id,
                                                        date_search_start=datetime(1985,1,1),
                                                        date_search_end=datetime(1985,2,4))
                #Original version
                else:
                    streamflow_array = qout_nc.get_qout(river_id,
                                                        time_index_start=20,
                                                        time_index_end=25)
                    
        """
        #get indices of where the streamflow data is
        riverid_index_list_subset = None
        if river_id_array is not None:
            if hasattr(river_id_array, "__len__") \
            and not isinstance(river_id_array, str):
                river_id_array = [river_id_array]
            riverid_index_list_subset = self.get_subset_riverid_index_list(river_id_array)[0]

        return self.get_qout_index(riverid_index_list_subset,
                                   date_search_start,
                                   date_search_end,
                                   time_index_start,
                                   time_index_end,
                                   time_index,
                                   time_index_array)
                       


    def get_qout_index(self, river_index_array=None,
                       date_search_start=None,
                       date_search_end=None,
                       time_index_start=None,
                       time_index_end=None,
                       time_index=None,
                       time_index_array=None):
        """
        This method extracts streamflow data by river index
        It allows for extracting single or multiple river streamflow arrays
        It has options to extract by date or by date index
        """
        if river_index_array is not None:
            if isinstance(river_index_array, list) or type(river_index_array).__module__ == np.array:
                if len(river_index_array) == 1:
                    river_index_array = river_index_array[0]
        if time_index_array is None:
            time_index_array = self.get_time_index_range(date_search_start,
                                                         date_search_end,
                                                         time_index_start,
                                                         time_index_end,
                                                         time_index)

        qout_variable = self.qout_nc.variables[self.q_var_name]
        qout_dimensions = qout_variable.dimensions
        streamflow_array = []
        if qout_dimensions[0].lower() == 'time' and qout_dimensions[1].lower() == self.river_id_dimension.lower():
            if time_index_array is not None and river_index_array is not None:
                streamflow_array = qout_variable[time_index_array,river_index_array].transpose()
            elif time_index_array is not None:
                streamflow_array = qout_variable[time_index_array,:].transpose()
            elif river_index_array is not None:
                streamflow_array = qout_variable[:,river_index_array].transpose()
            else:
                streamflow_array = qout_variable[:].transpose()
        elif qout_dimensions[1].lower() == 'time' and qout_dimensions[0].lower() == self.river_id_dimension.lower():
            if time_index_array is not None and river_index_array is not None:
                streamflow_array = qout_variable[river_index_array, time_index_array]
            elif time_index_array is not None:
                streamflow_array = qout_variable[:, time_index_array]
            elif river_index_array is not None:
                streamflow_array = qout_variable[river_index_array, :]
            else:
                streamflow_array = qout_variable[:]
        else:
            raise Exception( "Invalid RAPID Qout file dimensions ...")
        return streamflow_array

    def get_daily_qout_index(self, 
                             river_index_array, 
                             daily_time_index_array=None, 
                             steps_per_group=1, 
                             mode="mean"):
        """
        Gets the daily time series from RAPID output from river ID index.
        
        Parameters:
            river_index_array(list or int): A single river index or an array of river indices.
            daily_time_index_array(Optional[list of time indices]): This is a list of indices from the get_daily_time_index_array function.
            steps_per_group(Optional[int]): This is how many time steps per day. This is for the old version of RAPID Qout.
            mode(Optional[str]): You can get the daily average "mean" or the maximum "max".
            
        Returns:
            numpy.array: This is a 1D or 2D array or a single value depending on your input search. 
        
        This example demonstrates how to get daily streamflow averages as an
        array::

        Example::
        
            from RAPIDpy import RAPIDDataset
    
            path_to_rapid_qout = '/path/to/Qout.nc'
            river_id = 500
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                river_index = qout_nc.get_river_index(river_id)

                #CF-Compliant version
                if qout_nc.is_time_variable_valid():
                    streamflow_array = qout_nc.get_daily_qout_index(river_index)
                #Original version
                else:
                    streamflow_array = qout_nc.get_daily_qout_index(river_index,
                                                                    steps_per_group=8, #average 8 timesteps together for 1 day
                                                                    )
        """
        axis = None
        if mode=="mean":
            calc = np.mean
        elif mode=="max":
    	    calc = np.amax
        else:
    	    raise Exception("Invalid calc mode ...")

        if self.is_time_variable_valid() and steps_per_group<=1:
            qout_arr = self.get_qout_index(river_index_array)
            qout_arr_dim_size = len(qout_arr.shape)
            
            if not daily_time_index_array:
                daily_time_index_array = self.get_daily_time_index_array()
            
            last_possible_index = daily_time_index_array[-1]
            #IF NOT END OF ALL TIME, MAKE SURE ENTIRE DAY IS CAPTURED
            if last_possible_index < self.size_time-1:
                remaining_time_index_array = self.get_daily_time_index_array(range(last_possible_index, self.size_time))
                if len(remaining_time_index_array) > 1:
                    last_possible_index = remaining_time_index_array[1]
                else:
                    last_possible_index = -1
                
            len_daily_time_array = len(daily_time_index_array)
            
            
            if qout_arr_dim_size > 1:
                daily_qout = np.zeros((qout_arr.shape[0], len_daily_time_array))
            else:
                daily_qout = np.zeros(len_daily_time_array)
                
            for idx in xrange(len_daily_time_array):
                time_index_start = daily_time_index_array[idx]
                if idx+1 < len_daily_time_array:
                    next_time_index = daily_time_index_array[idx+1]
                    if qout_arr_dim_size > 1:
                        daily_qout[:,idx] = calc(qout_arr[:,time_index_start:next_time_index], axis=1)
                    else:
                        daily_qout[idx] = calc(qout_arr[time_index_start:next_time_index])
                elif idx+1 == len_daily_time_array:
                    if time_index_start < self.size_time - 1:
                        if qout_arr_dim_size > 1:
                            daily_qout[:,idx] =  calc(qout_arr[:,time_index_start:last_possible_index], axis=1)
                        else:
                            daily_qout[idx] =  calc(qout_arr[time_index_start:last_possible_index])
                    else:
                        if qout_arr_dim_size > 1:
                            daily_qout[:,idx] =  qout_arr[:,time_index_start]
                        else:
                            daily_qout[idx] =  qout_arr[time_index_start]
            return daily_qout
            
        elif steps_per_group > 1:
            flow_data = self.get_qout_index(river_index_array)
            qout_arr_dim_size = len(flow_data.shape)
            axis = None
            if qout_arr_dim_size > 1:
                axis = 1
            daily_qout = []
            for step_index in xrange(0, len(flow_data), steps_per_group):
                if qout_arr_dim_size > 1:
                    flows_slice = flow_data[:,step_index:step_index + steps_per_group]
                else:
                    flows_slice = flow_data[step_index:step_index + steps_per_group]
                daily_qout.append(calc(flows_slice, axis=axis))
            return np.array(daily_qout, np.float32)
        else:
            raise Exception("Must have steps_per_group set to a value greater than one "
                            "due to non CF-Compliant Qout file ...")

    def get_daily_qout(self, 
                       river_id, 
                       daily_time_index_array=None,
                       steps_per_group=1, 
                       mode="mean"):
        """
        Retrieves the daily qout for a river ID from RAPID time series
        
        Parameters:
            river_index_array(list or int): A single river index or an array of river indices.
            daily_time_index_array(Optional[list of time indices]): This is a list of indices from the get_daily_time_index_array function.
            steps_per_group(Optional[int]): This is how many time steps per day. This is for the old version of RAPID Qout.
            mode(Optional[str]): You can get the daily average "mean" or the maximum "max". Defauls is "mean".
            
        Returns:
            numpy.array: This is a 1D or 2D array or a single value depending on your input search. 
        
        This example demonstrates how to get daily streamflow averages as an
        array::
        
            from RAPIDpy import RAPIDDataset
    
            path_to_rapid_qout = '/path/to/Qout.nc'
            river_id = 500
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                #CF-Compliant version
                if qout_nc.is_time_variable_valid():
                    streamflow_array = qout_nc.get_daily_qout(river_id)
                #Original version
                else:
                    #average 8 timesteps together for 1 day
                    streamflow_array = qout_nc.get_daily_qout(river_id,
                                                              steps_per_group=8, 
                                                              )
    
        """
        self.get_daily_qout_index(self.get_river_index(river_id),
                                  daily_time_index_array,
                                  steps_per_group, mode)

    def get_seasonal_monthly_average(self,river_id_array,
                                     month):
        """
        This function loops through a CF compliant rapid streamflow
        file to produce estimates for current streamflow based on
        the seasonal average over the data within the historical streamflow
        file.
        """
        if not self.is_time_variable_valid():
            raise Exception("ERROR: File must be CF 1.6 compliant with time dimension ...")

        time_indices = []
        for idx, var_time in enumerate(self.get_time_array(return_datetime=True)):
            if var_time.month == month:
                time_indices.append(idx)

        if not time_indices:
            raise Exception("ERROR: No time steps found within range ...")
        
        print("Extracting data ...")
        return np.mean(self.get_qout(river_id_array, time_index_array=time_indices), axis=1)


    def write_flows_to_csv(self, path_to_output_file,
                           river_index=None, 
                           river_id=None,
                           date_search_start=None,
                           date_search_end=None,
                           daily=False,
                           mode="mean"):
        """
        Write out RAPID output to CSV file.
        
        .. note:: Need either *reach\_id* or *reach\_index* parameter, but either can be used.
        
        Parameters:
            path_to_output_file(str): Path to the output csv file.
            river_index(Optional[datetime]): This is the index of the river in the file you want the streamflow for.
            river_id(Optional[datetime]): This is the river ID that you want the streamflow for.
            date_search_start(Optional[datetime]): This is a datetime object with the date of the minimum date for starting.
            date_search_end(Optional[datetime]): This is a datetime object with the date of the maximum date for ending.
            daily(Optional[boolean]): If True and the file is CF-Compliant, write out daily flows.
            mode(Optional[str]): You can get the daily average "mean" or the maximum "max". Defauls is "mean".

        Example writing entire time series to file:
        
        .. code:: python
        
            from RAPIDpy import RAPIDDataset
    
            river_id = 3624735
            path_to_rapid_qout = '/path/to/Qout.nc'

            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                #for writing entire time series to file
                qout_nc.write_flows_to_csv('/timeseries/Qout_3624735.csv', 
                                           river_id=river_id,
                                           )
                                      
                river_index = qout_nc.get_river_index(river_id)
                
                #if file is CF compliant, you can write out daily average
                qout_nc.write_flows_to_csv('/timeseries/Qout_daily.csv',
                                           river_index=river_index,
                                           daily=True,
                                           )
                                           
        Example writing entire time series as daily average to file:
        
        .. code:: python
        
            from RAPIDpy import RAPIDDataset
    
            river_id = 3624735
            path_to_rapid_qout = '/path/to/Qout.nc'

            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                #NOTE: Getting the river index is not necessary
                #this is just an example of how to use this                                      
                river_index = qout_nc.get_river_index(river_id)
                
                #if file is CF compliant, you can write out daily average
                qout_nc.write_flows_to_csv('/timeseries/Qout_daily.csv',
                                           river_index=river_index,
                                           daily=True,
                                           )

        Example writing entire time series as daily average to file:
        
        .. code:: python
        
            from datetime import datetime
            from RAPIDpy import RAPIDDataset
    
            river_id = 3624735
            path_to_rapid_qout = '/path/to/Qout.nc'

            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                #if file is CF compliant, you can filter by date
                qout_nc.write_flows_to_csv('/timeseries/Qout_daily_date_filter.csv',
                                           river_index=river_index,
                                           daily=True,
                                           date_search_start=datetime(2002, 8, 31),
                                           date_search_end=datetime(2002, 9, 15),
                                           mode="max"
                                           )        
        """
        if river_id != None:
            river_index = self.get_river_index(river_id)
        elif river_id == None and river_index == None:
            raise Exception("ERROR: Need reach id or reach index ...")

        #analyze and write
        if self.is_time_variable_valid():
            time_index_range = self.get_time_index_range(date_search_start=date_search_start,
                                                         date_search_end=date_search_end)
            with open_csv(path_to_output_file, 'w') as outcsv:
                writer = csv_writer(outcsv)
                if daily:
                    daily_time_index_array = self.get_daily_time_index_array(time_index_range)
                    daily_qout = self.get_daily_qout_index(river_index, daily_time_index_array, mode=mode)
                    time_array = self.get_time_array()
                    for idx, time_idx in enumerate(daily_time_index_array):
                        current_day = time.gmtime(time_array[time_idx])
                        #write last average
                        writer.writerow([time.strftime("%Y/%m/%d", current_day), "{0:.5f}".format(daily_qout[idx])])
                else:
                    qout_arr = self.get_qout_index(river_index, time_index_array=time_index_range)
                    time_array = self.get_time_array(time_index_array=time_index_range)
                    with open(path_to_output_file, 'w') as outcsv:
                        for index in xrange(len(qout_arr)):
                            var_time = time.gmtime(time_array[index])
                            writer.writerow([time.strftime("%Y/%m/%d %H:00", var_time), "{0:.5f}".format(qout_arr[index])])

        else:
            print("Valid time variable not found. Printing values only ...")
            qout_arr = self.get_qout_index(river_index)
            with open_csv(path_to_output_file, 'w') as outcsv:
                writer = csv_writer(outcsv)
                for index in xrange(len(qout_arr)):
                    writer.writerow([index, "{0:.5f}".format(qout_arr[index])])

    def write_flows_to_gssha_time_series_xys(self, 
                                             path_to_output_file,
                                             series_name,
                                             series_id,
                                             river_index=None, 
                                             river_id=None,
                                             date_search_start=None,
                                             date_search_end=None,
                                             daily=False, 
                                             mode="mean"):
        """
        Write out RAPID output to GSSHA WMS time series xys file.
        
        Parameters:
            path_to_output_file(str): Path to the output xys file.
            series_name(str): The name for the series.
            series_id(int): The ID to give the series.
            river_index(Optional[datetime]): This is the index of the river in the file you want the streamflow for.
            river_id(Optional[datetime]): This is the river ID that you want the streamflow for.
            date_search_start(Optional[datetime]): This is a datetime object with the date of the minimum date for starting.
            date_search_end(Optional[datetime]): This is a datetime object with the date of the maximum date for ending.
            daily(Optional[boolean]): If True and the file is CF-Compliant, write out daily flows.
            mode(Optional[str]): You can get the daily average "mean" or the maximum "max". Defauls is "mean".

        Example writing entire time series to file:
        
        .. code:: python
        
            from RAPIDpy import RAPIDDataset
    
            river_id = 3624735
            path_to_rapid_qout = '/path/to/Qout.nc'

            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                #for 
                qout_nc.write_flows_to_gssha_time_series_xys('/timeseries/Qout_3624735.xys',
                                                             series_name="RAPID_TO_GSSHA_{0}".format(river_id),
                                                             series_id=34,
                                                             river_id=river_id,
                                                             )
                                      
        Example writing entire time series as daily average to file:
        
        .. code:: python
        
            from RAPIDpy import RAPIDDataset
    
            river_id = 3624735
            path_to_rapid_qout = '/path/to/Qout.nc'

            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                #NOTE: Getting the river index is not necessary
                #this is just an example of how to use this                                      
                river_index = qout_nc.get_river_index(river_id)
                
                #if file is CF compliant, you can write out daily average
                qout_nc.write_flows_to_gssha_time_series_xys('/timeseries/Qout_daily.xys',
                                                             series_name="RAPID_TO_GSSHA_{0}".format(river_id),
                                                             series_id=34,
                                                             river_index=river_index,
                                                             daily=True,
                                                             )

        Example writing subset of time series as daily maximum to file:
        
        .. code:: python
        
            from datetime import datetime
            from RAPIDpy import RAPIDDataset
            
            river_id = 3624735
            path_to_rapid_qout = '/path/to/Qout.nc'

            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                #NOTE: Getting the river index is not necessary
                #this is just an example of how to use this                                      
                river_index = qout_nc.get_river_index(river_id)
                
                #if file is CF compliant, you can filter by date and get daily values
                qout_nc.write_flows_to_gssha_time_series_xys('/timeseries/Qout_daily_date_filter.xys',
                                                             series_name="RAPID_TO_GSSHA_{0}".format(river_id),
                                                             series_id=34,
                                                             river_index=river_index,
                                                             date_search_start=datetime(2002, 8, 31),
                                                             date_search_end=datetime(2002, 9, 15),
                                                             daily=True,
                                                             mode="max"
                                                             )        
        """
        if river_id != None:
            river_index = self.get_river_index(river_id)
        elif river_id == None and river_index == None:
            raise Exception("ERROR: Need reach id or reach index ...")

        #analyze and write
        if self.is_time_variable_valid():
            time_index_range = self.get_time_index_range(date_search_start=date_search_start,
                                                         date_search_end=date_search_end)
            with open_csv(path_to_output_file, 'w') as out_ts:
                if daily:
                    daily_time_index_array = self.get_daily_time_index_array(time_index_range)
                    daily_qout = self.get_daily_qout_index(river_index, daily_time_index_array, mode=mode)
                    out_ts.write("XYS {0} {1} \"{2}\"\r\n".format(series_id, len(daily_qout), series_name))
                    time_array = self.get_time_array()
                    for idx, time_idx in enumerate(daily_time_index_array):
                        date_str = time.strftime("%m/%d/%Y %I:%M:%S %p", time.gmtime(time_array[time_idx]))
                        out_ts.write("\"{0}\" {1:.5f}\n".format(date_str, daily_qout[idx]))
                else:
                    qout_arr = self.get_qout_index(river_index, time_index_array=time_index_range)
                    out_ts.write("XYS {0} {1} \"{2}\"\r\n".format(series_id, len(qout_arr), series_name))
                    time_array = self.get_time_array(time_index_array=time_index_range)
                    for index in xrange(len(qout_arr)):
                        date_str = time.strftime("%m/%d/%Y %I:%M:%S %p", time.gmtime(time_array[index]))
                        out_ts.write("\"{0}\" {1:.5f}\n".format(date_str, qout_arr[index]))
        else:
            raise IndexError("Valid time variable not found. Valid time variable required in Qout file to proceed ...")

    def write_flows_to_gssha_time_series_ihg(self, 
                                             path_to_output_file,
                                             connection_list,
                                             date_search_start=None,
                                             date_search_end=None,
                                             daily=False, 
                                             mode="mean"):
        """
        Write out RAPID output to GSSHA time series ihg file
        
        .. note:: See: http://www.gsshawiki.com/Surface_Water_Routing:Introducing_Dischage/Constituent_Hydrographs
        
        .. note:: GSSHA project card is CHAN_POINT_INPUT
        
        Parameters:
            path_to_output_file(str): Path to the output xys file.
            connection_list(list): List of dictionaries with link_id, node_id, baseflow, and rapid_rivid.
            date_search_start(Optional[datetime]): This is a datetime object with the date of the minimum date for starting.
            date_search_end(Optional[datetime]): This is a datetime object with the date of the maximum date for ending.
            daily(Optional[boolean]): If True and the file is CF-Compliant, write out daily flows.
            mode(Optional[str]): You can get the daily average "mean" or the maximum "max". Defauls is "mean".

        Example writing entire time series to file:
        
        .. code:: python
        
            from RAPIDpy import RAPIDDataset
    
            path_to_rapid_qout = '/path/to/Qout.nc'
            
            #list to connect the RAPID rivers to GSSHA rivers
            connection_list = [
                               {
                                 'link_id': 599,
                                 'node_id': 1,
                                 'baseflow': 0.0,
                                 'rapid_rivid': 80968,
                               },
                               {
                                 'link_id': 603,
                                 'node_id': 1,
                                 'baseflow': 0.0,
                                 'rapid_rivid': 80967,
                               },
                             ]
                         
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                #for writing entire time series to file
                qout_nc.write_flows_to_gssha_time_series_ihg('/timeseries/Qout_3624735.ihg',
                                                             connection_list,
                                                             )
                                      
        Example writing entire time series as daily average to file:
        
        .. code:: python
        
            from RAPIDpy import RAPIDDataset
    
            path_to_rapid_qout = '/path/to/Qout.nc'
            
            #list to connect the RAPID rivers to GSSHA rivers
            connection_list = [
                               {
                                 'link_id': 599,
                                 'node_id': 1,
                                 'baseflow': 0.0,
                                 'rapid_rivid': 80968,
                               },
                              ]
                         
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                #if file is CF compliant, you can write out daily average
                qout_nc.write_flows_to_gssha_time_series_ihg('/timeseries/Qout_3624735.ihg',
                                                             connection_list,
                                                             daily=True,
                                                             )
                                                             
        
        Example writing subset of time series as daily maximum to file:
        
        .. code:: python
        
            from datetime import datetime
            from RAPIDpy import RAPIDDataset
    
            path_to_rapid_qout = '/path/to/Qout.nc'
            
            #list to connect the RAPID rivers to GSSHA rivers
            connection_list = [
                               {
                                 'link_id': 599,
                                 'node_id': 1,
                                 'baseflow': 0.0,
                                 'rapid_rivid': 80968,
                               },
                             ]
                         
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                #if file is CF compliant, you can filter by date and get daily values
                qout_nc.write_flows_to_gssha_time_series_ihg('/timeseries/Qout_daily_date_filter.ihg',
                                                             connection_list,
                                                             date_search_start=datetime(2002, 8, 31),
                                                             date_search_end=datetime(2002, 9, 15),
                                                             daily=True,
                                                             mode="max"
                                                             )        
        """
        #analyze and write
        if self.is_time_variable_valid():
            time_index_range = self.get_time_index_range(date_search_start=date_search_start,
                                                         date_search_end=date_search_end)
            with open_csv(path_to_output_file, 'w') as out_ts:
                #####HEADER SECTION EXAMPLE:
                #NUMPT 3
                #POINT 1 599 0.0
                #POINT 1 603 0.0
                #POINT 1 605 0.0
                
                out_ts.write("NUMPT {0}\n".format(len(connection_list)))
                river_idx_list = []
                for connection in connection_list:
                    out_ts.write("POINT {0} {1} {2}\n".format(connection['node_id'],
                                                              connection['link_id'], 
                                                              connection['baseflow'],
                                                              ))
                    river_idx_list.append(self.get_river_index(int(connection['rapid_rivid'])))
                
                
                #####INFLOW SECTION EXAMPLE:
                #NRPDS 54
                #INPUT 2002 01 01 00 00 15.551210 12.765090 0.000000
                #INPUT 2002 01 02 00 00 15.480830 12.765090 0.000000
                #INPUT 2002 01 03 00 00 16.078910 12.765090 0.000000
                # ...
                if daily:
                    daily_time_index_array = self.get_daily_time_index_array(time_index_range)
                    out_ts.write("NRPDS {0}\n".format(len(daily_time_index_array)))
                    daily_qout_2d_array = self.get_daily_qout_index(river_idx_list, daily_time_index_array, mode=mode)
                    qout_num_dimensions = len(daily_qout_2d_array.shape)
                    time_array = self.get_time_array()
                    for idx, time_idx in enumerate(daily_time_index_array):
                        date_str = time.strftime("%Y %m %d %H %M", time.gmtime(time_array[time_idx]))
                        if qout_num_dimensions > 1:
                            qout_str = " ".join(["{0:.5f}".format(daily_qout) for daily_qout in daily_qout_2d_array[:, idx]])
                        else:
                            qout_str = "{0:.5f}".format(daily_qout_2d_array[idx])
                            
                        out_ts.write("INPUT {0} {1}\n".format(date_str, qout_str))
                else:
                    qout_2d_array = self.get_qout_index(river_idx_list, time_index_array=time_index_range)
                    qout_num_dimensions = len(qout_2d_array.shape)
                    num_time_steps = qout_2d_array.shape[0]
                    if qout_num_dimensions > 1:
                        num_time_steps = qout_2d_array.shape[1]
                    
                    out_ts.write("NRPDS {0}\n".format(num_time_steps))
                    
                    time_array = self.get_time_array(time_index_array=time_index_range)
                    for index in xrange(num_time_steps):
                        date_str = time.strftime("%Y %m %d %H %M", time.gmtime(time_array[index]))
                        if qout_num_dimensions > 1:
                            qout_str = " ".join(["{0:.5f}".format(daily_qout) for daily_qout in qout_2d_array[:, index]])
                        else:
                            qout_str = "{0:.5f}".format(qout_2d_array[index])
                        
                        out_ts.write("INPUT {0} {1}\n".format(date_str, qout_str))
        else:
            raise IndexError("Valid time variable not found. Valid time variable required in Qout file to proceed ...")
