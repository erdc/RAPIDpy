# -*- coding: utf-8 -*-
##
##  dataset.py
##  RAPIDpy
##
##  Created by Alan D Snow, 2016.
##  Copyright © 2016 Alan D Snow. All rights reserved.
##

import datetime
from netCDF4 import Dataset
import numpy as np
from numpy.ma import masked
from pytz import utc

#------------------------------------------------------------------------------
#Main Dataset Manager Class
#------------------------------------------------------------------------------
class RAPIDDataset(object):
    """
    This class is designed to access data from the RAPID Qout
    netCDF file
    """
    def __init__(self, filename):
        """
        Initialize the class with variables given by the user
        """
        self.qout_nc = Dataset(filename, mode='r')
        
        #determine river ID dimension
        if 'rivid' in self.qout_nc.dimensions:
            self.river_id_dimension = 'rivid'
        elif 'COMID' in self.qout_nc.dimensions:
            self.river_id_dimension = 'COMID'
        elif 'DrainLnID' in self.qout_nc.dimensions:
            self.river_id_dimension = 'DrainLnID'
        elif 'FEATUREID' in self.qout_nc.dimensions:
            self.river_id_dimension = 'FEATUREID'
        else:
            raise Exception('ERROR: Could not find river ID dimension.')
        self.size_river_id = len(self.qout_nc.dimensions[self.river_id_dimension])
        
        
        #determine time dimension
        if 'time' in self.qout_nc.dimensions:
            self.size_time = len(self.qout_nc.dimensions['time'])
        elif 'Time' in self.qout_nc.dimensions:
            self.size_time = len(self.qout_nc.dimensions['Time'])
        else:
            raise Exception('ERROR: Could not find time dimension.')

        #determin streamflow variable
        if 'Qout' in self.qout_nc.variables:
            self.q_var_name = 'Qout'
        elif 'm3_riv' in self.qout_nc.variables:
            self.q_var_name = 'm3_riv'
        else:
            raise Exception('ERROR: Could not find flow variable. Looked for Qout and m3_riv.')

        self.size_q_var = len(self.qout_nc.variables[self.q_var_name])

    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        self.qout_nc.close()

    def is_time_variable_valid(self):
        """
        This function returns whether or not the time variable
        is valid
        """
        time_var_valid = False
        if 'time' in self.qout_nc.variables.keys():
            if len(self.qout_nc.dimensions['time'])>0:
                if not (self.qout_nc.variables['time'][:] == masked).any():
                    try:
                        timestep = (datetime.datetime.utcfromtimestamp(self.qout_nc.variables['time'][0])-
                                    datetime.datetime.utcfromtimestamp(self.qout_nc.variables['time'][1])).seconds
                        
                        if timestep > 0:
                            time_var_valid = True
                    except ValueError:
                        pass
    
        return time_var_valid

    def get_time_array(self, datetime_simulation_start=None,
                       simulation_time_step_seconds=None,
                       return_datetime=False):
        """
        This method extracts or generates an array of time
        """
        time_array = []

        #CF-1.6 compliant file
        if self.is_time_variable_valid():
            time_array = self.qout_nc.variables['time'][:]
        #Original Qout file
        elif datetime_simulation_start is not None and simulation_time_step_seconds is not None:
            initial_time_seconds = (datetime_simulation_start.replace(tzinfo=utc)-
                                    datetime.datetime(1970,1,1, tzinfo=utc)).total_seconds()+simulation_time_step_seconds
            final_time_seconds = initial_time_seconds + self.size_time*simulation_time_step_seconds
            time_array = np.arange(initial_time_seconds, final_time_seconds, simulation_time_step_seconds)
        else:
            raise Exception("ERROR: This file does not contain the time variable."
                            " To get time array, add datetime_simulation_start"
                            " and simulation_time_step_seconds")
        
        if not return_datetime:
            return time_array
        else:
            return [datetime.datetime.utcfromtimestamp(t) for t in time_array]

    def get_time_index_range(self, date_search_start=None,
                             date_search_end=None,
                             time_index_start=None,
                             time_index_end=None,
                             time_index=None):
        """
        Generates a time index range based on datetimes
        """
        #get the range of time based on datetime range
        time_range = None
        if 'time' in self.qout_nc.variables and (date_search_start is not None or date_search_end is not None):
            print "Determining time range ({0} to {1})...".format(date_search_start, date_search_end)
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
        
        return time_range

    def get_river_id_array(self):
        """
        This method returns the river index array for this file
        """
        return self.qout_nc.variables[self.river_id_dimension][:]
    
    def get_river_index(self, river_id):
        """
        This method retrieves the river index in the netCDF
        dataset corresponding to the river ID
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
                print "ReachID", river_id, "not found in netCDF dataset. Skipping ..."
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
        This method extracts streamflow data by river id
        It allows for extracting single or multiple river streamflow arrays
        It has options to extract by date or by date index
        """
        #get indices of where the streamflow data is
        riverid_index_list_subset = None
        if river_id_array is not None:
            if not isinstance(river_id_array, list) \
            or not type(river_id_array).__module__ == np.array:
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
    
    def get_daily_qout(self, reach_index, steps_per_group=1):
        """
        Gets the daily time series from RAPID output
        """
        if self.is_time_variable_valid() and steps_per_group<=1:
            current_day = datetime.datetime.utcfromtimestamp(self.qout_nc.variables['time'][0])
            flow = 0.0
            num_days = 0
            qout_arr = self.get_qout_index(reach_index)
            daily_qout = []
            for idx, t in enumerate(self.get_time_array()):
                var_time = datetime.datetime.utcfromtimestamp(t)
                if current_day.day == var_time.day:
                    flow += qout_arr[idx]
                    num_days += 1
                else:
                    if num_days > 0:
                        #write last average
                        daily_qout.append(flow/num_days)
                    
                    #start new average
                    current_day = var_time
                    num_days = 1
                    flow = qout_arr[idx]
            
            return np.array(daily_qout, np.float32)
        elif steps_per_group > 1:
            flow_data = self.get_qout_index(reach_index)
            daily_qout = []
            for step_index in xrange(0, len(flow_data), steps_per_group):
                flows_slice = flow_data[step_index:step_index + steps_per_group]
                daily_qout.append(np.mean(flows_slice))
            return np.array(daily_qout, np.float32)
        else:
            print steps_per_group
            raise Exception("Must have steps_per_group set to a value greater than one "
                            "due to non CF-Compliant Qout file ...")

