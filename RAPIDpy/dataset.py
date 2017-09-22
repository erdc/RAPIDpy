# -*- coding: utf-8 -*-
"""
   dataset.py
   RAPIDpy

   Created by Alan D Snow, 2016.
   License: BSD-3-Clause
"""
from csv import writer as csv_writer
import datetime

from netCDF4 import Dataset, num2date
import numpy as np
from numpy.ma import is_masked
import pandas as pd
from past.builtins import xrange  # pylint: disable=redefined-builtin
from pytz import utc

from .helper_functions import log, open_csv


# -----------------------------------------------------------------------------
# Helper Function
# -----------------------------------------------------------------------------
def compare_qout_files(dataset1_path, dataset2_path):
    """
    This function compares the output of RAPID Qout and tells you where
    they are different.
    """
    qout_same = False

    d1 = RAPIDDataset(dataset1_path)
    d2 = RAPIDDataset(dataset2_path)

    if len(d1.get_river_id_array()) != len(d2.get_river_id_array()):
        log("Length of COMID/rivid input not the same.",
            "ERROR")

    if not (d1.get_river_id_array() == d2.get_river_id_array()).all():
        log("COMID/rivid order is different in each dataset."
            " Reordering data for comparison.",
            "WARNING")

        d2_reordered_river_index_list = []
        for rivid in d1.get_river_id_array():
            reordered_index = np.where(d2.get_river_id_array() == rivid)[0][0]
            d2_reordered_river_index_list.append(reordered_index)
        d2_reordered_qout = d2.get_qout_index(d2_reordered_river_index_list)
    else:
        d2_reordered_qout = d2.get_qout()

    # get where the files are different
    d1_qout = d1.get_qout()
    where_diff = np.where(d1_qout != d2_reordered_qout)
    un_where_diff = np.unique(where_diff[0])

    # if different, check to see how different
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
                decimal_test = -1
            except AssertionError as ex:
                if decimal_test <= 1:
                    log(ex, "WARNING")
                decimal_test -= 1

        log("Number of different timeseries: {0}".format(len(un_where_diff)),
            "INFO")
        log("COMID idexes where different: {0}".format(un_where_diff),
            "INFO")
        log("COMID idexes where different: {0}".format(un_where_diff),
            "INFO")
        index = un_where_diff[0]
        log("Dataset 1 example. COMID index: "
            "{0}".format(d1.get_qout_index(index)),
            "INFO")
        log("Dataset 2 example. COMID index: "
            "{0}".format(d2_reordered_qout[index, :]),
            "INFO")

    else:
        qout_same = True
        log("Output Qout data is the same.",
            "INFO")

    d1.close()
    d2.close()
    return qout_same


# ------------------------------------------------------------------------------
# Main Dataset Manager Class
# ------------------------------------------------------------------------------
class RAPIDDataset(object):
    """
    This class is designed to access data from the RAPID Qout
    NetCDF file.

    Attributes
    ----------
    filename: str
        Path to the RAPID Qout NetCDF file.
    river_id_dimension: str, optional
        Name of the river ID dimension. Default is to search through
        a pre-defined list.
    river_id_variable: str, optional
        Name of the river ID variable. Default is to search through
        a pre-defined list.
    streamflow_variable: str, optional
        Name of the streamflow varaible. Default is to search through
        a pre-defined list.
    datetime_simulation_start: :obj:`datetime.datetime`, optional
        This is a datetime object with the date of the simulation start time.
    simulation_time_step_seconds: int, optional
        This is the time step of the simulation output in seconds.
    out_tzinfo: tzinfo, optional
        Time zone to output data as. The dates will be converted from UTC
        to the time zone input. Default is UTC.


    Example::

        from RAPIDpy import RAPIDDataset

        path_to_rapid_qout = '/path/to/Qout.nc'
        with RAPIDDataset(path_to_rapid_qout) as qout_nc:
            #USE FUNCTIONS TO ACCESS DATA HERE

    """
    # pylint: disable=too-many-instance-attributes
    def __init__(self, filename,
                 river_id_dimension="",
                 river_id_variable="",
                 streamflow_variable="",
                 datetime_simulation_start=None,
                 simulation_time_step_seconds=None,
                 out_tzinfo=None):
        """
        Initialize the class with variables given by the user
        """
        self.qout_nc = Dataset(filename, mode='r')

        # determine river ID dimension
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
                raise IndexError('Could not find river ID dimension.')
        elif river_id_dimension not in self.qout_nc.dimensions:
            raise IndexError('Could not find river ID dimension:'
                             ' {0}.'.format(river_id_dimension))

        self.size_river_id = len(self.qout_nc
                                     .dimensions[self.river_id_dimension])

        variable_keys = self.qout_nc.variables.keys()

        # determine streamflow variable
        self.q_var_name = streamflow_variable
        if not streamflow_variable:
            if 'Qout' in variable_keys:
                self.q_var_name = 'Qout'
            elif 'streamflow' in variable_keys:
                self.q_var_name = 'streamflow'
            elif 'm3_riv' in variable_keys:
                self.q_var_name = 'm3_riv'
            else:
                raise IndexError('ERROR: Could not find flow variable.'
                                 ' Looked for Qout, streamflow, and m3_riv.')
        elif streamflow_variable not in variable_keys:
            raise IndexError('Could not find flow variable.'
                             ' Looked for {0}.'.format(streamflow_variable))

        self.size_q_var = len(self.qout_nc.variables[self.q_var_name])

        # determine time dimension
        if 'time' in self.qout_nc.dimensions:
            self.size_time = len(self.qout_nc.dimensions['time'])
        elif 'Time' in self.qout_nc.dimensions:
            self.size_time = len(self.qout_nc.dimensions['Time'])
        else:
            raise IndexError('Could not find time dimension.')

        # determine river ID variable
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
                log('Could not find river ID variable'
                    ' in {0}.'.format(variable_keys),
                    "WARNING")
        elif river_id_variable not in variable_keys:
            log('Could not find river ID variable:'
                ' {0}.'.format(river_id_variable),
                "WARNING")

        self.out_tzinfo = out_tzinfo
        self.datetime_simulation_start = datetime_simulation_start
        self.simulation_time_step_seconds = simulation_time_step_seconds

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        """Close the dataset."""
        self.qout_nc.close()

    def _is_legacy_time_valid(self):
        """
        This determines whether or not legacy time is set correctly.

        Returns
        -------
        boolean:
            True if the legacy time is setup correctly, otherwise false.
        """
        return self.datetime_simulation_start is not None and \
            self.simulation_time_step_seconds is not None

    def is_time_variable_valid(self):
        """
        This function returns whether or not the time variable
        is valid.

        Returns
        -------
        boolean
            True if the time variable is valid, otherwise false.


        Example::

            from RAPIDpy import RAPIDDataset

            path_to_rapid_qout = '/path/to/Qout.nc'
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                if qout_nc.is_time_variable_valid():
                    #DO WORK HERE

        """
        # pylint: disable=len-as-condition
        time_var_valid = False
        if 'time' in self.qout_nc.variables.keys():
            if len(self.qout_nc.dimensions['time']) > 0:
                if not is_masked(self.qout_nc.variables['time'][:]):
                    try:
                        timestep = (datetime.datetime
                                    .utcfromtimestamp(
                                        self.qout_nc.variables['time'][1]
                                    ) -
                                    datetime.datetime
                                    .utcfromtimestamp(
                                        self.qout_nc.variables['time'][0]
                                    )).total_seconds()
                        if timestep > 0:
                            time_var_valid = True
                    except ValueError:
                        pass

        return time_var_valid

    def raise_time_valid(self):
        """Raise ValueError if time not valid"""
        if not (self.is_time_variable_valid() or self._is_legacy_time_valid()):
            raise IndexError("Valid time variable not found. Valid time"
                             " variable required in Qout file to proceed ...")

    def get_time_array(self,
                       datetime_simulation_start=None,
                       simulation_time_step_seconds=None,
                       return_datetime=False,
                       time_index_array=None):
        """
        This method extracts or generates an array of time.
        The new version of RAPID output has the time array stored.
        However, the old version requires the user to know when the
        simulation began and the time step of the output.

        Parameters
        ----------
        datetime_simulation_start: :obj:`datetime.datetime`, optional
            The start datetime of the simulation. Only required if the time
            variable is not included in the file.
        simulation_time_step_seconds: int, optional
            The time step of the file in seconds. Only required if the time
            variable is not included in the file.
        return_datetime: bool, optional
            If true, it converts the data to a list of datetime objects.
            Default is False.
        time_index_array: list or :obj:`numpy.array`, optional
            This is used to extract the datetime values by index from the main
            list. This can be from the *get_time_index_range* function.

        Returns
        -------
        list:
            An array of integers representing seconds since Jan 1, 1970 UTC
            or datetime objects if *return_datetime* is set to True.

        These examples demonstrates how to retrieve or generate a time array
        to go along with your RAPID streamflow series.


        CF-Compliant Qout File Example:

        .. code:: python

            from RAPIDpy import RAPIDDataset

            path_to_rapid_qout = '/path/to/Qout.nc'
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                #retrieve integer timestamp array
                time_array = qout_nc.get_time_array()

                #or, to get datetime array
                time_datetime = qout_nc.get_time_array(return_datetime=True)


        Legacy Qout File Example:

        .. code:: python

            from RAPIDpy import RAPIDDataset

            path_to_rapid_qout = '/path/to/Qout.nc'
            with RAPIDDataset(path_to_rapid_qout,
                              datetime_simulation_start=datetime(1980, 1, 1),
                              simulation_time_step_seconds=3 * 3600)\
                    as qout_nc:

                #retrieve integer timestamp array
                time_array = qout_nc.get_time_array()

                #or, to get datetime array
                time_datetime = qout_nc.get_time_array(return_datetime=True)

        """
        # Original Qout file
        if datetime_simulation_start is not None:
            self.datetime_simulation_start = datetime_simulation_start
        if simulation_time_step_seconds is not None:
            self.simulation_time_step_seconds = simulation_time_step_seconds

        epoch = datetime.datetime(1970, 1, 1, tzinfo=utc)
        time_units = "seconds since {0}".format(epoch)

        # CF-1.6 compliant file
        if self.is_time_variable_valid():
            time_array = self.qout_nc.variables['time'][:]
            if self.qout_nc.variables['time'].units:
                time_units = self.qout_nc.variables['time'].units

        # Original Qout file
        elif self._is_legacy_time_valid():
            initial_time_seconds = ((self.datetime_simulation_start
                                    .replace(tzinfo=utc) - epoch)
                                    .total_seconds() +
                                    self.simulation_time_step_seconds)
            final_time_seconds = (initial_time_seconds +
                                  self.size_time *
                                  self.simulation_time_step_seconds)
            time_array = np.arange(initial_time_seconds,
                                   final_time_seconds,
                                   self.simulation_time_step_seconds)
        else:
            raise ValueError("This file does not contain the time"
                             " variable. To get time array, add"
                             " datetime_simulation_start and"
                             " simulation_time_step_seconds")

        if time_index_array is not None:
            time_array = time_array[time_index_array]

        if return_datetime:
            time_array = num2date(time_array, time_units)

            if self.out_tzinfo is not None:
                for i in xrange(len(time_array)):
                    # convert time to output timezone
                    time_array[i] = utc.localize(time_array[i]) \
                                       .astimezone(self.out_tzinfo) \
                                       .replace(tzinfo=None)

        return time_array

    def get_time_index_range(self,
                             date_search_start=None,
                             date_search_end=None,
                             time_index_start=None,
                             time_index_end=None,
                             time_index=None):
        """
        Generates a time index range based on time bounds given.
        This is useful for subset data extraction.

        Parameters
        ----------
        date_search_start: :obj:`datetime.datetime`, optional
            This is a datetime object with the date of the minimum date for
            starting.
        date_search_end: :obj:`datetime.datetime`, optional
            This is a datetime object with the date of the maximum date
            for ending.
        time_index_start: int, optional
            This is the index of the start of the time array subset.
            Useful for the old file version.
        time_index_end: int, optional
            This is the index of the end of the time array subset.
            Useful for the old file version.
        time_index: int, optional
            This is the index of time to return in the case that your
            code only wants one index. Used internally.

        Returns
        -------
        :obj:`numpy.array`:
            This is an array of time indices used to extract a subset of data.


        CF-Compliant Qout File Example:

        .. code:: python

            from datetime import datetime
            from RAPIDpy import RAPIDDataset

            path_to_rapid_qout = '/path/to/Qout.nc'
            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                time_index_range = qout_nc.get_time_index_range(
                    date_search_start=datetime(1980, 1, 1),
                    date_search_end=datetime(1980, 12, 11))


        Legacy Qout File Example:

        .. code:: python

            from datetime import datetime
            from RAPIDpy import RAPIDDataset

            path_to_rapid_qout = '/path/to/Qout.nc'
            with RAPIDDataset(path_to_rapid_qout,
                              datetime_simulation_start=datetime(1980, 1, 1),
                              simulation_time_step_seconds=3600) as qout_nc:

                time_index_range = qout_nc.get_time_index_range(
                    date_search_start=datetime(1980, 1, 1),
                    date_search_end=datetime(1980, 12, 11))

        """
        # get the range of time based on datetime range
        time_range = None
        if ((self.is_time_variable_valid() or self._is_legacy_time_valid()) and
                (date_search_start is not None or
                 date_search_end is not None)):

            log("Determining time range ({0} to {1})"
                "...".format(date_search_start, date_search_end),
                "INFO")
            time_array = self.get_time_array()
            if date_search_start is not None:
                date_search_start_utc = date_search_start
                if self.out_tzinfo is not None:
                    date_search_start_utc = self.out_tzinfo \
                                                .localize(date_search_start) \
                                                .astimezone(utc) \
                                                .replace(tzinfo=None)
                seconds_start = (date_search_start_utc -
                                 datetime.datetime(1970, 1, 1)).total_seconds()
                time_range = np.where(time_array >= seconds_start)[0]

            if date_search_end is not None:
                date_search_end_utc = date_search_end
                if self.out_tzinfo is not None:
                    date_search_end_utc = self.out_tzinfo \
                                              .localize(date_search_end) \
                                              .astimezone(utc) \
                                              .replace(tzinfo=None)

                seconds_end = (date_search_end_utc -
                               datetime.datetime(1970, 1, 1)).total_seconds()
                if time_range is not None:
                    time_range = np.intersect1d(time_range,
                                                np.where(time_array <=
                                                         seconds_end)[0])
                else:
                    time_range = np.where(time_array <= seconds_end)[0]

        # get the range of time based on time index range
        elif time_index_start is not None or time_index_end is not None:
            if time_index_start is None:
                time_index_start = 0
            if time_index_end is None:
                time_index_end = self.size_time
            time_range = range(time_index_start, time_index_end)

        # get only one time step
        elif time_index is not None:
            time_range = [time_index]
        # return all
        else:
            time_range = range(self.size_time)

        return time_range

    def get_river_id_array(self):
        """
        This method returns the river ID array for this file.

        Returns
        -------
        :obj:`numpy.array`:
            An array of the river ID's


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

        Parameters
        ----------
        river_id: int
            The ID of the river segment.

        Returns
        -------
        int:
            The index of the river ID's in the file.


        Example::

            from RAPIDpy import RAPIDDataset

            path_to_rapid_qout = '/path/to/Qout.nc'
            river_id = 53458

            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                river_index = qout_nc.get_river_index(river_id)

        """
        try:
            return np.where(self.get_river_id_array() == river_id)[0][0]
        except IndexError:
            raise IndexError("ERROR: River ID {0} not found in dataset "
                             "...".format(river_id))

    def get_subset_riverid_index_list(self, river_id_list):
        """
        Gets the subset riverid_list from the netcdf file
        Optional returns include the list of valid river ids in the dataset
        as well as a list of missing rive rids

        Parameters
        ----------
        river_id_list: list or :obj:`numpy.array`
            Array of river ID's for the river segments you want the index of.

        Returns
        -------
        :obj:`numpy.array`
            A sorted array of the river index in the NetCDF file that
            were found.
        :obj:`numpy.array`
            A sorted array of the river IDs that were found.
        list
            An array of the missing river ids.

        """
        netcdf_river_indices_list = []
        valid_river_ids = []
        missing_river_ids = []
        for river_id in river_id_list:
            # get where streamids are in netcdf file
            try:
                netcdf_river_indices_list \
                    .append(self.get_river_index(river_id))
                valid_river_ids.append(river_id)
            except IndexError:
                log("ReachID {0} not found in netCDF dataset."
                    " Skipping ...".format(river_id),
                    "WARNING")
                missing_river_ids.append(river_id)

        np_valid_river_indices_list = np.array(netcdf_river_indices_list)
        np_valid_river_ids = np.array(valid_river_ids)
        sorted_indexes = np.argsort(np_valid_river_indices_list)

        return(np_valid_river_indices_list[sorted_indexes],
               np_valid_river_ids[sorted_indexes],
               np.array(missing_river_ids))

    def get_qout(self,
                 river_id_array=None,
                 date_search_start=None,
                 date_search_end=None,
                 time_index_start=None,
                 time_index_end=None,
                 time_index=None,
                 time_index_array=None,
                 daily=False,
                 pd_filter=None,
                 filter_mode="mean",
                 as_dataframe=False):
        """
        This method extracts streamflow data by a single river ID
        or by a river ID array. It has options to extract by date
        or by date index.

        Parameters
        ----------
        river_id_array: :obj:`numpy.array` or list or int, optional
            A single river ID or an array of river IDs.
        date_search_start: :obj:`datetime.datetime`, optional
            This is a datetime object with the date of the minimum date
            for starting.
        date_search_end: :obj:`datetime.datetime`, optional
            This is a datetime object with the date of the maximum date
            for ending.
        time_index_start: int, optional
            This is the index of the start of the time array subset.
            Useful for the old file version.
        time_index_end: int, optional
            This is the index of the end of the time array subset.
            Useful for the old file version.
        time_index: int, optional
            This is the index of time to return in the case that your
            code only wants one index. Used internally.
        time_index_array: list or :obj:`numpy.array`, optional
            This is used to extract the vales only for particular dates.
            This can be from the *get_time_index_range* function.
        daily: bool, optional
            If true, this will convert qout to daily average.
        pd_filter: str, optional
            This is a valid pandas resample frequency filter.
        filter_mode: str, optional
            You can get the daily average "mean" or the maximum "max".
            Default is "mean".
        as_dataframe: bool, optional
            Return as a pandas dataframe object. Default is False.


        Returns
        -------
        qout_array: :obj:`numpy.array`
            This is a 1D or 2D array or a single value depending on your
            input search.


        This example demonstrates how to retrieve the streamflow associated
        with the reach you are interested in::

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
                streamflow_array = qout_nc.get_qout(
                    river_id,
                    date_search_start=datetime(1985,1,1),
                    date_search_end=datetime(1985,2,4))

        """
        # get indices of where the streamflow data is
        riverid_index_list_subset = None
        if river_id_array is not None:
            if not hasattr(river_id_array, "__len__"):
                river_id_array = [river_id_array]
            riverid_index_list_subset = \
                self.get_subset_riverid_index_list(river_id_array)[0]

        return self.get_qout_index(riverid_index_list_subset,
                                   date_search_start,
                                   date_search_end,
                                   time_index_start,
                                   time_index_end,
                                   time_index,
                                   time_index_array,
                                   daily,
                                   pd_filter,
                                   filter_mode,
                                   as_dataframe)

    def get_qout_index(self,
                       river_index_array=None,
                       date_search_start=None,
                       date_search_end=None,
                       time_index_start=None,
                       time_index_end=None,
                       time_index=None,
                       time_index_array=None,
                       daily=False,
                       pd_filter=None,
                       filter_mode="mean",
                       as_dataframe=False):
        """
        This method extracts streamflow data by river index.
        It allows for extracting single or multiple river streamflow arrays
        It has options to extract by date or by date index.

        See: :meth:`RAPIDpy.RAPIDDataset.get_qout`
        """
        if river_index_array is not None:
            if hasattr(river_index_array, "__len__"):
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
        if qout_dimensions[0].lower() == 'time' and \
                qout_dimensions[1].lower() == self.river_id_dimension.lower():
            if time_index_array is not None and river_index_array is not None:
                streamflow_array = qout_variable[time_index_array,
                                                 river_index_array].transpose()
            elif time_index_array is not None:
                streamflow_array = qout_variable[time_index_array, :] \
                                   .transpose()
            elif river_index_array is not None:
                streamflow_array = qout_variable[:, river_index_array] \
                                   .transpose()
            else:
                streamflow_array = qout_variable[:].transpose()
        elif qout_dimensions[1].lower() == 'time' and \
                qout_dimensions[0].lower() == self.river_id_dimension.lower():
            if time_index_array is not None and river_index_array is not None:
                streamflow_array = qout_variable[river_index_array,
                                                 time_index_array]
            elif time_index_array is not None:
                streamflow_array = qout_variable[:, time_index_array]
            elif river_index_array is not None:
                streamflow_array = qout_variable[river_index_array, :]
            else:
                streamflow_array = qout_variable[:]
        else:
            raise Exception("Invalid RAPID Qout file dimensions ...")

        if daily:
            pd_filter = "D"

        if pd_filter is not None or as_dataframe:
            time_array = self.get_time_array(return_datetime=True,
                                             time_index_array=time_index_array)
            qout_df = pd.DataFrame(streamflow_array.T, index=time_array)

            if pd_filter is not None:
                qout_df = qout_df.resample(pd_filter)
                if filter_mode == "mean":
                    qout_df = qout_df.mean()
                elif filter_mode == "max":
                    qout_df = qout_df.max()
                else:
                    raise Exception("Invalid filter_mode ...")

            if as_dataframe:
                return qout_df

            streamflow_array = qout_df.as_matrix().T

            if streamflow_array.ndim > 0 and streamflow_array.shape[0] == 1:
                streamflow_array = streamflow_array[0]

        return streamflow_array

    def write_flows_to_csv(self, path_to_output_file,
                           river_index=None,
                           river_id=None,
                           date_search_start=None,
                           date_search_end=None,
                           daily=False,
                           filter_mode="mean"):
        """
        Write out RAPID output to CSV file.

        .. note:: Need either *reach_id* or *reach_index* parameter,
                  but either can be used.

        Parameters
        ----------
        path_to_output_file: str
            Path to the output csv file.
        river_index: :obj:`datetime.datetime`, optional
            This is the index of the river in the file you want the
            streamflow for.
        river_id: :obj:`datetime.datetime`, optional
            This is the river ID that you want the streamflow for.
        date_search_start: :obj:`datetime.datetime`, optional
            This is a datetime object with the date of the minimum date
            for starting.
        date_search_end: :obj:`datetime.datetime`, optional
            This is a datetime object with the date of the maximum date
            for ending.
        daily: bool, optional
            If True and the file is CF-Compliant, write out daily flows.
        filter_mode: str, optional
            You can get the daily average "mean" or the maximum "max".
            Default is "mean".


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


                #if file is CF compliant, you can write out daily average

                #NOTE: Getting the river index is not necessary
                #this is just an example of how to use this
                river_index = qout_nc.get_river_index(river_id)
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
                # if file is CF compliant, you can filter by date
                qout_nc.write_flows_to_csv(
                    '/timeseries/Qout_daily_date_filter.csv',
                    river_id=river_id,
                    daily=True,
                    date_search_start=datetime(2002, 8, 31),
                    date_search_end=datetime(2002, 9, 15),
                    filter_mode="max"
                )
        """
        if river_id is not None:
            river_index = self.get_river_index(river_id)
        elif river_id is None and river_index is None:
            raise ValueError("Need reach id or reach index ...")

        # analyze and write
        if self.is_time_variable_valid() or self._is_legacy_time_valid():
            qout_df = self.get_qout_index(river_index,
                                          date_search_start=date_search_start,
                                          date_search_end=date_search_end,
                                          daily=daily,
                                          filter_mode=filter_mode,
                                          as_dataframe=True)

            qout_df.to_csv(path_to_output_file, header=False)

        else:
            log("Valid time variable not found. Printing values only ...",
                "WARNING")
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
                                             filter_mode="mean"):
        """
        Write out RAPID output to GSSHA WMS time series xys file.

        Parameters
        ----------
        path_to_output_file: str
            Path to the output xys file.
        series_name: str
            The name for the series.
        series_id: int
            The ID to give the series.
        river_index: :obj:`datetime.datetime`, optional
            This is the index of the river in the file you want the
            streamflow for.
        river_id: :obj:`datetime.datetime`, optional
            This is the river ID that you want the streamflow for.
        date_search_start: :obj:`datetime.datetime`, optional
            This is a datetime object with the date of the minimum date for
            starting.
        date_search_end: :obj:`datetime.datetime`, optional
            This is a datetime object with the date of the maximum date for
            ending.
        daily: bool, optional
            If True and the file is CF-Compliant, write out daily flows.
        filter_mode: str, optional
            You can get the daily average "mean" or the maximum "max".
            Defauls is "mean".


        Example writing entire time series to file:

        .. code:: python

            from RAPIDpy import RAPIDDataset

            river_id = 3624735
            path_to_rapid_qout = '/path/to/Qout.nc'

            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                qout_nc.write_flows_to_gssha_time_series_xys(
                    '/timeseries/Qout_{0}.xys'.format(river_id),
                    series_name="RAPID_TO_GSSHA_{0}".format(river_id),
                    series_id=34,
                    river_id=river_id)


        Example writing entire time series as daily average to file:

        .. code:: python

            from RAPIDpy import RAPIDDataset

            river_id = 3624735
            path_to_rapid_qout = '/path/to/Qout.nc'

            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                # NOTE: Getting the river index is not necessary
                # this is just an example of how to use this
                river_index = qout_nc.get_river_index(river_id)

                # if file is CF compliant, you can write out daily average
                qout_nc.write_flows_to_gssha_time_series_xys(
                    '/timeseries/Qout_daily.xys',
                    series_name="RAPID_TO_GSSHA_{0}".format(river_id),
                    series_id=34,
                    river_index=river_index,
                    daily=True)


        Example writing subset of time series as daily maximum to file:

        .. code:: python

            from datetime import datetime
            from RAPIDpy import RAPIDDataset

            river_id = 3624735
            path_to_rapid_qout = '/path/to/Qout.nc'

            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                # NOTE: Getting the river index is not necessary
                # this is just an example of how to use this
                river_index = qout_nc.get_river_index(river_id)

                # if file is CF compliant, you can filter by date and
                # get daily values
                qout_nc.write_flows_to_gssha_time_series_xys(
                    '/timeseries/Qout_daily_date_filter.xys',
                    series_name="RAPID_TO_GSSHA_{0}".format(river_id),
                    series_id=34,
                    river_index=river_index,
                    date_search_start=datetime(2002, 8, 31),
                    date_search_end=datetime(2002, 9, 15),
                    daily=True,
                    filter_mode="max")

        """
        if river_id is not None:
            river_index = self.get_river_index(river_id)
        elif river_id is None and river_index is None:
            raise ValueError(" Need reach id or reach index ...")

        self.raise_time_valid()

        # analyze and write
        qout_df = self.get_qout_index(river_index,
                                      date_search_start=date_search_start,
                                      date_search_end=date_search_end,
                                      daily=daily,
                                      filter_mode=filter_mode,
                                      as_dataframe=True)

        with open_csv(path_to_output_file, 'w') as out_ts:
            out_ts.write("XYS {0} {1} \"{2}\"\r\n".format(series_id,
                                                          len(qout_df.index),
                                                          series_name))
            for index, pd_row in qout_df.iterrows():
                date_str = index.strftime("%m/%d/%Y %I:%M:%S %p")
                out_ts.write("\"{0}\" {1:.5f}\n".format(date_str,
                                                        pd_row[0]))

    def write_flows_to_gssha_time_series_ihg(self,
                                             path_to_output_file,
                                             connection_list_file,
                                             date_search_start=None,
                                             date_search_end=None,
                                             daily=False,
                                             filter_mode="mean"):
        # pylint: disable=line-too-long
        """
        Write out RAPID output to GSSHA time series ihg file

        .. note:: See: http://www.gsshawiki.com/Surface_Water_Routing:Introducing_Dischage/Constituent_Hydrographs

        .. note:: GSSHA project card is CHAN_POINT_INPUT

        Parameters
        ----------
        path_to_output_file: str
            Path to the output xys file.
        connection_list_file: str
            CSV file with link_id, node_id, baseflow, and rapid_rivid header
            and rows with data.
        date_search_start: :obj:`datetime.datetime`, optional
            This is a datetime object with the date of the minimum date
            for starting.
        date_search_end: :obj:`datetime.datetime`, optional
            This is a datetime object with the date of the maximum date
            for ending.
        daily: bool, optional
            If True and the file is CF-Compliant, write out daily flows.
        filter_mode: str, optional
            You can get the daily average "mean" or the maximum "max".
            Defauls is "mean".


        Example connection list file::

            link_id, node_id, baseflow, rapid_rivid
            599, 1, 0.0, 80968
            603, 1, 0.0, 80967


        Example writing entire time series to file:

        .. code:: python

            from RAPIDpy import RAPIDDataset

            path_to_rapid_qout = '/path/to/Qout.nc'
            connection_list_file = '/path/to/connection_list_file.csv'

            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                #for writing entire time series to file
                qout_nc.write_flows_to_gssha_time_series_ihg(
                    '/timeseries/Qout_3624735.ihg',
                    connection_list_file)


        Example writing entire time series as daily average to file:

        .. code:: python

            from RAPIDpy import RAPIDDataset

            path_to_rapid_qout = '/path/to/Qout.nc'
            connection_list_file = '/path/to/connection_list_file.csv'

            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                # if file is CF compliant, you can write out daily average
                qout_nc.write_flows_to_gssha_time_series_ihg(
                    '/timeseries/Qout_3624735.ihg',
                    connection_list_file,
                    daily=True)


        Example writing subset of time series as daily maximum to file:

        .. code:: python

            from datetime import datetime
            from RAPIDpy import RAPIDDataset

            path_to_rapid_qout = '/path/to/Qout.nc'
            connection_list_file = '/path/to/connection_list_file.csv'

            with RAPIDDataset(path_to_rapid_qout) as qout_nc:
                # if file is CF compliant, you can filter by
                # date and get daily values
                qout_nc.write_flows_to_gssha_time_series_ihg(
                    '/timeseries/Qout_daily_date_filter.ihg',
                    connection_list_file,
                    date_search_start=datetime(2002, 8, 31),
                    date_search_end=datetime(2002, 9, 15),
                    daily=True,
                    filter_mode="max")
        """  # noqa
        self.raise_time_valid()

        # analyze and write
        with open_csv(path_to_output_file, 'w') as out_ts:
            # HEADER SECTION EXAMPLE:
            # NUMPT 3
            # POINT 1 599 0.0
            # POINT 1 603 0.0
            # POINT 1 605 0.0

            connection_list = np.loadtxt(connection_list_file,
                                         skiprows=1, ndmin=1,
                                         delimiter=',',
                                         usecols=(0, 1, 2, 3),
                                         dtype={'names': ('link_id',
                                                          'node_id',
                                                          'baseflow',
                                                          'rapid_rivid'),
                                                'formats': ('i8', 'i8',
                                                            'f4', 'i8')
                                                },
                                         )

            out_ts.write("NUMPT {0}\n".format(connection_list.size))

            river_idx_list = []
            for connection in connection_list:
                out_ts.write("POINT {0} {1} {2}\n"
                             "".format(connection['node_id'],
                                       connection['link_id'],
                                       connection['baseflow'],
                                       ),
                             )
                river_idx_list.append(
                    self.get_river_index(connection['rapid_rivid'])
                    )

            # INFLOW SECTION EXAMPLE:
            # NRPDS 54
            # INPUT 2002 01 01 00 00 15.551210 12.765090 0.000000
            # INPUT 2002 01 02 00 00 15.480830 12.765090 0.000000
            # INPUT 2002 01 03 00 00 16.078910 12.765090 0.000000
            # ...
            qout_df = self.get_qout_index(
                river_idx_list,
                date_search_start=date_search_start,
                date_search_end=date_search_end,
                daily=daily,
                filter_mode=filter_mode,
                as_dataframe=True)

            out_ts.write("NRPDS {0}\n".format(len(qout_df.index)))

            for index, pd_row in qout_df.iterrows():
                date_str = index.strftime("%Y %m %d %H %M")
                qout_str = " ".join(["{0:.5f}".format(pd_row[column])
                                     for column in qout_df])
                out_ts.write("INPUT {0} {1}\n".format(date_str, qout_str))
