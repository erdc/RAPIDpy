# -*- coding: utf-8 -*-
"""
   CreateInflowFileFromGriddedRunoff.py
   RAPIDpy

   Created by Alan D. Snow, 2016
   License: BSD-3-Clause
"""
from abc import abstractmethod
import csv
from datetime import datetime
import os

from netCDF4 import Dataset
import numpy as np
from pytz import utc
from past.builtins import xrange  # pylint: disable=redefined-builtin

# local
from ..helper_functions import open_csv


class CreateInflowFileFromGriddedRunoff(object):
    """Create Inflow File From Gridded Runoff

    Base class for creating RAPID NetCDF input
    of water inflow based on land surface model
    runoff and previously created weight table.
    """
    land_surface_model_name = "land surface model"
    header_wt = ['rivid', 'area_sqm', 'lon_index', 'lat_index', 'npoints']
    runoff_vars = []

    def __init__(self):
        self.dict_list = []
        self.count = 0
        self.size_stream_id = 0
        self.simulation_time_step_seconds = 0
        self.error_messages = [
            "Missing Variable 'time'",
            "Incorrect dimensions in the input {} runoff file."
            .format(self.land_surface_model_name),
            "Incorrect variables in the input {} runoff file."
            .format(self.land_surface_model_name),
            "Incorrect time variable in the input {} runoff file"
            .format(self.land_surface_model_name),
            "Incorrect number of columns in the weight table",
            "No or incorrect header in the weight table",
            "Incorrect sequence of rows in the weight table"
        ]

    def read_in_weight_table(self, in_weight_table):
        """
        Read in weight table
        """
        print("Reading the weight table...")
        with open_csv(in_weight_table, "r") as csvfile:
            reader = csv.reader(csvfile)
            header_row = next(reader)
            # check number of columns in the weight table
            if len(header_row) < len(self.header_wt):
                raise Exception(self.error_messages[4])
            # check header
            if header_row[1:len(self.header_wt)] != self.header_wt[1:]:
                raise Exception(self.error_messages[5])

        self.dict_list = \
            np.loadtxt(
                in_weight_table,
                delimiter=",",
                usecols=(0, 1, 2, 3, 4),
                skiprows=1,
                dtype={
                    'names': (self.header_wt[0],
                              self.header_wt[1],
                              self.header_wt[2],
                              self.header_wt[3],
                              self.header_wt[4]),
                    'formats': ('i8', 'f8', 'i8', 'i8', 'i8')
                },
            )

        self.count = self.dict_list.shape[0]
        self.size_stream_id = \
            len(np.unique(np.array(self.dict_list[self.header_wt[0]],
                                   dtype=np.int32)))

    @staticmethod
    def _write_lat_lon(data_out_nc, rivid_lat_lon_z_file):
        """Add latitude and longitude each netCDF feature
        Lookup table is a CSV file with rivid, Lat, Lon, columns.
        Columns must be in that order and these must be the first
        three columns.
        """
        # only add if user adds
        if rivid_lat_lon_z_file and os.path.exists(rivid_lat_lon_z_file):
            # get list of COMIDS
            lookup_table = np.loadtxt(
                rivid_lat_lon_z_file,
                delimiter=",",
                usecols=(0, 1, 2),
                skiprows=1,
                dtype={
                    'names': ('rivid', 'lat', 'lon'),
                    'formats': ('i8', 'f8', 'f8'),
                },
            )

            # Get relevant arrays while we update them
            nc_rivids = data_out_nc.variables['rivid'][:]
            lats = data_out_nc.variables['lat'][:]
            lons = data_out_nc.variables['lon'][:]

            lat_min = None
            lat_max = None
            lon_min = None
            lon_max = None

            # Process each row in the lookup table
            for nc_index, nc_rivid in enumerate(nc_rivids):
                try:
                    lookup_index = \
                        np.where(lookup_table['rivid'] == nc_rivid)[0][0]
                except Exception:
                    raise Exception('rivid {0} misssing in '
                                    'comid_lat_lon_z file'.format(nc_rivid))

                lat = float(lookup_table['lat'][lookup_index])
                lats[nc_index] = lat
                if lat_min is None or lat < lat_min:
                    lat_min = lat
                if lat_max is None or lat > lat_max:
                    lat_max = lat

                lon = float(lookup_table['lon'][lookup_index])
                lons[nc_index] = lon
                if lon_min is None or lon < lon_min:
                    lon_min = lon
                if lon_max is None or lon > lon_max:
                    lon_max = lon

            # Overwrite netCDF variable values
            data_out_nc.variables['lat'][:] = lats
            data_out_nc.variables['lon'][:] = lons

            # Update metadata
            if lat_min is not None:
                data_out_nc.geospatial_lat_min = lat_min
            if lat_max is not None:
                data_out_nc.geospatial_lat_max = lat_max
            if lon_min is not None:
                data_out_nc.geospatial_lon_min = lon_min
            if lon_max is not None:
                data_out_nc.geospatial_lon_max = lon_max
        else:
            print('No comid_lat_lon_z file. Not adding values ...')

    def generateOutputInflowFile(self,
                                 out_nc,
                                 start_datetime_utc,
                                 number_of_timesteps,
                                 simulation_time_step_seconds,
                                 in_rapid_connect_file,
                                 in_rivid_lat_lon_z_file,
                                 land_surface_model_description,
                                 modeling_institution
                                 ):
        """
        Generate inflow file for RAPID
        """
        self.simulation_time_step_seconds = simulation_time_step_seconds

        # Create output inflow netcdf data
        print("Generating inflow file ...")
        data_out_nc = Dataset(out_nc, "w", format="NETCDF3_CLASSIC")
        rivid_list = np.loadtxt(in_rapid_connect_file,
                                delimiter=",",
                                ndmin=1,
                                usecols=(0,),
                                dtype=int)
        # create dimensions
        data_out_nc.createDimension('time', number_of_timesteps)
        data_out_nc.createDimension('rivid', len(rivid_list))
        data_out_nc.createDimension('nv', 2)
        # create variables
        # m3_riv
        m3_riv_var = data_out_nc.createVariable('m3_riv', 'f4',
                                                ('time', 'rivid'),
                                                fill_value=0)
        m3_riv_var.long_name = 'accumulated external water volume ' \
                               'inflow upstream of each river reach'
        m3_riv_var.units = 'm3'
        m3_riv_var.coordinates = 'lon lat'
        m3_riv_var.grid_mapping = 'crs'
        m3_riv_var.cell_methods = "time: sum"
        data_out_nc.close()

        try:
            data_out_nc = Dataset(out_nc, "a", format="NETCDF3_CLASSIC")
            # rivid
            rivid_var = data_out_nc.createVariable('rivid', 'i4',
                                                   ('rivid',))
            rivid_var.long_name = 'unique identifier for each river reach'
            rivid_var.units = '1'
            rivid_var.cf_role = 'timeseries_id'
            rivid_var[:] = rivid_list

            # time
            time_var = data_out_nc.createVariable('time', 'i4',
                                                  ('time',))
            time_var.long_name = 'time'
            time_var.standard_name = 'time'
            time_var.units = 'seconds since 1970-01-01 00:00:00+00:00'
            time_var.axis = 'T'
            time_var.calendar = 'gregorian'
            time_var.bounds = 'time_bnds'

            initial_time_seconds = \
                (start_datetime_utc.replace(tzinfo=utc) -
                 datetime(1970, 1, 1, tzinfo=utc)).total_seconds()
            final_time_seconds = \
                initial_time_seconds + number_of_timesteps\
                * simulation_time_step_seconds
            time_array = np.arange(initial_time_seconds, final_time_seconds,
                                   simulation_time_step_seconds)
            time_var[:] = time_array

            # time_bnds
            time_bnds_var = data_out_nc.createVariable('time_bnds', 'i4',
                                                       ('time', 'nv',))
            for time_index, time_element in enumerate(time_array):
                time_bnds_var[time_index, 0] = time_element
                time_bnds_var[time_index, 1] = \
                    time_element + simulation_time_step_seconds

            # longitude
            lon_var = data_out_nc.createVariable('lon', 'f8', ('rivid',),
                                                 fill_value=-9999.0)
            lon_var.long_name = \
                'longitude of a point related to each river reach'
            lon_var.standard_name = 'longitude'
            lon_var.units = 'degrees_east'
            lon_var.axis = 'X'

            # latitude
            lat_var = data_out_nc.createVariable('lat', 'f8', ('rivid',),
                                                 fill_value=-9999.0)
            lat_var.long_name = \
                'latitude of a point related to each river reach'
            lat_var.standard_name = 'latitude'
            lat_var.units = 'degrees_north'
            lat_var.axis = 'Y'

            crs_var = data_out_nc.createVariable('crs', 'i4')
            crs_var.grid_mapping_name = 'latitude_longitude'
            crs_var.epsg_code = 'EPSG:4326'  # WGS 84
            crs_var.semi_major_axis = 6378137.0
            crs_var.inverse_flattening = 298.257223563

            # add global attributes
            data_out_nc.Conventions = 'CF-1.6'
            data_out_nc.title = 'RAPID Inflow from {0}'\
                .format(land_surface_model_description)
            data_out_nc.history = 'date_created: {0}'\
                .format(datetime.utcnow().replace(tzinfo=utc))
            data_out_nc.featureType = 'timeSeries'
            data_out_nc.institution = modeling_institution

            # write lat lon data
            self._write_lat_lon(data_out_nc, in_rivid_lat_lon_z_file)

            # close file
            data_out_nc.close()
        except RuntimeError:
            print("File size too big to add data beforehand."
                  " Performing conversion after ...")

    def get_conversion_factor(self, in_nc, num_nc_files):
        """get conversion_factor"""
        data_in_nc = Dataset(in_nc)

        # convert from kg/m^2 (i.e. mm) to m
        conversion_factor = 0.001

        # ECMWF units are in m
        if data_in_nc.variables[self.runoff_vars[0]] \
                .getncattr("units") == "m":
            conversion_factor = 1

        # ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/GLDAS_V1/README.GLDAS.pdf
        if "s" in data_in_nc.variables[self.runoff_vars[0]] \
                .getncattr("units"):
            # that means kg/m^2/s in GLDAS v1 that is 3-hr avg,
            # so multiply by 3 hr (ex. 3*3600). Assumed same
            # for others (ex. 1*3600).
            # If combining files, need to take average of these,
            # so divide by number of files
            conversion_factor *= \
                self.simulation_time_step_seconds / \
                num_nc_files
        data_in_nc.close()

        return conversion_factor

    @abstractmethod
    def data_validation(self, in_nc):
        """Check the necessary dimensions and variables in the input
        netcdf data"""
        pass

    def execute(self, nc_file_list, index_list, in_weight_table,
                out_nc, grid_type, mp_lock):

        """The source code of the tool."""
        if not os.path.exists(out_nc):
            raise Exception("Outfile has not been created. "
                            "You need to run: generateOutputInflowFile "
                            "function ...")

        if len(nc_file_list) != len(index_list):
            raise Exception("ERROR: Number of runoff files not equal to "
                            "number of indices ...")

        demo_file_list = nc_file_list[0]
        if not isinstance(nc_file_list[0], list):
            demo_file_list = [demo_file_list]

        self.data_validation(demo_file_list[0])
        self.read_in_weight_table(in_weight_table)

        conversion_factor = self.get_conversion_factor(demo_file_list[0],
                                                       len(demo_file_list))

        # get indices of subset of data
        lon_ind_all = [int(i) for i in self.dict_list[self.header_wt[2]]]
        lat_ind_all = [int(j) for j in self.dict_list[self.header_wt[3]]]

        # Obtain a subset of  runoff data based on the indices in the
        # weight table
        min_lon_ind_all = min(lon_ind_all)
        max_lon_ind_all = max(lon_ind_all)
        min_lat_ind_all = min(lat_ind_all)
        max_lat_ind_all = max(lat_ind_all)
        lon_slice = slice(min_lon_ind_all, max_lon_ind_all + 1)
        lat_slice = slice(min_lat_ind_all, max_lat_ind_all + 1)
        index_new = []

        # combine inflow data
        for nc_file_array_index, nc_file_array in enumerate(nc_file_list):

            index = index_list[nc_file_array_index]

            if not isinstance(nc_file_array, list):
                nc_file_array = [nc_file_array]

            data_subset_all = None
            for nc_file in nc_file_array:
                # Validate the netcdf dataset
                self.data_validation(nc_file)

                # Read the netcdf dataset
                data_in_nc = Dataset(nc_file)

                # Calculate water inflows
                runoff_dimension_size = \
                    len(data_in_nc.variables[self.runoff_vars[0]].dimensions)
                if runoff_dimension_size == 2:
                    # obtain subset of surface and subsurface runoff
                    data_subset_runoff = \
                        data_in_nc.variables[self.runoff_vars[0]][
                            lat_slice, lon_slice]
                    for var_name in self.runoff_vars[1:]:
                        data_subset_runoff += \
                            data_in_nc.variables[var_name][
                                lat_slice, lon_slice]

                    # get runoff dims
                    len_time_subset = 1
                    len_lat_subset = data_subset_runoff.shape[0]
                    len_lon_subset = data_subset_runoff.shape[1]

                    # reshape the runoff
                    data_subset_runoff = data_subset_runoff.reshape(
                        len_lat_subset * len_lon_subset)

                elif runoff_dimension_size == 3:
                    # obtain subset of surface and subsurface runoff
                    data_subset_runoff = \
                        data_in_nc.variables[self.runoff_vars[0]][
                            :, lat_slice, lon_slice]
                    for var_name in self.runoff_vars[1:]:
                        data_subset_runoff += \
                            data_in_nc.variables[var_name][
                                :, lat_slice, lon_slice]

                    # get runoff dims
                    len_time_subset = data_subset_runoff.shape[0]
                    len_lat_subset = data_subset_runoff.shape[1]
                    len_lon_subset = data_subset_runoff.shape[2]
                    # reshape the runoff
                    data_subset_runoff = \
                        data_subset_runoff.reshape(
                            len_time_subset,
                            (len_lat_subset * len_lon_subset))

                data_in_nc.close()

                if not index_new:
                    # compute new indices based on the data_subset_surface
                    for r in range(0, self.count):
                        ind_lat_orig = lat_ind_all[r]
                        ind_lon_orig = lon_ind_all[r]
                        index_new.append(
                            (ind_lat_orig - min_lat_ind_all) * len_lon_subset
                            + (ind_lon_orig - min_lon_ind_all))

                # obtain a new subset of data
                if runoff_dimension_size == 2:
                    data_subset_new = data_subset_runoff[index_new]
                elif runoff_dimension_size == 3:
                    data_subset_new = data_subset_runoff[:, index_new]

                # FILTER DATA
                try:
                    # set masked values to zero
                    data_subset_new = data_subset_new.filled(fill_value=0)
                except AttributeError:
                    pass
                # set negative values to zero
                data_subset_new[data_subset_new < 0] = 0

                # combine data
                if data_subset_all is None:
                    data_subset_all = data_subset_new
                else:
                    data_subset_all = np.add(data_subset_all, data_subset_new)

            if runoff_dimension_size == 3 and len_time_subset > 1:
                inflow_data = np.zeros((len_time_subset, self.size_stream_id))
            else:
                inflow_data = np.zeros(self.size_stream_id)

            pointer = 0
            for stream_index in xrange(self.size_stream_id):
                npoints = int(self.dict_list[self.header_wt[4]][pointer])
                # Check if all npoints points correspond to the same streamID
                if len(set(self.dict_list[self.header_wt[0]][
                           pointer: (pointer + npoints)])) != 1:
                    print("ROW INDEX {0}".format(pointer))
                    print("COMID {0}".format(
                        self.dict_list[self.header_wt[0]][pointer]))
                    raise Exception(self.error_messages[2])

                area_sqm_npoints = \
                    np.array([float(k) for k in
                              self.dict_list[self.header_wt[1]][
                              pointer: (pointer + npoints)]])

                # assume data is incremental
                if runoff_dimension_size == 3:
                    data_goal = data_subset_all[:, pointer:(pointer + npoints)]
                else:
                    data_goal = data_subset_all[pointer:(pointer + npoints)]

                if grid_type == 't255':
                    # A) ERA Interim Low Res (T255) - data is cumulative
                    # from time 3/6/9/12
                    # (time zero not included, so assumed to be zero)
                    ro_first_half = \
                        np.concatenate([data_goal[0:1, ],
                                        np.subtract(data_goal[1:4, ],
                                                    data_goal[0:3, ])])
                    # from time 15/18/21/24
                    # (time restarts at time 12, assumed to be zero)
                    ro_second_half = \
                        np.concatenate([data_goal[4:5, ],
                                        np.subtract(data_goal[5:, ],
                                                    data_goal[4:7, ])])
                    ro_stream = \
                        np.multiply(
                            np.concatenate([ro_first_half, ro_second_half]),
                            area_sqm_npoints)

                else:
                    ro_stream = data_goal * area_sqm_npoints * \
                                conversion_factor

                # filter nan
                ro_stream[np.isnan(ro_stream)] = 0

                if ro_stream.any():
                    if runoff_dimension_size == 3 and len_time_subset > 1:
                        inflow_data[:, stream_index] = ro_stream.sum(axis=1)
                    else:
                        inflow_data[stream_index] = ro_stream.sum()

                pointer += npoints

            # only one process is allowed to write at a time to netcdf file
            mp_lock.acquire()
            data_out_nc = Dataset(out_nc, "a", format="NETCDF3_CLASSIC")
            if runoff_dimension_size == 3 and len_time_subset > 1:
                data_out_nc.variables['m3_riv'][
                    index*len_time_subset:(index+1)*len_time_subset, :] = \
                    inflow_data
            else:
                data_out_nc.variables['m3_riv'][index] = inflow_data
            data_out_nc.close()
            mp_lock.release()
