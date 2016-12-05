# -*- coding: utf-8 -*-
##
##  CreateInflowFileFromGriddedRunoff.py
##  RAPIDpy
##
##  Created by Alan D Snow.
##  Copyright © 2016 Alan D Snow. All rights reserved.
##  BSD 3-Clause

import csv
from datetime import datetime
import netCDF4 as NET
import numpy as np
import os
from pytz import utc

#local
from ..helper_functions import open_csv

class CreateInflowFileFromGriddedRunoff(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.header_wt = ['rivid', 'area_sqm', 'lon_index', 'lat_index', 'npoints']

    def readInWeightTable(self, in_weight_table):
        """
        Read in weight table
        """
        print("Reading the weight table...")
        with open_csv(in_weight_table, "r") as csvfile:
            reader = csv.reader(csvfile)
            header_row = next(reader)
            #check number of columns in the weight table
            if len(header_row) < len(self.header_wt):
                raise Exception(self.errorMessages[4])
            #check header
            if header_row[1:len(self.header_wt)] != self.header_wt[1:]:
                raise Exception(self.errorMessages[5])

        self.dict_list = np.loadtxt(in_weight_table, 
                                    delimiter=",", 
                                    usecols=(0,1,2,3,4),
                                    skiprows=1,
                                    dtype={'names': (self.header_wt[0], self.header_wt[1], self.header_wt[2], self.header_wt[3], self.header_wt[4]),
                                           'formats': ('i8', 'f8', 'i8', 'i8', 'i8')},
                                    )
                                    
        self.count = self.dict_list.shape[0]
        self.size_streamID = len(np.unique(np.array(self.dict_list[self.header_wt[0]], dtype=np.int32)))
        
    def _write_lat_lon(self, data_out_nc, rivid_lat_lon_z_file):
        """Add latitude and longitude each netCDF feature
        Lookup table is a CSV file with rivid, Lat, Lon, columns.
        Columns must be in that order and these must be the first three columns.
        """
        #only add if user adds
        if rivid_lat_lon_z_file and os.path.exists(rivid_lat_lon_z_file):
            #get list of COMIDS
            lookup_table = np.loadtxt(rivid_lat_lon_z_file, 
                                      delimiter=",", 
                                      usecols=(0,1,2),
                                      skiprows=1,
                                      dtype={'names': ('rivid', 'lat', 'lon'),
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
                    lookup_index = np.where(lookup_table['rivid'] == nc_rivid)[0][0]
                except Exception:
                    raise Exception('rivid {0} misssing in comid_lat_lon_z file'.format(nc_rivid))
        
                lat = float(lookup_table['lat'][lookup_index])
                lats[nc_index] = lat
                if (lat_min) is None or lat < lat_min:
                    lat_min = lat
                if (lat_max) is None or lat > lat_max:
                    lat_max = lat
        
                lon = float(lookup_table['lon'][lookup_index])
                lons[nc_index] = lon
                if (lon_min) is None or lon < lon_min:
                    lon_min = lon
                if (lon_max) is None or lon > lon_max:
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
                                 out_nc, #file generated for inflows
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

        # Create output inflow netcdf data
        print("Generating inflow file ...")
        data_out_nc = NET.Dataset(out_nc, "w", format="NETCDF3_CLASSIC")
        rivid_list = np.loadtxt(in_rapid_connect_file, 
                                delimiter=",",
                                ndmin=1, 
                                usecols=(0,), 
                                dtype=int)
        #create dimensions
        data_out_nc.createDimension('time', number_of_timesteps)
        data_out_nc.createDimension('rivid', len(rivid_list))
        data_out_nc.createDimension('nv', 2)
        #create variables
        #m3_riv
        m3_riv_var = data_out_nc.createVariable('m3_riv', 'f4', 
                                                ('time', 'rivid'),
                                                fill_value=0)
        m3_riv_var.long_name = 'accumulated external water volume inflow upstream of each river reach'
        m3_riv_var.units = 'm3'
        m3_riv_var.coordinates = 'lon lat'
        m3_riv_var.grid_mapping = 'crs'
        m3_riv_var.cell_methods = "time: sum"
        data_out_nc.close()
        
        try:
            data_out_nc = NET.Dataset(out_nc, "a", format="NETCDF3_CLASSIC")
            #rivid
            rivid_var = data_out_nc.createVariable('rivid', 'i4', 
                                                  ('rivid',))
            rivid_var.long_name = 'unique identifier for each river reach'
            rivid_var.units = '1'
            rivid_var.cf_role = 'timeseries_id'
            
            rivid_var[:] = rivid_list
    
            #time
            time_var = data_out_nc.createVariable('time', 'i4',
                                                  ('time',))
            time_var.long_name = 'time'
            time_var.standard_name = 'time'
            time_var.units = 'seconds since 1970-01-01 00:00:00+00:00'
            time_var.axis = 'T'
            time_var.calendar = 'gregorian'
            time_var.bounds = 'time_bnds'
    
            initial_time_seconds = (start_datetime_utc.replace(tzinfo=utc)- \
                                    datetime(1970,1,1, tzinfo=utc)).total_seconds()
            final_time_seconds = initial_time_seconds + number_of_timesteps*simulation_time_step_seconds
            time_array = np.arange(initial_time_seconds, final_time_seconds, simulation_time_step_seconds)
            time_var[:] = time_array
    
            #time_bnds
            time_bnds_var = data_out_nc.createVariable('time_bnds', 'i4',
                                                        ('time', 'nv',))
            for time_index, time_element in enumerate(time_array):
                time_bnds_var[time_index, 0] = time_element
                time_bnds_var[time_index, 1] = time_element+simulation_time_step_seconds
    
            #longitude
            lon_var = data_out_nc.createVariable('lon', 'f8', ('rivid',),
                                                fill_value=-9999.0)
            lon_var.long_name = 'longitude of a point related to each river reach'
            lon_var.standard_name = 'longitude'
            lon_var.units = 'degrees_east'
            lon_var.axis = 'X'
    
            #latitude
            lat_var = data_out_nc.createVariable('lat', 'f8', ('rivid',),
                                                fill_value=-9999.0)
            lat_var.long_name = 'latitude of a point related to each river reach'
            lat_var.standard_name = 'latitude'
            lat_var.units = 'degrees_north'
            lat_var.axis = 'Y'
                                       
            crs_var = data_out_nc.createVariable('crs', 'i4')
            crs_var.grid_mapping_name = 'latitude_longitude'
            crs_var.epsg_code = 'EPSG:4326'  # WGS 84
            crs_var.semi_major_axis = 6378137.0
            crs_var.inverse_flattening = 298.257223563
            
            #add global attributes
            data_out_nc.Conventions = 'CF-1.6'
            data_out_nc.title = 'RAPID Inflow from {0}'.format(land_surface_model_description)
            data_out_nc.history = 'date_created: {0}'.format(datetime.utcnow().replace(tzinfo=utc))
            data_out_nc.featureType = 'timeSeries'
            data_out_nc.institution = modeling_institution
            
            #write lat lon data
            self._write_lat_lon(data_out_nc, in_rivid_lat_lon_z_file)
            
            #close file
            data_out_nc.close()
        except RuntimeError:
            print("File size too big to add data beforehand. Performing conversion after ...")
            pass