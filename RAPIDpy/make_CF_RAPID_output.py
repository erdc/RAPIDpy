# -*- coding: utf-8 -*-
##
##  make_CF_output.py
##  RAPIDpy
##
##  Created by Alan D Snow & Tim Whitaker, 2015.
##  Copyright © 2015 Alan D Snow & Tim Whitaker. All rights reserved.
##

"""Copies data from RAPID netCDF output to a CF-compliant netCDF file.
Code originated from Tim Whitaker at University of Texas. The code was
modified by Alan Snow at US Army ERDC.
 
Remarks:
    A new netCDF file is created with data from RAPID [1] simulation model
    output. The result follows CF conventions [2] with additional metadata
    prescribed by the NODC timeSeries Orthogonal template [3] for time series
    at discrete point feature locations.

    This script was created for the National Flood Interoperability Experiment,
    and so metadata in the result reflects that.

Requires:
    netcdf4-python - https://github.com/Unidata/netcdf4-python

Inputs:
    Lookup CSV table with COMID, Lat, Lon, and Elev_m columns. Columns must
    be in that order and these must be the first four columns. The order of
    COMIDs in the table must match the order of features in the netCDF file.


///////////////////////////////////////////////////
netcdf result_2014100520141101 {
dimensions:
    Time = UNLIMITED ; // (224 currently)
    COMID = 61818 ;
variables:
    float Qout(Time, COMID) ;
///////////////////////////////////////////////////

Outputs:
    CF-compliant netCDF file of RAPID results, named with original filename
    with "_CF" appended to the filename. File is written to 'output' folder.

    Input netCDF file is archived or deleted, based on 'archive' config
    parameter.

Usage:
    import the script, e.g., import ConvertRAPIDOutputToCF as cf.


References:
    [1] http://rapid-hub.org/
    [2] http://cfconventions.org/
    [3] http://www.nodc.noaa.gov/data/formats/netcdf/v1.1/
"""


from datetime import datetime, timedelta
import os
from netCDF4 import Dataset
import numpy as np

#local
from helper_functions import csv_to_list, remove_files

def log(message, severity, print_debug=True):
    """Logs, prints, or raises a message.

    Arguments:
        message -- message to report
        severity -- string of one of these values:
            CRITICAL|ERROR|WARNING|INFO|DEBUG
    """

    print_me = ['WARNING', 'INFO', 'DEBUG']
    if severity in print_me:
        if severity == 'DEBUG':
            if print_debug:
                print severity, message
        else:
                print severity, message
    else:
        raise Exception("%s: %s" % (severity, message))


class ConvertRAPIDOutputToCF(object):
    """
    Class to convert RAPID output to be CF compliant        
    """
    def __init__(self, rapid_output_file, #location of timeseries output file(s)
                       start_datetime, #time of the start of the simulation time
                       time_step, #time step(s) of simulation in seconds
                       qinit_file="", #RAPID qinit file
                       comid_lat_lon_z_file="", #path to comid_lat_lon_z file
                       rapid_connect_file="", #path to RAPID connect file
                       project_name="Default RAPID Project", #name of your project
                       output_id_dim_name='rivid', #name of ID dimension in output file, typically COMID or FEATUREID
                       output_flow_var_name='Qout', #name of streamflow variable in output file, typically Qout or m3_riv
                       print_debug=False
                       ):

       if not isinstance(rapid_output_file, list): 
           self.rapid_output_file_list = [rapid_output_file]
       else:
           self.rapid_output_file_list = rapid_output_file
       self.start_datetime = start_datetime
       
       if not isinstance(time_step, list): 
           self.time_step_array = [time_step]
       else:
           self.time_step_array = time_step
           
       self.qinit_file = qinit_file
       self.comid_lat_lon_z_file = comid_lat_lon_z_file
       self.rapid_connect_file = rapid_connect_file
       self.project_name = project_name
       self.output_id_dim_name = output_id_dim_name
       self.output_flow_var_name = output_flow_var_name
       self.print_debug = print_debug
       self.cf_compliant_file = '%s_CF.nc' % os.path.splitext(self.rapid_output_file_list[0])[0]

    def _validate_raw_nc(self):
        """Checks that raw netCDF file has the right dimensions and variables.
    
        Arguments:
            nc -- netCDF dataset object representing raw RAPID output
    
        Returns:
            name of ID dimension,
            length of time dimension,
            name of flow variable
    
        Remarks: Raises exception if file doesn't validate.
        """

        self.raw_nc_list = []
        total_time_len = 0  
        initial_file = True
        id_dim_name_list = []
        id_len_list = []
        q_var_name_list = []
        self.time_len_array = []
        for rapid_output_file in self.rapid_output_file_list:
            raw_nc = Dataset(rapid_output_file)
            dims = raw_nc.dimensions
            if 'rivid' in dims:
                id_dim_name = 'rivid'
            elif 'COMID' in dims:
                id_dim_name = 'COMID'
            elif 'FEATUREID' in dims:
                id_dim_name = 'FEATUREID'
            else:
                msg = 'Could not find ID dimension. Looked for rivid, COMID and FEATUREID.'
                raise Exception(msg)
            id_dim_name_list.append(id_dim_name)

            id_len_list.append(len(dims[id_dim_name]))

            if 'time' not in dims:
                msg = 'Could not find time dimension. Looked for Time.'
                raise Exception(msg)
                
            time_len = len(dims['time'])
            if initial_file:
                time_len += 1 #add one for the first flow value RAPID
                            #does not include
                initial_file = False
            total_time_len += time_len
            self.time_len_array.append(time_len)
            
            variables = raw_nc.variables
        
            if 'Qout' in variables:
                q_var_name = 'Qout'
            elif 'm3_riv' in variables:
                q_var_name = 'm3_riv'
            else:
                log('Could not find flow variable. Looked for Qout and m3_riv.',
                    'ERROR')

            var_dims = variables[q_var_name].dimensions
            id_var_name = None
            if 'rivid' in var_dims:
                id_var_name = 'rivid'
            elif 'COMID' in var_dims:
                id_var_name = 'COMID'
            elif 'FEATUREID' in var_dims:
                id_var_name = 'FEATUREID'
    
            if id_var_name is not None and id_var_name != id_dim_name:
                msg = ('ID dimension name (' + id_dim_name + ') does not equal ID ' +
                       'variable name (' + id_var_name + ').')
                log(msg, 'WARNING')

            q_var_name_list.append(q_var_name)
            
            #make sure all id_dim_names same
            if not all(x == id_dim_name_list[0] for x in id_dim_name_list):
                log('ID dimension name not same for all files.', 'ERROR')
            #make sure all q var_names same       
            if not all(x == id_len_list[0] for x in id_len_list):
                log('ID dimension length not same for all files.', 'ERROR')
            #make sure all q var_names same       
            if not all(x == q_var_name_list[0] for x in q_var_name_list):
                log('Q variable name not same for all files.', 'ERROR')
            
            self.raw_nc_list.append(raw_nc)
            
        return id_dim_name, id_len_list[0], total_time_len, q_var_name


    def _initialize_output(self, time_len, id_len):
        """Creates netCDF file with CF dimensions and variables, but no data.
    
        Arguments:
            filename -- full path and filename for output netCDF file
            id_dim_name -- name of Id dimension and variable, e.g., COMID
            time_len -- (integer) length of time dimension (number of time steps)
            id_len -- (integer) length of Id dimension (number of time series)
            time_step_seconds -- (integer) number of seconds per time step
        """

        log('Initializing new file %s' % self.cf_compliant_file, 'INFO')
        
        self.cf_nc = Dataset(self.cf_compliant_file, 'w', format='NETCDF3_CLASSIC')
    
        # Create global attributes
        log('    globals', 'DEBUG', self.print_debug)
        self.cf_nc.featureType = 'timeSeries'
        self.cf_nc.Metadata_Conventions = 'Unidata Dataset Discovery v1.0'
        self.cf_nc.Conventions = 'CF-1.6'
        self.cf_nc.cdm_data_type = 'Station'
        self.cf_nc.nodc_template_version = (
            'NODC_NetCDF_TimeSeries_Orthogonal_Template_v1.1')
        self.cf_nc.standard_name_vocabulary = ('NetCDF Climate and Forecast (CF) ' +
                                          'Metadata Convention Standard Name ' +
                                          'Table v28')
        self.cf_nc.title = 'RAPID Result'
        self.cf_nc.summary = ("Results of RAPID river routing simulation. Each river " +
                         "reach (i.e., feature) is represented by a point " +
                         "feature at its midpoint, and is identified by the " +
                         "reach's unique NHDPlus COMID identifier.")
        self.cf_nc.time_coverage_resolution = 'point'
        self.cf_nc.geospatial_lat_min = 0.0
        self.cf_nc.geospatial_lat_max = 0.0
        self.cf_nc.geospatial_lat_units = 'degrees_north'
        self.cf_nc.geospatial_lat_resolution = 'midpoint of stream feature'
        self.cf_nc.geospatial_lon_min = 0.0
        self.cf_nc.geospatial_lon_max = 0.0
        self.cf_nc.geospatial_lon_units = 'degrees_east'
        self.cf_nc.geospatial_lon_resolution = 'midpoint of stream feature'
        self.cf_nc.geospatial_vertical_min = 0.0
        self.cf_nc.geospatial_vertical_max = 0.0
        self.cf_nc.geospatial_vertical_units = 'm'
        self.cf_nc.geospatial_vertical_resolution = 'midpoint of stream feature'
        self.cf_nc.geospatial_vertical_positive = 'up'
        self.cf_nc.project = self.project_name
        self.cf_nc.processing_level = 'Raw simulation result'
        self.cf_nc.keywords_vocabulary = ('NASA/Global Change Master Directory ' +
                                     '(GCMD) Earth Science Keywords. Version ' +
                                     '8.0.0.0.0')
        self.cf_nc.keywords = 'DISCHARGE/FLOW'
        self.cf_nc.comment = 'Result time step(s) (seconds): ' + str(self.time_step_array)
    
        timestamp = datetime.utcnow().isoformat() + 'Z'
        self.cf_nc.date_created = timestamp
        self.cf_nc.history = (timestamp + '; added time, lat, lon, z, crs variables; ' +
                         'added metadata to conform to NODC_NetCDF_TimeSeries_' +
                         'Orthogonal_Template_v1.1')
    
        # Create dimensions
        log('    dimming', 'DEBUG', self.print_debug)
        self.cf_nc.createDimension('time', time_len)
        self.cf_nc.createDimension(self.output_id_dim_name, id_len)
    
        # Create variables
        log('    timeSeries_var', 'DEBUG', self.print_debug)
        timeSeries_var = self.cf_nc.createVariable(self.output_id_dim_name, 'i4', 
                                                   (self.output_id_dim_name,))
        timeSeries_var.long_name = (
            'Unique NHDPlus COMID identifier for each river reach feature')
        timeSeries_var.cf_role = 'timeseries_id'
    
        log('    time_var', 'DEBUG', self.print_debug)
        time_var = self.cf_nc.createVariable('time', 'i4', ('time',))
        time_var.long_name = 'time'
        time_var.standard_name = 'time'
        time_var.units = 'seconds since 1970-01-01 00:00:00 0:00'
        time_var.axis = 'T'
    
        log('    lat_var', 'DEBUG', self.print_debug)
        lat_var = self.cf_nc.createVariable('lat', 'f8', (self.output_id_dim_name,),
                                       fill_value=-9999.0)
        lat_var.long_name = 'latitude'
        lat_var.standard_name = 'latitude'
        lat_var.units = 'degrees_north'
        lat_var.axis = 'Y'
    
        log('    lon_var', 'DEBUG', self.print_debug)
        lon_var = self.cf_nc.createVariable('lon', 'f8', (self.output_id_dim_name,),
                                       fill_value=-9999.0)
        lon_var.long_name = 'longitude'
        lon_var.standard_name = 'longitude'
        lon_var.units = 'degrees_east'
        lon_var.axis = 'X'
    
        log('    z_var', 'DEBUG', self.print_debug)
        z_var = self.cf_nc.createVariable('z', 'f8', (self.output_id_dim_name,),
                                     fill_value=-9999.0)
        z_var.long_name = ('Elevation referenced to the North American ' +
                           'Vertical Datum of 1988 (NAVD88)')
        z_var.standard_name = 'surface_altitude'
        z_var.units = 'm'
        z_var.axis = 'Z'
        z_var.positive = 'up'
    
        log('    crs_var', 'DEBUG', self.print_debug)
        crs_var = self.cf_nc.createVariable('crs', 'i4')
        crs_var.grid_mapping_name = 'latitude_longitude'
        crs_var.epsg_code = 'EPSG:4269'  # NAD83, which is what NHD uses.
        crs_var.semi_major_axis = 6378137.0
        crs_var.inverse_flattening = 298.257222101

    def _write_comid_lat_lon_z(self):
        """Add latitude, longitude, and z values for each netCDF feature
    
        Arguments:
            cf_nc -- netCDF Dataset object to be modified
            lookup_filename -- full path and filename for lookup table
            id_var_name -- name of Id variable
    
        Remarks:
            Lookup table is a CSV file with COMID, Lat, Lon, and Elev_m columns.
            Columns must be in that order and these must be the first four columns.
        """
        if self.comid_lat_lon_z_file:
            #get list of COMIDS
            lookup_table = csv_to_list(self.comid_lat_lon_z_file )
            lookup_comids = np.array([int(float(row[0])) for row in lookup_table[1:]])
        
            # Get relevant arrays while we update them
            nc_comids = self.cf_nc.variables[self.output_id_dim_name][:]
            lats = self.cf_nc.variables['lat'][:]
            lons = self.cf_nc.variables['lon'][:]
            zs = self.cf_nc.variables['z'][:]
        
            lat_min = None
            lat_max = None
            lon_min = None
            lon_max = None
            z_min = None
            z_max = None
        
            # Process each row in the lookup table
            for nc_index, nc_comid in enumerate(nc_comids):
                try:
                    lookup_index = np.where(lookup_comids == nc_comid)[0][0] + 1
                except Exception:
                    log('COMID %s misssing in comid_lat_lon_z file' % nc_comid,
                        'ERROR')
        
                lat = float(lookup_table[lookup_index][1])
                lats[nc_index] = lat
                if (lat_min) is None or lat < lat_min:
                    lat_min = lat
                if (lat_max) is None or lat > lat_max:
                    lat_max = lat
        
                lon = float(lookup_table[lookup_index][2])
                lons[nc_index] = lon
                if (lon_min) is None or lon < lon_min:
                    lon_min = lon
                if (lon_max) is None or lon > lon_max:
                    lon_max = lon
        
                z = float(lookup_table[lookup_index][3])
                zs[nc_index] = z
                if (z_min) is None or z < z_min:
                    z_min = z
                if (z_max) is None or z > z_max:
                    z_max = z
        
            # Overwrite netCDF variable values
            self.cf_nc.variables['lat'][:] = lats
            self.cf_nc.variables['lon'][:] = lons
            self.cf_nc.variables['z'][:] = zs
        
            # Update metadata
            if lat_min is not None:
                self.cf_nc.geospatial_lat_min = lat_min
            if lat_max is not None:
                self.cf_nc.geospatial_lat_max = lat_max
            if lon_min is not None:
                self.cf_nc.geospatial_lon_min = lon_min
            if lon_max is not None:
                self.cf_nc.geospatial_lon_max = lon_max
            if z_min is not None:
                self.cf_nc.geospatial_vertical_min = z_min
            if z_max is not None:
                self.cf_nc.geospatial_vertical_max = z_max
        else:
            log('No comid_lat_lon_z file. Not adding values ...', 'INFO')
            
    def _generate_time_values(self):
        """
        Generates time values for out nc file
        """
        # Populate time values
        log('writing times', 'INFO')
        total_seconds = 0
        d1970 = datetime(1970, 1, 1)
        time_array = []

        for index, time_len in enumerate(self.time_len_array):
            if index > 0:
                total_seconds += self.time_step_array[index]
                
            secs_start = int((self.start_datetime - d1970 + timedelta(seconds=total_seconds)).total_seconds())
            time_delta = self.time_step_array[index] * time_len
            total_seconds += time_delta
            if index < len(self.time_len_array)-1:
                #push back total seconds to time of next run
                total_seconds -= self.time_step_array[index]
            secs_end = secs_start + time_delta
                
            time_array.append(np.arange(secs_start, secs_end, self.time_step_array[index]))

        end_date = (self.start_datetime + timedelta(seconds=total_seconds))
        self.cf_nc.variables['time'][:] = np.concatenate(time_array)
        self.cf_nc.time_coverage_start = self.start_datetime.isoformat() + 'Z'
        self.cf_nc.time_coverage_end = end_date.isoformat() + 'Z'

    def _copy_streamflow_values(self, input_flow_var_name):
        """
        Copies streamflow values from raw output to CF file
        """
        q_var = self.cf_nc.createVariable(
            self.output_flow_var_name, 'f4', (self.output_id_dim_name, 'time'))
        q_var.long_name = 'Discharge'
        q_var.units = 'm^3/s'
        q_var.coordinates = 'time lat lon z'
        q_var.grid_mapping = 'crs'
        q_var.source = ('Generated by the Routing Application for Parallel ' +
                        'computatIon of Discharge (RAPID) river routing model.')
        q_var.references = 'http://rapid-hub.org/'
        q_var.comment = ('lat, lon, and z values taken at midpoint of river ' +
                         'reach feature')
        log('Copying streamflow values', 'INFO')
        begin_time_step_index = 1
        end_time_step_index = -1
        for raw_nc_index, raw_nc in enumerate(self.raw_nc_list):
            if raw_nc_index == 0:
                end_time_step_index = self.time_len_array[raw_nc_index]
            else:
                end_time_step_index = begin_time_step_index + self.time_len_array[raw_nc_index]
            q_var[:,begin_time_step_index:end_time_step_index] = raw_nc.variables[input_flow_var_name][:].transpose()
            begin_time_step_index = end_time_step_index
        
        #add initial flow to RAPID output file
        if self.qinit_file and self.rapid_connect_file:
            lookup_table = csv_to_list(self.rapid_connect_file)
            lookup_comids = np.array([int(float(row[0])) for row in lookup_table])
            
            init_flow_table = csv_to_list(self.qinit_file)
            
            for index, comid in enumerate(self.cf_nc.variables[self.output_id_dim_name][:]):
                try:
                    lookup_index = np.where(lookup_comids == comid)[0][0]
                except Exception:
                    log('COMID %s misssing in rapid_connect file' % comid,
                        'ERROR')
                q_var[index,0] = float(init_flow_table[lookup_index][0])
        else:
            for index, comid in enumerate(self.cf_nc.variables[self.output_id_dim_name][:]):
                q_var[index,0] = 0
                

    def convert(self):
        """
        Copies data from RAPID netCDF output to a CF-compliant netCDF file.
        """
   
        try:
            log('Processing %s ...' % self.rapid_output_file_list[0], 'INFO')
            time_start_conversion = datetime.utcnow()

            # Validate the raw netCDF file
            log('validating input netCDF file', 'INFO')
            input_id_dim_name, id_len, time_len, input_flow_var_name = self._validate_raw_nc()

            # Initialize the output file (create dimensions and variables)
            log('initializing output', 'INFO')
            self._initialize_output(time_len, id_len)

            self._generate_time_values()
            
            # Populate comid, lat, lon, z
            log('writing comid lat lon z', 'INFO')
            lookup_start = datetime.now()
            self.cf_nc.variables[self.output_id_dim_name][:] = self.raw_nc_list[0].variables[input_id_dim_name][:]
            
            self._write_comid_lat_lon_z()
            duration = str((datetime.now() - lookup_start).total_seconds())
            log('Lookup Duration (s): ' + duration, 'INFO')

            # Create a variable for streamflow. This is big, and slows down
            # previous steps if we do it earlier.
            log('Creating streamflow variable', 'INFO')
            self._copy_streamflow_values(input_flow_var_name)
            
            #close files
            for raw_nc in self.raw_nc_list:
                raw_nc.close()
            self.cf_nc.close()
            
            #delete original RAPID output
            remove_files(*self.rapid_output_file_list)

            #rename nc compliant file to original name
            os.rename(self.cf_compliant_file, self.rapid_output_file_list[0])
            log('Time to process %s' % (datetime.utcnow()-time_start_conversion), 'INFO')
        except Exception, e:
            #delete cf RAPID output
            remove_files(self.cf_compliant_file)
            log('Conversion Error %s' % e, 'ERROR')
