# -*- coding: utf-8 -*-
##
##  merge.py
##  RAPIDpy
##
##  Created by Tim Whitaker, 2015.
##  Modified by Alan D Snow, 2015-2016
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


from datetime import datetime
import os
from netCDF4 import Dataset
import numpy as np
from pytz import utc

#local
from ..dataset import RAPIDDataset
from ..helper_functions import csv_to_list, remove_files, log

#in Python 3 xrange is now range
try:
    xrange
except NameError:
    xrange = range
    pass

class ConvertRAPIDOutputToCF(object):
    """
    Class to convert RAPID output to be CF compliant. You can also use this to 
    combine consecutive RAPID output files into one file.

    Attributes:
        rapid_output_file(str or list): Path to a single RAPID Qout file or a list of RAPID Qout files.
        start_datetime(datetime): Datetime object with the time of the start of the simulation.
        time_step(int or list): Time step of simulation in seconds if single Qout file or a list of time steps corresponding to each Qout file in the *rapid_output_file*.
        qinit_file(Optional[str]): Path to the Qinit file for the simulation. If used, it will use the values in the file for the flow at simulation time zero.
        comid_lat_lon_z_file(Optional[str]): Path to comid_lat_lon_z file. If included, the spatial information will be added to the output NetCDF file.
        rapid_connect_file(Optional[str]): Path to RAPID connect file. This is required if *qinit_file* is added.
        project_name(Optional[str]): Name of your project in the output file. Default is "Default RAPID Project".
        output_id_dim_name(Optional[str]): Name of the output river ID dimension name. Default is 'rivid'.
        output_flow_var_name(Optional[str]): Name of streamflow variable in output file, typically 'Qout' or 'm3_riv'. Default is 'Qout'.
        print_debug(Optional[bool]): If True, the debug output will be printed to the console. Default is False.       

    .. warning:: This code replaces the first file with the combined output and
                 deletes the second file. BACK UP YOUR FILES!!!!
    
    Example:
    
    .. code:: python
    
        import datetime
        from RAPIDpy.postprocess import ConvertRAPIDOutputToCF
        
        file1 = "/path/to/Qout_1980to1981.nc"
        file2 = "/path/to/Qout_1981to1982.nc"
        
        cv = ConvertRAPIDOutputToCF(rapid_output_file=[file1, file2],
                                    start_datetime=datetime.datetime(2005,1,1),
                                    time_step=[3*3600, 3*3600],
                                    project_name="NLDAS(VIC)-RAPID historical flows by US Army ERDC",
                                    )
        cv.convert()
    """
    def __init__(self, rapid_output_file, 
                       start_datetime,
                       time_step,
                       qinit_file="",
                       comid_lat_lon_z_file="",
                       rapid_connect_file="",
                       project_name="Default RAPID Project",
                       output_id_dim_name='rivid',
                       output_flow_var_name='Qout',
                       print_debug=False
                       ):

       if not isinstance(rapid_output_file, list): 
           self.rapid_output_file_list = [rapid_output_file]
       else:
           self.rapid_output_file_list = rapid_output_file
       self.start_datetime = start_datetime.replace(tzinfo=utc)
       
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
        total_time_len = 1 #add one for the first flow value RAPID
                           #does not include
        id_len_list = []
        for rapid_output_file in self.rapid_output_file_list:
            qout_nc = RAPIDDataset(rapid_output_file)
            id_len_list.append(qout_nc.size_river_id)
            total_time_len += qout_nc.size_time
            self.raw_nc_list.append(qout_nc)
    
        #make sure river id lists are the same
        for id_len_undex in range(1, len(id_len_list)):
            if id_len_list[id_len_undex] != id_len_list[0]:
                raise Exception("ERROR: River ID size is different in one of the files ...")
        
        for raw_nc_index in range(1, len(self.raw_nc_list)):
            if not (self.raw_nc_list[raw_nc_index].get_river_id_array() == self.raw_nc_list[0].get_river_id_array()).all():
                raise Exception("ERROR: River IDs are different in files ...")

        return id_len_list[0], total_time_len


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
        
        #only add if user adds
        if self.comid_lat_lon_z_file and os.path.exists(self.comid_lat_lon_z_file):
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
            crs_var.epsg_code = 'EPSG:4326'  # WGS 84
            crs_var.semi_major_axis = 6378137.0
            crs_var.inverse_flattening = 298.257223563

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
        #only add if user adds
        if self.comid_lat_lon_z_file and os.path.exists(self.comid_lat_lon_z_file):
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
        d1970 = datetime(1970, 1, 1, tzinfo=utc)
        time_array = [[int((self.start_datetime - d1970).total_seconds())]]
        
        datetime_nc_start_simulation = self.start_datetime
        for raw_nc_index, raw_nc in enumerate(self.raw_nc_list):
            
            raw_nc_time = raw_nc.get_time_array(datetime_simulation_start=datetime_nc_start_simulation,
                                                simulation_time_step_seconds=self.time_step_array[raw_nc_index])
            
            time_array.append(raw_nc_time)
            datetime_nc_start_simulation = datetime.utcfromtimestamp(raw_nc_time[-1])
                
        self.cf_nc.variables['time'][:] = np.concatenate(time_array)
        end_date = datetime.utcfromtimestamp(self.cf_nc.variables['time'][-1])
        self.cf_nc.time_coverage_start = self.start_datetime.isoformat() + 'Z'
        self.cf_nc.time_coverage_end = end_date.isoformat() + 'Z'

    def _copy_streamflow_values(self):
        """
        Copies streamflow values from raw output to CF file
        """
        log('Creating streamflow variable', 'INFO')
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
        master_begin_time_step_index = 1
        master_end_time_step_index = len(self.cf_nc.dimensions['time'])
        
        #to reduce RAM, copy by chunks
        max_2d_dimension = 1000000000 #~8GB Max
        for raw_nc_index, raw_nc in enumerate(self.raw_nc_list):
            max_time_step_size = min(raw_nc.size_time, max(1, int(float(max_2d_dimension)/float(raw_nc.size_river_id))))
            raw_nc_begin_time_step_index = 0
            raw_nc_end_time_step_index = raw_nc.size_time
            for raw_nc_time_index in xrange(0, raw_nc.size_time, max_time_step_size):
                time_interval_size = max(1, min(raw_nc.size_time-raw_nc_time_index, max_time_step_size))

                raw_nc_end_time_step_index = raw_nc_begin_time_step_index + time_interval_size
                master_end_time_step_index = master_begin_time_step_index + time_interval_size
                
                q_var[:,master_begin_time_step_index:master_end_time_step_index] = raw_nc.get_qout(time_index_start=raw_nc_begin_time_step_index,
                                                                                                   time_index_end=raw_nc_end_time_step_index)
                
                master_begin_time_step_index = master_end_time_step_index
                raw_nc_begin_time_step_index = raw_nc_end_time_step_index

        log('Adding initial streamflow values', 'INFO')
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
            id_len, time_len = self._validate_raw_nc()

            # Initialize the output file (create dimensions and variables)
            log('initializing output', 'INFO')
            self._initialize_output(time_len, id_len)

            self._generate_time_values()
            
            #copy river ids over
            self.cf_nc.variables[self.output_id_dim_name][:] = self.raw_nc_list[0].get_river_id_array()

            # Populate comid, lat, lon, z
            log('writing comid lat lon z', 'INFO')
            lookup_start = datetime.now()
            self._write_comid_lat_lon_z()
            duration = str((datetime.now() - lookup_start).total_seconds())
            log('Lookup Duration (s): ' + duration, 'INFO')

            # Create a variable for streamflow. This is big, and slows down
            # previous steps if we do it earlier.
            self._copy_streamflow_values()
            
            #close files
            for raw_nc in self.raw_nc_list:
                raw_nc.close()
            self.cf_nc.close()
            
            #delete original RAPID output
            remove_files(*self.rapid_output_file_list)

            #rename nc compliant file to original name
            os.rename(self.cf_compliant_file, self.rapid_output_file_list[0])
            log('Time to process %s' % (datetime.utcnow()-time_start_conversion), 'INFO')
        except Exception:
            #delete cf RAPID output
            remove_files(self.cf_compliant_file)
            #log('Conversion Error %s' % e, 'ERROR')
            raise
