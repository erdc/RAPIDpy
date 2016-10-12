# -*- coding: utf-8 -*-
##
##  CreateInflowFileFromLDASRunoff.py
##  RAPIDpy
##
##  Created by Alan D. Snow (adapted from CreateInflowFileFromECMWFRunoff.py).
##  Copyright © 2015-2016 Alan D Snow. All rights reserved.
##  License: BSD-3 Clause

import netCDF4 as NET
import numpy as NUM
import os
from past.builtins import xrange

from .CreateInflowFileFromGriddedRunoff import CreateInflowFileFromGriddedRunoff

class CreateInflowFileFromLDASRunoff(CreateInflowFileFromGriddedRunoff):
    def __init__(self, lat_dim, #"g0_lat_0", 
                       lon_dim, #"g0_lon_1", 
                       lat_var, #"g0_lat_0", 
                       lon_var, #"g0_lon_1", 
                       surface_runoff_var, #"Qs_GDS0_SFC_ave1h",
                       subsurface_runoff_var, #"Qsb_GDS0_SFC_ave1h",
                       time_step_seconds, #=3*3600
                       snowmelt_runoff_var=""
                       ):
        """Define the attributes to look for"""
        self.dims_oi = [lon_dim, lat_dim]
        self.vars_oi = [lon_var, lat_var, surface_runoff_var, subsurface_runoff_var]
        
        if snowmelt_runoff_var:
            self.vars_oi += [snowmelt_runoff_var]
            
        self.length_time = {"Hourly": 1}
        self.time_step_seconds = time_step_seconds
        self.errorMessages = ["Missing Variable 'time'",
                              "Incorrect dimensions in the input runoff file.",
                              "Incorrect variables in the input runoff file.",
                              "Incorrect time variable in the input runoff file",
                              "Incorrect number of columns in the weight table",
                              "No or incorrect header in the weight table",
                              "Incorrect sequence of rows in the weight table"]
                              
        super(CreateInflowFileFromLDASRunoff, self).__init__()


    def dataValidation(self, in_nc):
        """Check the necessary dimensions and variables in the input netcdf data"""
        data_nc = NET.Dataset(in_nc)
        for dim in self.dims_oi:
            if dim not in data_nc.dimensions.keys():
                data_nc.close()
                raise Exception(self.errorMessages[1])

        for var in self.vars_oi:
            if var not in data_nc.variables.keys():
                data_nc.close()
                raise Exception(self.errorMessages[2])

        data_nc.close()
        return


    def execute(self, nc_file_list, index_list, in_weight_table, 
                out_nc, grid_type, mp_lock):
                
        """The source code of the tool."""
        if not os.path.exists(out_nc):
            print("ERROR: Outfile has not been created. You need to run: generateOutputInflowFile function ...")
            raise Exception("ERROR: Outfile has not been created. You need to run: generateOutputInflowFile function ...")
            
        if len(nc_file_list) != len(index_list):
            print("ERROR: Number of runoff files not equal to number of indices ...")
            raise Exception("ERROR: Number of runoff files not equal to number of indices ...")
        
        self.readInWeightTable(in_weight_table)

        #get indices of subset of data
        lon_ind_all = [int(i) for i in self.dict_list[self.header_wt[2]]]
        lat_ind_all = [int(j) for j in self.dict_list[self.header_wt[3]]]

        # Obtain a subset of  runoff data based on the indices in the weight table
        min_lon_ind_all = min(lon_ind_all)
        max_lon_ind_all = max(lon_ind_all)
        min_lat_ind_all = min(lat_ind_all)
        max_lat_ind_all = max(lat_ind_all)

        index_new = []
        conversion_factor = None

        #combine inflow data
        for nc_file_array_index, nc_file_array in enumerate(nc_file_list):

            index = index_list[nc_file_array_index]
            
            if not isinstance(nc_file_array, list): 
                nc_file_array = [nc_file_array]
            else:
                nc_file_array = nc_file_array
                
            data_subset_all = None

            for nc_file in nc_file_array:
                # Validate the netcdf dataset
                self.dataValidation(nc_file)

                #self.dataIdentify(nc_file, vars_oi_index)

                ''' Read the netcdf dataset'''
                data_in_nc = NET.Dataset(nc_file)

                '''Calculate water inflows'''
                #print("Calculating water inflows for {0} {1} ...".format(os.path.basename(nc_file) , grid_type))
                runoff_dimension_size = len(data_in_nc.variables[self.vars_oi[2]].dimensions)
                if runoff_dimension_size == 2:
                    #obtain subset of surface and subsurface runoff
                    data_subset_runoff = data_in_nc.variables[self.vars_oi[2]][min_lat_ind_all:max_lat_ind_all+1, min_lon_ind_all:max_lon_ind_all+1] + \
                                         data_in_nc.variables[self.vars_oi[3]][min_lat_ind_all:max_lat_ind_all+1, min_lon_ind_all:max_lon_ind_all+1]
                                         
                    #add snowmelt if applicable
                    if len(self.vars_oi) == 5:
                        #multiply by 0.1 for assumed infiltration
                        data_subset_runoff += 0.1 * data_in_nc.variables[self.vars_oi[4]][min_lat_ind_all:max_lat_ind_all+1, min_lon_ind_all:max_lon_ind_all+1]
                    
                    #get runoff dims
                    len_time_subset = 1
                    len_lat_subset = data_subset_runoff.shape[0]
                    len_lon_subset = data_subset_runoff.shape[1]

                    #reshape the runoff
                    data_subset_runoff = data_subset_runoff.reshape(len_lat_subset * len_lon_subset)

                elif runoff_dimension_size == 3:
                    #obtain subset of surface and subsurface runoff
                    data_subset_runoff = data_in_nc.variables[self.vars_oi[2]][:, min_lat_ind_all:max_lat_ind_all+1, min_lon_ind_all:max_lon_ind_all+1] + \
                                         data_in_nc.variables[self.vars_oi[3]][:, min_lat_ind_all:max_lat_ind_all+1, min_lon_ind_all:max_lon_ind_all+1]
                    #add snowmelt if applicable
                    if len(self.vars_oi) == 5:
                        #multiply by 0.1 for assumed infiltration
                        data_subset_runoff += 0.1 * data_in_nc.variables[self.vars_oi[4]][:, min_lat_ind_all:max_lat_ind_all+1, min_lon_ind_all:max_lon_ind_all+1]
                    #get runoff dims
                    len_time_subset = data_subset_runoff.shape[0]
                    len_lat_subset = data_subset_runoff.shape[1]
                    len_lon_subset = data_subset_runoff.shape[2]
                    #reshape the runoff
                    data_subset_runoff = data_subset_runoff.reshape(len_time_subset, (len_lat_subset * len_lon_subset))

                if conversion_factor == None: 
                    #get conversion_factor
                    conversion_factor = 0.001 #convert from kg/m^2 (i.e. mm) to m
                    if "s" in data_in_nc.variables[self.vars_oi[2]].getncattr("units"):
                        #that means kg/m^2/s in GLDAS v1 that is 3-hr avg, so multiply
                        #by 3 hr (ex. 3*3600). Assumed same for others (ex. 1*3600).
                        #ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/GLDAS_V1/README.GLDAS.pdf
                        #If combining files, need to take average of these, so divide by number of files
                        conversion_factor *= self.time_step_seconds/len(nc_file_array)

                data_in_nc.close()

                 
                if not index_new:
                    # compute new indices based on the data_subset_surface
                    for r in range(0,self.count):
                        ind_lat_orig = lat_ind_all[r]
                        ind_lon_orig = lon_ind_all[r]
                        index_new.append((ind_lat_orig - min_lat_ind_all)*len_lon_subset + (ind_lon_orig - min_lon_ind_all))

                #obtain a new subset of data
                if runoff_dimension_size == 2:
                    data_subset_new = data_subset_runoff[index_new]
                elif runoff_dimension_size == 3:
                    data_subset_new = data_subset_runoff[:, index_new]
                
                #FILTER DATA
                try:
                    #set masked values to zero
                    data_subset_new = data_subset_new.filled(fill_value=0)
                except AttributeError:
                    pass
                #set negative values to zero
                data_subset_new[data_subset_new<0] = 0

                #combine data
                if data_subset_all is None:
                    data_subset_all = data_subset_new
                else:
                    data_subset_all = NUM.add(data_subset_all, data_subset_new)
                    
            if runoff_dimension_size == 3 and len_time_subset > 1:
                inflow_data = NUM.zeros((len_time_subset, self.size_streamID))
            else:
                inflow_data = NUM.zeros(self.size_streamID)
            
            pointer = 0
            for stream_index in xrange(self.size_streamID):
                npoints = int(self.dict_list[self.header_wt[4]][pointer])
                # Check if all npoints points correspond to the same streamID
                if len(set(self.dict_list[self.header_wt[0]][pointer : (pointer + npoints)])) != 1:
                    print("ROW INDEX {0}".format(pointer))
                    print("COMID {0}".format(self.dict_list[self.header_wt[0]][pointer]))
                    raise Exception(self.errorMessages[2])

                area_sqm_npoints = NUM.array([float(k) for k in self.dict_list[self.header_wt[1]][pointer : (pointer + npoints)]])

                #assume data is incremental
                if runoff_dimension_size == 3:
                    data_goal = data_subset_all[:, pointer:(pointer + npoints)]
                else:
                    data_goal = data_subset_all[pointer:(pointer + npoints)]
                    
                ro_stream = data_goal * area_sqm_npoints * conversion_factor
                #filter nan
                ro_stream[NUM.isnan(ro_stream)] = 0
                
                if ro_stream.any():
                    if runoff_dimension_size == 3 and len_time_subset > 1:
                        inflow_data[:,stream_index] = ro_stream.sum(axis=1)
                    else:
                        inflow_data[stream_index] = ro_stream.sum()
                
                pointer += npoints
            #only one process is allowed to write at a time to netcdf file
            mp_lock.acquire()
            data_out_nc = NET.Dataset(out_nc, "a", format = "NETCDF3_CLASSIC")
            if runoff_dimension_size == 3 and len_time_subset > 1:
                data_out_nc.variables['m3_riv'][index*len_time_subset:(index+1)*len_time_subset,:] = inflow_data
            else:
                data_out_nc.variables['m3_riv'][index] = inflow_data
            data_out_nc.close()
            mp_lock.release()

