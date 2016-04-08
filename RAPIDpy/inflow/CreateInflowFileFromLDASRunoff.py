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

from CreateInflowFileFromGriddedRunoff import CreateInflowFileFromGriddedRunoff

class CreateInflowFileFromLDASRunoff(CreateInflowFileFromGriddedRunoff):
    def __init__(self, lat_dim="g0_lat_0", 
                       lon_dim="g0_lon_1", 
                       lat_var="g0_lat_0", 
                       lon_var="g0_lon_1", 
                       surface_runoff_var="Qs_GDS0_SFC_ave1h",
                       subsurface_runoff_var="Qsb_GDS0_SFC_ave1h",
                       time_step_seconds=3*3600):
        """Define the attributes to look for"""
        self.header_wt = ['rivid', 'area_sqm', 'lon_index', 'lat_index', 'npoints']
        self.dims_oi = [lon_dim, lat_dim]
        self.vars_oi = [lon_var, lat_var, surface_runoff_var, subsurface_runoff_var]
        self.length_time = {"Hourly": 1}
        self.time_step_seconds = time_step_seconds
        self.errorMessages = ["Missing Variable 'time'",
                              "Incorrect dimensions in the input runoff file.",
                              "Incorrect variables in the input runoff file.",
                              "Incorrect time variable in the input runoff file",
                              "Incorrect number of columns in the weight table",
                              "No or incorrect header in the weight table",
                              "Incorrect sequence of rows in the weight table"]


    def dataValidation(self, in_nc):
        """Check the necessary dimensions and variables in the input netcdf data"""
        data_nc = NET.Dataset(in_nc)
        for dim in self.dims_oi:
            if dim not in data_nc.dimensions.keys():
                data_nc.close()
                raise Exception(self.errorMessages[1])

        for var in self.vars_oi:
            if var not in data_nc.variables.keys():
                print var
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
        lon_ind_all = [long(i) for i in self.dict_list[self.header_wt[2]]]
        lat_ind_all = [long(j) for j in self.dict_list[self.header_wt[3]]]

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
                
            data_subset_surface_all = None
            data_subset_subsurface_all = None

            for nc_file in nc_file_array:
                # Validate the netcdf dataset
                self.dataValidation(nc_file)

                #self.dataIdentify(nc_file, vars_oi_index)

                ''' Read the netcdf dataset'''
                data_in_nc = NET.Dataset(nc_file)

                '''Calculate water inflows'''
                print("Calculating water inflows for {0} {1} ...".format(os.path.basename(nc_file) , grid_type))
                data_subset_surface_runoff = data_in_nc.variables[self.vars_oi[2]][min_lat_ind_all:max_lat_ind_all+1, min_lon_ind_all:max_lon_ind_all+1]
                data_subset_subsurface_runoff = data_in_nc.variables[self.vars_oi[3]][min_lat_ind_all:max_lat_ind_all+1, min_lon_ind_all:max_lon_ind_all+1]
                #check surface runoff dims
                len_lat_subset_surface = data_subset_surface_runoff.shape[0]
                len_lon_subset_surface = data_subset_surface_runoff.shape[1]
                #check subsurface runoff dims
                len_lat_subset_subsurface = data_subset_surface_runoff.shape[0]
                len_lon_subset_subsurface = data_subset_surface_runoff.shape[1]
                #make sure they are the same
                if len_lat_subset_surface != len_lat_subset_subsurface:
                    data_in_nc.close()
                    raise Exception("Surface and subsurface lat lengths do not agree ...")
                if len_lon_subset_surface != len_lon_subset_subsurface:
                    data_in_nc.close()
                    raise Exception("Surface and subsurface lon lengths do not agree ...")

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

                #reshape the runoff
                data_subset_surface_runoff = data_subset_surface_runoff.reshape(len_lat_subset_surface * len_lon_subset_surface)
                data_subset_subsurface_runoff = data_subset_subsurface_runoff.reshape(len_lat_subset_subsurface * len_lon_subset_subsurface)
                 
                if not index_new:
                    # compute new indices based on the data_subset_surface
                    for r in range(0,self.count-1):
                        ind_lat_orig = lat_ind_all[r]
                        ind_lon_orig = lon_ind_all[r]
                        index_new.append((ind_lat_orig - min_lat_ind_all)*len_lon_subset_surface + (ind_lon_orig - min_lon_ind_all))

                #obtain a new subset of data
                data_subset_surface_new = data_subset_surface_runoff[index_new]
                data_subset_subsurface_new = data_subset_subsurface_runoff[index_new]
                
                #FILTER DATA
                #set negative values to zero
                data_subset_surface_new[data_subset_surface_new<0] = 0
                data_subset_subsurface_new[data_subset_subsurface_new<0] = 0
                try:
                    #set masked values to zero
                    data_subset_surface_new = data_subset_surface_new.filled(fill_value=0)
                    data_subset_subsurface_new = data_subset_subsurface_new.filled(fill_value=0)
                except AttributeError:
                    pass

                #combine data
                if data_subset_surface_all is None:
                    data_subset_surface_all = data_subset_surface_new
                else:
                    data_subset_surface_all = NUM.add(data_subset_surface_all, data_subset_surface_new)
                    
                if data_subset_subsurface_all is None:
                    data_subset_subsurface_all = data_subset_subsurface_new
                else:
                    data_subset_subsurface_all = NUM.add(data_subset_subsurface_all, data_subset_subsurface_new)

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
                data_goal_surface = data_subset_surface_all[pointer:(pointer + npoints)]
                data_goal_subsurface = data_subset_subsurface_all[pointer:(pointer + npoints)]
                ro_stream = NUM.add(data_goal_surface, data_goal_subsurface) * area_sqm_npoints * conversion_factor
                #filter nan
                ro_stream = ro_stream[~NUM.isnan(ro_stream)]
                
                if ro_stream.any():
                    inflow_data[stream_index] = ro_stream.sum()
                
                pointer += npoints
            #only one process is allowed to write at a time to netcdf file
            mp_lock.acquire()
            data_out_nc = NET.Dataset(out_nc, "a", format = "NETCDF3_CLASSIC")
            data_out_nc.variables['m3_riv'][index] = inflow_data
            data_out_nc.close()
            mp_lock.release()

