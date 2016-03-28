# -*- coding: utf-8 -*-
##
##  CreateInflowFileFromHighResECMWFRunoff.py
##  spt_lsm_autorapid_process
##
##  Created by Alan D. Snow (adapted from CreateInflowFileFromECMWFRunoff.py).
##  Copyright © 2015-2016 Alan D Snow. All rights reserved.
##  License: BSD-3 Clause

import csv
import os
import netCDF4 as NET
import numpy as NUM
import re

class CreateInflowFileFromHighResECMWFRunoff(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create Inflow File From ECMWF Runoff"
        self.description = ("Creates RAPID NetCDF input of water inflow "
                            "based on ECMWF high resoulution runoff results "
                            "and previously created weight table.")
        self.canRunInBackground = False
        self.header_wt = ['StreamID', 'area_sqm', 'lon_index', 'lat_index', 'npoints']
        self.dims_oi = ['lon', 'lat', 'time']
        self.vars_oi = ["lon", "lat", "time", "RO"]
        self.length_time = {"HighRes": 125}
        self.length_time_opt = {"HighRes-1hr": 91, "HighRes-3hr": 49, "HighRes-6hr": 41}
        self.errorMessages = ["Missing Variable 'time'",
                              "Incorrect dimensions in the input ECMWF runoff file.",
                              "Incorrect variables in the input ECMWF runoff file.",
                              "Incorrect time variable in the input ECMWF runoff file",
                              "Incorrect number of columns in the weight table",
                              "No or incorrect header in the weight table",
                              "Incorrect sequence of rows in the weight table"]


    def dataValidation(self, in_nc):
        """Check the necessary dimensions and variables in the input netcdf data"""
        data_nc = NET.Dataset(in_nc)

        dims = data_nc.dimensions.keys()
        if dims != self.dims_oi:
            raise Exception(self.errorMessages[1])

        vars = data_nc.variables.keys()
        if vars != self.vars_oi:
            raise Exception(self.errorMessages[2])

        return


    def dataIdentify(self, in_nc):
        """Check if the data is Ensemble 1-51 (low resolution) or 52 (high resolution)"""
        data_nc = NET.Dataset(in_nc)
        name_time = self.vars_oi[2]
        time = data_nc.variables[name_time][:]
        diff = NUM.unique(NUM.diff(time))
        data_nc.close()
        time_interval_highres = NUM.array([1.0,3.0,6.0],dtype=float)
        if (diff == time_interval_highres).all():
            return "HighRes"
        else:
            return None

    def execute(self, in_sorted_nc_files, in_weight_table, out_nc, in_time_interval="6hr"):
        """The source code of the tool."""

        ''' Read the weight table '''
        print "Reading the weight table..."
        dict_list = {self.header_wt[0]:[], self.header_wt[1]:[], self.header_wt[2]:[],
                     self.header_wt[3]:[], self.header_wt[4]:[]}
        streamID = ""
        with open(in_weight_table, "rb") as csvfile:
            reader = csv.reader(csvfile)
            count = 0
            for row in reader:
                if count == 0:
                    #check number of columns in the weight table
                    if len(row) < len(self.header_wt):
                        raise Exception(self.errorMessages[4])
                    #check header
                    if row[1:len(self.header_wt)] != self.header_wt[1:]:
                        raise Exception(self.errorMessages[5])
                    streamID = row[0]
                    count += 1
                else:
                    for i in xrange(5):
                       dict_list[self.header_wt[i]].append(row[i])
                    count += 1

        size_streamID = len(set(dict_list[self.header_wt[0]]))
        size_time = len(in_sorted_nc_files) * 12

        # Create output inflow netcdf data
        # data_out_nc = NET.Dataset(out_nc, "w") # by default format = "NETCDF4"
        data_out_nc = NET.Dataset(out_nc, "w", format = "NETCDF3_CLASSIC")
        dim_Time = data_out_nc.createDimension('Time', size_time)
        dim_RiverID = data_out_nc.createDimension('rivid', size_streamID)
        var_m3_riv = data_out_nc.createVariable('m3_riv', 'f4', ('Time', streamID))
        data_temp = NUM.empty(shape = [size_time, size_streamID])

        lon_ind_all = [long(i) for i in dict_list[self.header_wt[2]]]
        lat_ind_all = [long(j) for j in dict_list[self.header_wt[3]]]

        # Obtain a subset of  runoff data based on the indices in the weight table
        min_lon_ind_all = min(lon_ind_all)
        max_lon_ind_all = max(lon_ind_all)
        min_lat_ind_all = min(lat_ind_all)
        max_lat_ind_all = max(lat_ind_all)

        index_pointer = 0
        for file_index, in_nc in enumerate(in_sorted_nc_files):
            # Validate the netcdf dataset
            self.dataValidation(in_nc)

            # identify if the input netcdf data is the High Resolution data with three different time intervals
            id_data = self.dataIdentify(in_nc)
            if id_data is None:
                raise Exception(self.errorMessages[3])

            ''' Read the netcdf dataset'''
            data_in_nc = NET.Dataset(in_nc)
            time = data_in_nc.variables[self.vars_oi[2]][:]

            # Check the size of time variable in the netcdf data
            if len(time) != self.length_time[id_data]:
                raise Exception(self.errorMessages[3])

            '''Calculate water inflows'''
            print "Calculating water inflows for", os.path.basename(in_nc), "..."
            data_subset_all = data_in_nc.variables[self.vars_oi[3]][:, min_lat_ind_all:max_lat_ind_all+1, min_lon_ind_all:max_lon_ind_all+1]
            len_time_subset_all = data_subset_all.shape[0]
            len_lat_subset_all = data_subset_all.shape[1]
            len_lon_subset_all = data_subset_all.shape[2]
            data_subset_all = data_subset_all.reshape(len_time_subset_all, (len_lat_subset_all * len_lon_subset_all))


            # compute new indices based on the data_subset_all
            index_new = []
            for r in range(0,count-1):
                ind_lat_orig = lat_ind_all[r]
                ind_lon_orig = lon_ind_all[r]
                index_new.append((ind_lat_orig - min_lat_ind_all)*len_lon_subset_all + (ind_lon_orig - min_lon_ind_all))

            # obtain a new subset of data
            data_subset_new = data_subset_all[:,index_new]

            # start compute inflow
            pointer = 0
            for s in range(0, size_streamID):
                npoints = int(dict_list[self.header_wt[4]][pointer])
                # Check if all npoints points correspond to the same streamID
                if len(set(dict_list[self.header_wt[0]][pointer : (pointer + npoints)])) != 1:
                    print "ROW INDEX", pointer
                    print "COMID", dict_list[self.header_wt[0]][pointer]
                    raise Exception(self.errorMessages[2])

                area_sqm_npoints = [float(k) for k in dict_list[self.header_wt[1]][pointer : (pointer + npoints)]]
                area_sqm_npoints = NUM.array(area_sqm_npoints)
                area_sqm_npoints = area_sqm_npoints.reshape(1, npoints)
                data_goal = data_subset_new[:, pointer:(pointer + npoints)]

                ''''IMPORTANT NOTE: runoff variable in ECMWF dataset is cumulative through time'''
                if "HighRes" in id_data:
                #For data with High Resolution, from Hour 0 to 90 (the first 91 time points) are of 1 hr time interval,
                # then from Hour 90 to 144 (19 time points) are of 3 hour time interval, and from Hour 144 to 240 (15 time points)
                # are of 6 hour time interval
                    # get hourly incremental time series for first 12 hours
                    ro_stream = NUM.subtract(data_goal[1:13,], data_goal[:12,]) * area_sqm_npoints


                num_data_points = len(ro_stream)
                data_temp[index_pointer:index_pointer + num_data_points,s] = ro_stream.sum(axis = 1)

                pointer += npoints

            index_pointer += num_data_points

        '''Write inflow data'''
        print "Writing inflow data..."
        var_m3_riv[:] = data_temp
        # close the input and output netcdf datasets
        data_in_nc.close()
        data_out_nc.close()

        return

if __name__ == "__main__":
    calc = CreateInflowFileFromHighResECMWFRunoff()
    in_nc_directory='/Users/Alan/Documents/RESEARCH/RAPID/ECMWF/'
    nc_files = sorted([os.path.join(in_nc_directory, filename) for filename in os.listdir(in_nc_directory) \
                  if re.search(r'.*\.52\.205\.runoff\.grib\.runoff\.netcdf', filename, re.IGNORECASE)])
    calc.execute(in_sorted_nc_files=nc_files,
                 in_weight_table='/Users/Alan/Documents/RESEARCH/RAPID/input/nfie_texas_gulf_region/rapid_updated/weight_high_res.csv',
                 out_nc='/Users/Alan/Documents/RESEARCH/RAPID/input/nfie_texas_gulf_region/rapid_updated/m3_high_res.nc')
