# -*- coding: utf-8 -*-
##
##  helper_functions.py
##  RAPIDpy
##
##  Created by Alan D Snow, 2015.
##  Copyright © 2015 Alan D Snow. All rights reserved.
##
from csv import reader as csvreader
from csv import writer as csvwriter
from netCDF4 import Dataset
from numpy import where, unique
from numpy.ma import masked
from numpy.testing import assert_almost_equal
from os import remove
import time

#------------------------------------------------------------------------------
# HELPER FUNCTIONS
#------------------------------------------------------------------------------
def csv_to_list(csv_file, delimiter=','):
    """
    Reads in a CSV file and returns the contents as list,
    where every row is stored as a sublist, and each element
    in the sublist represents 1 cell in the table.

    """
    with open(csv_file, 'rb') as csv_con:
        reader = csvreader(csv_con, delimiter=delimiter)
        return list(reader)

    
def remove_files(*args):
    """
    This function removes all files input as arguments
    """
    for arg in args:
        try:
            remove(arg)
        except OSError:
            pass

def get_rapid_timeseries(nc_file_handle, reach_index, id_dim_name, out_var):
    """
    This function get's a time series from RAPID output for a specific
    reach id
    """
    qout_dims = nc_file_handle.variables[out_var].dimensions
    if qout_dims[0] == id_dim_name:
        return nc_file_handle.variables[out_var][reach_index, :]
    else:
        return nc_file_handle.variables[out_var][:, reach_index]

def compare_qout_files(dataset1_path, dataset2_path, Qout_var="Qout"):
    """
    This function compares the output of RAPID Qout and tells you where they are different.
    """
    qout_same = False
    
    d1 = Dataset(dataset1_path)
    dims1 = d1.dimensions
    id_dim_name1 = 'COMID'
    if 'rivid' in dims1:
        id_dim_name1 = 'rivid'
    d2 = Dataset(dataset2_path)
    dims2 = d2.dimensions
    id_dim_name2 = 'COMID'
    if 'rivid' in dims2:
        id_dim_name2 = 'rivid'

    if len(d1.variables[id_dim_name1][:]) != len(d2.variables[id_dim_name2][:]):
        raise Exception("Length of COMID/rivid input not the same.")

    if not (d1.variables[id_dim_name1][:] == d2.variables[id_dim_name2][:]).all():
        print "WARNING: COMID/rivid order is different in each dataset. Reordering data for comparison."
        
        d2_comid_list = d2.variables[id_dim_name2][:]
        d2_reordered_comid_list = []
        for comid in d1.variables[id_dim_name1][:]:
            d2_reordered_comid_list.append(where(d2_comid_list==comid)[0][0])
        qout_dimensions = d2.variables[Qout_var].dimensions
        if qout_dimensions[0].lower() == 'time' and \
           qout_dimensions[1].lower() == id_dim_name2.lower():
            d2_reordered_qout = d2.variables[Qout_var][:,d2_reordered_comid_list]
        elif qout_dimensions[1].lower() == 'time' and \
             qout_dimensions[0].lower() == id_dim_name2.lower():
            d2_reordered_qout = d2.variables[Qout_var][d2_reordered_comid_list,:]
        else:
            raise Exception("Invalid RAPID Qout file.")
    else:
        d2_reordered_qout = d2.variables[Qout_var][:]
        
    #get where the files are different
    where_diff = where(d1.variables[Qout_var][:] != d2_reordered_qout)
    un_where_diff = unique(where_diff[0])
    
    #if different, check to see how different
    if un_where_diff.any():
        decimal_test = 7
        while decimal_test > 0:
            try:
                assert_almost_equal(d1.variables[Qout_var][:], 
                                    d2_reordered_qout, 
                                    decimal=decimal_test)
                print "\nALMOST EQUAL to", decimal_test, "decimal places.\n"
                qout_same = True
                decimal_test=-1
            except AssertionError as ex:
                if decimal_test <= 1:
                    print ex
                decimal_test-=1
                pass
        print "Number of different timeseries:", len(un_where_diff)
        print "COMID idexes where different:"
        print un_where_diff
        index = un_where_diff[0]
        print "Dataset 1 example. COMID index:", index
        print d1.variables[Qout_var][index, :]
        print "Dataset 2 example. COMID index:", index
        print d2_reordered_qout[index, :]
    
    else:
        qout_same = True
        print "Output Qout data is the same."

    d1.close()
    d2.close()
    return qout_same

def write_flows_to_csv(path_to_rapid_qout_file, path_to_output_file, 
                       reach_index=None, reach_id=None, daily=False, 
                       out_var='Qout'):
    """
        Write out RAPID output to CSV file
    """
    data_nc = Dataset(path_to_rapid_qout_file)
    
    dims = data_nc.dimensions
    id_dim_name = 'COMID'
    if 'rivid' in dims:
        id_dim_name = 'rivid'

    if reach_id != None:
        reach_ids = data_nc.variables[id_dim_name][:]
        reach_index = where(reach_ids==reach_id)[0][0]

    nc_vars = data_nc.variables.keys()
    
    time_var_valid = False
    if 'time' in nc_vars:
        if len(data_nc.dimensions['time'])>0:
            if not (data_nc.variables['time'][:] == masked).any():
                try:
                    time.gmtime(data_nc.variables['time'][0])
                    time_var_valid = True
                except ValueError:
                    pass

    #analyze and write
    if time_var_valid:
        if daily:
            current_day = time.gmtime(data_nc.variables['time'][0])
            flow = 0
            num_days = 0
            
            qout_arr = get_rapid_timeseries(data_nc, reach_index, id_dim_name, out_var)
            with open(path_to_output_file, 'w') as outcsv:
                writer = csvwriter(outcsv)
                for idx, t in enumerate(data_nc.variables['time'][:]):
                    var_time = time.gmtime(t)
                    if current_day.tm_yday == var_time.tm_yday:
                        flow += qout_arr[idx]
                        num_days += 1
                    else:
                        if num_days > 0:
                            #write last average
                            writer.writerow([time.strftime("%Y/%m/%d", current_day), flow/num_days])
                        
                        #start new average
                        current_day = var_time
                        num_days = 1
                        flow = qout_arr[idx]
        else:
            qout = get_rapid_timeseries(data_nc, reach_index, id_dim_name, out_var)
            time_array = data_nc.variables['time'][:]
            with open(path_to_output_file, 'w') as outcsv:
                writer = csvwriter(outcsv)
                for index in xrange(len(qout)):
                    var_time = time.gmtime(time_array[index])
                    writer.writerow([time.strftime("%Y/%m/%d %H:00", var_time), qout[index]])

    else:
        print "Valid time variable not found. Printing values only ..."
        qout = get_rapid_timeseries(data_nc, reach_index, id_dim_name, out_var)
        with open(path_to_output_file, 'w') as outcsv:
            writer = csvwriter(outcsv)
            for index in xrange(len(qout)):
                writer.writerow([index, qout[index]])

    data_nc.close()
    return time_var_valid