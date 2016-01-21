# -*- coding: utf-8 -*-
##
##  helper_functions.py
##  RAPIDpy
##
##  Created by Alan D Snow, 2015.
##  Copyright © 2015 Alan D Snow. All rights reserved.
##
from csv import reader as csvreader
from netCDF4 import Dataset
import numpy as np
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
            d2_reordered_comid_list.append(np.where(d2_comid_list==comid)[0][0])
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
    where_diff = np.where(d1.variables[Qout_var][:] != d2_reordered_qout)
    un_where_diff = np.unique(where_diff[0])
    
    #if different, check to see how different
    if un_where_diff.any():
        decimal_test = 7
        while decimal_test > 0:
            try:
                np.testing.assert_almost_equal(d1.variables[Qout_var][:], 
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

def write_flows_to_csv(path_to_file, ind=None, reach_id=None, daily=False, out_var='Qout'):
    """
        Write out RAPID output to CSV file
    """
    data_nc = Dataset(path_to_file)
    
    if reach_id != None:
        dims = data_nc.dimensions
        id_dim_name = 'COMID'
        if 'rivid' in dims:
            id_dim_name = 'rivid'
        reach_ids = data_nc.variables[id_dim_name][:]
        ind = np.where(reach_ids==reach_id)[0][0]
        print ind, reach_id
    nc_vars = data_nc.variables.keys()

    #analyze and write
    if 'time' in nc_vars:
        if daily:
            current_day = datetime.datetime.fromtimestamp(data_nc.variables['time'][0], tz=utc)
            flow = 0
            num_days = 0
            qout_arr = data_nc.variables[out_var][ind, :]
            
            with open(os.path.join(os.path.dirname(path_to_file), "daily_flows.csv"), 'w') as outcsv:
                writer = csvwriter(outcsv)
                for idx, t in enumerate(data_nc.variables['time'][:]):
                    var_time = datetime.datetime.fromtimestamp(t, tz=utc)
                    if current_day.day == var_time.day:
                        flow += qout_arr[idx]
                        num_days += 1
                    else:
                        if num_days > 0:
                            #write last average
                            writer.writerow([current_day.strftime("%Y/%m/%d"), flow/num_days])
                        
                        #start new average
                        current_day = var_time
                        num_days = 1
                        flow = qout_arr[idx]
        else:
            qout = data_nc.variables[out_var][ind, :]
            time = data_nc.variables['time'][:]
            with open(os.path.join(os.path.dirname(path_to_file), "flows.csv"), 'w') as outcsv:
                writer = csvwriter(outcsv)
                for index in xrange(len(qout)):
                    var_time = datetime.datetime.fromtimestamp(time[index], tz=utc)
                    writer.writerow([var_time.strftime("%Y/%m/%d %H:00"), qout[index]])

    else:
        qout = data_nc.variables[out_var][:, ind]
        with open(os.path.join(os.path.dirname(path_to_file), "flows.csv"), 'w') as outcsv:
            writer = csvwriter(outcsv)
            for index in xrange(len(qout)):
                writer.writerow([index, qout[index]])

    data_nc.close()
