# -*- coding: utf-8 -*-
##
##  helper_functions.py
##  RAPIDpy
##
##  Created by Alan D Snow, 2015.
##  Copyright © 2015 Alan D Snow. All rights reserved.
##
import csv
from numpy import where, unique
from numpy.testing import assert_almost_equal
from numpy import array as np_array
from numpy import float32 as np_float32
from os import remove

#local imports
from dataset import RAPIDDataset
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
        if len(delimiter) > 1:
            dialect = csv.Sniffer().sniff(csv_con.read(1024), delimiters=delimiter)
            csv_con.seek(0)
            reader = csv.reader(csv_con, dialect)
        else:
            reader = csv.reader(csv_con, delimiter=delimiter)
        return list(reader)

def compare_csv_decimal_files(file1, file2, header=True):
    """
    This function compares two csv files
    """
    with open(file1, 'rb') as fh1, \
         open(file2, 'rb') as fh2:
        csv1 = csv.reader(fh1)
        csv2 = csv.reader(fh2)
        files_equal = True
        if header:
            files_equal = (csv1.next() == csv2.next()) #header
        while files_equal:
            try:
                try:
                    assert_almost_equal(np_array(csv1.next(), dtype=np_float32),
                                        np_array(csv2.next(), dtype=np_float32),
                                        decimal=2)
                except AssertionError:
                    files_equal = False
                    break
                    pass
            except StopIteration:
                break
                pass
    return files_equal
    
def compare_csv_timeseries_files(file1, file2, header=True):
    """
    This function compares two csv files
    """
    with open(file1, 'rb') as fh1, \
         open(file2, 'rb') as fh2:
        csv1 = csv.reader(fh1)
        csv2 = csv.reader(fh2)
        files_equal = True
        if header:
            files_equal = (csv1.next() == csv2.next()) #header
        while files_equal:
            try:
                try:
                    row1 = csv1.next()
                    row2 = csv2.next()
                    files_equal = row1[0] == row2[0] #check dates
                    assert_almost_equal(np_array(row1[1:], dtype=np_float32),
                                        np_array(row2[1:], dtype=np_float32),
                                        decimal=2)
                except AssertionError:
                    files_equal = False
                    break
                    pass
            except StopIteration:
                break
                pass
    return files_equal

def remove_files(*args):
    """
    This function removes all files input as arguments
    """
    for arg in args:
        try:
            remove(arg)
        except OSError:
            pass

def compare_qout_files(dataset1_path, dataset2_path, Qout_var="Qout"):
    """
    This function compares the output of RAPID Qout and tells you where they are different.
    """
    qout_same = False
    
    d1 = RAPIDDataset(dataset1_path)
    d2 = RAPIDDataset(dataset2_path)

    if len(d1.get_river_id_array()) != len(d2.get_river_id_array()):
        raise Exception("Length of COMID/rivid input not the same.")

    if not (d1.get_river_id_array() == d2.get_river_id_array()).all():
        print "WARNING: COMID/rivid order is different in each dataset. Reordering data for comparison."
        
        d2_reordered_reach_index_list = []
        for comid in d1.get_river_id_array():
            d2_reordered_reach_index_list.append(where(d2.get_river_id_array()==comid)[0][0])
        d2_reordered_qout = d2.get_qout_index(d2_reordered_reach_index_list)
    else:
        d2_reordered_qout = d2.get_qout()
        
    #get where the files are different
    d1_qout = d1.get_qout()
    where_diff = where(d1_qout != d2_reordered_qout)
    un_where_diff = unique(where_diff[0])
    
    #if different, check to see how different
    if un_where_diff.any():
        decimal_test = 7
        while decimal_test > 0:
            try:
                assert_almost_equal(d1_qout,
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
        print d1.get_qout_index(index)
        print "Dataset 2 example. COMID index:", index
        print d2_reordered_qout[index, :]
    
    else:
        qout_same = True
        print "Output Qout data is the same."

    d1.close()
    d2.close()
    return qout_same