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
from sys import version_info

#local imports
from .dataset import RAPIDDataset

#------------------------------------------------------------------------------
# HELPER FUNCTIONS
#------------------------------------------------------------------------------
def open_csv(csv_file, mode='r'):
    """
    Get mode depending on Python version
    Based on: http://stackoverflow.com/questions/29840849/writing-a-csv-file-in-python-that-works-for-both-python-2-7-and-python-3-3-in
    """
    if version_info[0] == 2:  # Not named on 2.6
        access = '{0}b'.format(mode)
        kwargs = {}
    else:
        access = '{0}t'.format(mode)
        kwargs = {'newline':''}
        
    return open(csv_file, access, **kwargs)
        
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
                print("{0}: {1}".format(severity, message))
        else:
                print("{0}: {1}".format(severity, message))
    else:
        raise Exception("{0}: {1}".format(severity, message))

def csv_to_list(csv_file, delimiter=','):
    """
    Reads in a CSV file and returns the contents as list,
    where every row is stored as a sublist, and each element
    in the sublist represents 1 cell in the table.
    """
    with open_csv(csv_file) as csv_con:
        if len(delimiter) > 1:
            dialect = csv.Sniffer().sniff(csv_con.read(1024), delimiters=delimiter)
            csv_con.seek(0)
            reader = csv.reader(csv_con, dialect)
        else:
            reader = csv.reader(csv_con, delimiter=delimiter)
        return list(reader)

def get_rivid_list_from_file(in_rapid_connect):
    """
    Gets the first row of rivids in rapid connect file or riv_bas_id file
    """
    rapid_connect_rivid_list = []
    with open_csv(in_rapid_connect) as csvfile:
        reader = csv.reader(csvfile, delimiter=",")
        for row in reader:
            rapid_connect_rivid_list.append(row[0])
    return np_array(rapid_connect_rivid_list, dtype=int)

def compare_csv_decimal_files(file1, file2, header=True):
    """
    This function compares two csv files
    """
    with open_csv(file1) as fh1, \
         open_csv(file2) as fh2:
        csv1 = csv.reader(fh1)
        csv2 = csv.reader(fh2)
        files_equal = True
        if header:
            files_equal = (next(csv1) == next(csv2)) #header
        while files_equal:
            try:
                try:
                    assert_almost_equal(np_array(next(csv1), dtype=np_float32),
                                        np_array(next(csv2), dtype=np_float32),
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
    with open_csv(file1) as fh1, \
         open_csv(file2) as fh2:
        csv1 = csv.reader(fh1)
        csv2 = csv.reader(fh2)
        files_equal = True
        if header:
            files_equal = (next(csv1) == next(csv2)) #header
        while files_equal:
            try:
                try:
                    row1 = next(csv1)
                    row2 = next(csv2)
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
        log("Length of COMID/rivid input not the same.",
            "ERROR")

    if not (d1.get_river_id_array() == d2.get_river_id_array()).all():
        log("COMID/rivid order is different in each dataset. Reordering data for comparison.",
            "WARNING")
        
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
                log("ALMOST EQUAL to {0} decimal places.".format(decimal_test),
                    "INFO")
                qout_same = True
                decimal_test=-1
            except AssertionError as ex:
                if decimal_test <= 1:
                    print(ex)
                decimal_test-=1
                pass
        log("Number of different timeseries: {0}".format(len(un_where_diff)),
            "INFO")
        log("COMID idexes where different: {0}".format(un_where_diff),
            "INFO")
        log("COMID idexes where different: {0}".format(un_where_diff),
            "INFO")
        index = un_where_diff[0]
        log("Dataset 1 example. COMID index: {0}".format(d1.get_qout_index(index)),
            "INFO")
        log("Dataset 2 example. COMID index: {0}".format(d2_reordered_qout[index, :]),
            "INFO")
    
    else:
        qout_same = True
        log("Output Qout data is the same.",
            "INFO")

    d1.close()
    d2.close()
    return qout_same