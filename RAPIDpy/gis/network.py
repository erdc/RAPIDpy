# -*- coding: utf-8 -*-
##
##  centroid.py
##  RAPIDpy
##
##  Created by Alan D Snow.
##  Based on RAPID_Toolbox for ArcMap
##  Copyright © 2016 Alan D Snow. All rights reserved.
##  License: BSD 3-Clause

from csv import writer as csv_writer
import numpy as np

try:
    from osgeo import ogr
except Exception:
    raise Exception("You need the gdal python package to run this tool ...")


def StreamIDNextDownIDToConnectivity(stream_id_array,
                                     next_down_id_array,
                                     out_csv_file):
    """
    Creates RAPID connect file from stream_id array and next down id array
    """
    list_all = []
    max_count_upstream = 0

    for hydroid in np.sort(stream_id_array):
        # find the HydroID of the upstreams
        list_upstreamID = stream_id_array[next_down_id_array==hydroid]
        # count the total number of the upstreams
        count_upstream = len(list_upstreamID)
        if count_upstream > max_count_upstream:
            max_count_upstream = count_upstream
        nextDownID = next_down_id_array[stream_id_array==hydroid][0]
#THIS IS REMOVED DUE TO THE FACT THAT THERE ARE STREAMS WITH ID OF ZERO
#        # replace the nextDownID with 0 if it equals to -1 (no next downstream)
#        if nextDownID == -1:
#            nextDownID = 0
        # append the list of Stream HydroID, NextDownID, Count of Upstream ID, and  HydroID of each Upstream into a larger list
        list_all.append(np.concatenate([np.array([hydroid,nextDownID,count_upstream]),list_upstreamID]).astype(int))

    with open(out_csv_file,'wb') as csvfile:
        connectwriter = csv_writer(csvfile)
        for row_list in list_all:
            out = np.concatenate([row_list, np.array([0 for i in xrange(max_count_upstream - row_list[2])])])
            connectwriter.writerow(out.astype(int))

def CreateNetworkConnectivity(in_drainage_line,
                              stream_id,
                              next_down_id,
                              out_csv_file,
                              file_geodatabase=None):
    """
    Creates Network Connectivity input CSV file for RAPID
    based on the Drainage Line feature class with HydroID and NextDownID fields
    """    
    if file_geodatabase:
        gdb_driver = ogr.GetDriverByName("OpenFileGDB")
        ogr_file_geodatabase = gdb_driver.Open(file_geodatabase, 0)
        ogr_drainage_line_shapefile_lyr = ogr_file_geodatabase.GetLayer(in_drainage_line)
    else:
        ogr_drainage_line_shapefile = ogr.Open(in_drainage_line)
        ogr_drainage_line_shapefile_lyr = ogr_drainage_line_shapefile.GetLayer()

    stream_id_array = []
    next_down_id_array = []
    for drainage_line_feature in ogr_drainage_line_shapefile_lyr:
        stream_id_array.append(drainage_line_feature.GetField(stream_id))
        next_down_id_array.append(drainage_line_feature.GetField(next_down_id))
    
    stream_id_array = np.array(stream_id_array, dtype=np.int32)
    next_down_id_array = np.array(next_down_id_array, dtype=np.int32)

    StreamIDNextDownIDToConnectivity(stream_id_array,
                                     next_down_id_array,
                                     out_csv_file)
                                     
def CreateNetworkConnectivityTauDEMTree(network_connectivity_tree_file,
                                        out_csv_file):
    """
    Creates Network Connectivity input CSV file for RAPID
    based on the TauDEM network connectivity tree file
    """    
    stream_id_array = []
    next_down_id_array = []
    with open(network_connectivity_tree_file, "rb") as csvfile:
        for row in csvfile:
            split_row = row.split()
            stream_id_array.append(split_row[0].strip()) #link number
            next_down_id_array.append(split_row[3].strip()) #next downstream link number

    stream_id_array = np.array(stream_id_array, dtype=np.int32)
    next_down_id_array = np.array(next_down_id_array, dtype=np.int32)

    StreamIDNextDownIDToConnectivity(stream_id_array,
                                     next_down_id_array,
                                     out_csv_file)

def CreateSubsetFile(in_drainage_line,
                     in_stream_id, 
                     out_csv_file,
                     file_geodatabase=None):
    """
    Creates River Basin ID subset input CSV file for RAPID
    based on the Drainage Line feature class with HydroID and NextDownID fields"
    """    
    
    if file_geodatabase:
        gdb_driver = ogr.GetDriverByName("OpenFileGDB")
        ogr_file_geodatabase = gdb_driver.Open(file_geodatabase, 0)
        ogr_drainage_line_shapefile_lyr = ogr_file_geodatabase.GetLayer(in_drainage_line)
    else:
        ogr_drainage_line_shapefile = ogr.Open(in_drainage_line)
        ogr_drainage_line_shapefile_lyr = ogr_drainage_line_shapefile.GetLayer()

    ogr_drainage_line_definition = ogr_drainage_line_shapefile_lyr.GetLayerDefn()
    
    orig_field_names = []
    for idx in xrange(ogr_drainage_line_definition.GetFieldCount()):
        orig_field_names.append(ogr_drainage_line_definition.GetFieldDefn(idx).GetName())
    
    upper_field_names = [field.upper() for field in orig_field_names]
    sort_field = None
    
    #Sort by HYDROSEQ order if the option exists
    if 'HYDROSEQ' in upper_field_names:
        #with this method, smaller is downstream
        sort_field = orig_field_names[upper_field_names.index('HYDROSEQ')]
        print("Sorting by {0}".format(sort_field))

    hydroseq_list = []
    hydroid_list = []
    '''The script line below makes sure that rows in the subset file are
       arranged in descending order of NextDownID of stream segements'''
    for drainage_line_feature in ogr_drainage_line_shapefile_lyr:
        hydroid_list.append(drainage_line_feature.GetField(in_stream_id))
        if sort_field:
            hydroseq_list.append(drainage_line_feature.GetField(sort_field))

    hydroid_list = np.array(hydroid_list, dtype=np.int32)
    if hydroseq_list:
        hydroseq_list = np.array(hydroseq_list, dtype=np.int32)
        sort_order = hydroseq_list.argsort()[::-1]
        hydroid_list = hydroid_list[sort_order]
    else:
        hydroid_list = np.sort(hydroid_list)
          

    with open(out_csv_file,'wb') as csvfile:
        connectwriter = csv_writer(csvfile)
        for hydroid in hydroid_list:
            connectwriter.writerow([hydroid])
