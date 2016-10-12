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
from past.builtins import xrange

try:
    from osgeo import gdal, ogr
except Exception:
    raise Exception("You need the gdal python package to run this tool ...")

# Enable GDAL/OGR exceptions
gdal.UseExceptions()

#local
from ..helper_functions import open_csv

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

    with open_csv(out_csv_file,'w') as csvfile:
        connectwriter = csv_writer(csvfile)
        for row_list in list_all:
            out = np.concatenate([row_list, np.array([0 for i in xrange(max_count_upstream - row_list[2])])])
            connectwriter.writerow(out.astype(int))

def CreateNetworkConnectivity(in_drainage_line,
                              river_id,
                              next_down_id,
                              out_connectivity_file,
                              file_geodatabase=None):
    """
    Creates Network Connectivity input CSV file for RAPID
    based on the Drainage Line shapefile with river ID and next downstream ID fields

    Args:
        in_drainage_line(str): Path to the stream network (i.e. Drainage Line) shapefile.
        river_id(str): The name of the field with the river ID (Ex. 'HydroID', 'COMID', or 'LINKNO').
        next_down_id(str): The name of the field with the river ID of the next downstream river segment (Ex. 'NextDownID' or 'DSLINKNO').
        out_connectivity_file(str): The path to the output connectivity file.
        file_geodatabase(Optional[str]): Path to the file geodatabase. If you use this option, in_drainage_line is the name of the stream network feature class. (WARNING: Not always stable with GDAL.)
    
    Example::
    
        from RAPIDpy.gis.network import CreateNetworkConnectivity
        #------------------------------------------------------------------------------
        #main process
        #------------------------------------------------------------------------------
        if __name__ == "__main__":
            CreateNetworkConnectivity(in_drainage_line='/path/to/drainageline.shp',
                                      river_id='LINKNO',
                                      next_down_id='DSLINKNO',
                                      out_connectivity_file='/path/to/rapid_connect.csv',
                                     )
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
        stream_id_array.append(drainage_line_feature.GetField(river_id))
        next_down_id_array.append(drainage_line_feature.GetField(next_down_id))
    
    stream_id_array = np.array(stream_id_array, dtype=np.int32)
    next_down_id_array = np.array(next_down_id_array, dtype=np.int32)

    StreamIDNextDownIDToConnectivity(stream_id_array,
                                     next_down_id_array,
                                     out_connectivity_file)
                                     
def CreateNetworkConnectivityTauDEMTree(network_connectivity_tree_file,
                                        out_csv_file):
    """
    Creates Network Connectivity input CSV file for RAPID
    based on the TauDEM network connectivity tree file
    """    
    stream_id_array = []
    next_down_id_array = []
    with open_csv(network_connectivity_tree_file, "r") as csvfile:
        for row in csvfile:
            split_row = row.split()
            stream_id_array.append(split_row[0].strip()) #link number
            next_down_id_array.append(split_row[3].strip()) #next downstream link number

    stream_id_array = np.array(stream_id_array, dtype=np.int32)
    next_down_id_array = np.array(next_down_id_array, dtype=np.int32)

    StreamIDNextDownIDToConnectivity(stream_id_array,
                                     next_down_id_array,
                                     out_csv_file)
                                     
def CreateNetworkConnectivityNHDPlus(in_drainage_line,
                                     out_connectivity_file,
                                     file_geodatabase=None):
    """
    Creates Network Connectivity input CSV file for RAPID
    based on the NHDPlus drainage lines with COMID, FROMNODE, TONODE, and DIVERGENCE fields.
    
    Args:
        in_drainage_line(str): Path to the stream network (i.e. Drainage Line) shapefile.
        out_connectivity_file(str): The path to the output connectivity file.
        file_geodatabase(Optional[str]): Path to the file geodatabase. If you use this option, in_drainage_line is the name of the stream network feature class. (WARNING: Not always stable with GDAL.)
    
    Example::
    
        from RAPIDpy.gis.network import CreateNetworkConnectivityNHDPlus
        #------------------------------------------------------------------------------
        #main process
        #------------------------------------------------------------------------------
        if __name__ == "__main__":
            CreateNetworkConnectivityNHDPlus(in_drainage_line='/path/to/drainageline.shp',
                                             out_connectivity_file='/path/to/rapid_connect.csv',
                                             )
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

    def get_field_name_index(upper_field_name, upper_field_names):
        """
        returns index of field name
        """
        try:
            return upper_field_names.index(upper_field_name)
        except ValueError:
            raise IndexError("{0} not found in shapefile ..".format(upper_field_name))

    rivid_field = orig_field_names[get_field_name_index('COMID', upper_field_names)]
    fromnode_field = orig_field_names[get_field_name_index('FROMNODE', upper_field_names)]    
    tonode_field = orig_field_names[get_field_name_index('TONODE', upper_field_names)]    
    divergence_field = orig_field_names[get_field_name_index('DIVERGENCE', upper_field_names)]
        
    number_of_features = ogr_drainage_line_shapefile_lyr.GetFeatureCount()
    rivid_list = np.zeros(number_of_features, dtype=np.int32)
    fromnode_list = np.zeros(number_of_features, dtype=np.int32)
    tonode_list = np.zeros(number_of_features, dtype=np.int32)
    divergence_list = np.zeros(number_of_features, dtype=np.int32)
    for feature_idx, catchment_feature in enumerate(ogr_drainage_line_shapefile_lyr):
        rivid_list[feature_idx] = catchment_feature.GetField(rivid_field)
        fromnode_list[feature_idx] = catchment_feature.GetField(fromnode_field)
        tonode_list[feature_idx] = catchment_feature.GetField(tonode_field)
        divergence_list[feature_idx] = catchment_feature.GetField(divergence_field)

    #-------------------------------------------------------------------------------
    #Compute connectivity (based on: https://github.com/c-h-david/rrr/blob/master/src/rrr_riv_tot_gen_all_nhdplus.py)
    #-------------------------------------------------------------------------------
    fromnode_list[fromnode_list==0] = -9999
    #Some NHDPlus v1 reaches have FLOWDIR='With Digitized' but no info in VAA table
    
    fromnode_list[divergence_list==2] = -9999
    #Virtually disconnect the upstream node of all minor divergences
    divergence_list = [] #don't need this anymore

    next_down_id_list = np.zeros(number_of_features, dtype=np.int32)
    for rivid_index, rivid in enumerate(rivid_list):
        try:
            next_down_id_list[rivid_index] = rivid_list[np.where(fromnode_list==tonode_list[rivid_index])[0][0]]
        except IndexError:
            next_down_id_list[rivid_index] = -1 #this is an outlet
            pass
    #determine the downstream reach for each reach
        
    #empty unecessary lists
    fromnode_list = []
    tonode_list = []
    
    StreamIDNextDownIDToConnectivity(rivid_list,
                                     next_down_id_list,
                                     out_connectivity_file)
                                     
def CreateSubsetFile(in_drainage_line,
                     river_id, 
                     out_riv_bas_id_file,
                     file_geodatabase=None):
    """
    Creates River Basin ID subset input CSV file for RAPID
    based on the Drainage Line shapefile with river ID and next downstream ID fields

    Args:
        in_drainage_line(str): Path to the stream network (i.e. Drainage Line) shapefile.
        river_id(str): The name of the field with the river ID (Ex. 'HydroID', 'COMID', or 'LINKNO').
        out_riv_bas_id_file(str): The path to the output river basin ID subset file.
        file_geodatabase(Optional[str]): Path to the file geodatabase. If you use this option, in_drainage_line is the name of the stream network feature class. (WARNING: Not always stable with GDAL.)
    
    Example::
    
        from RAPIDpy.gis.network import CreateSubsetFile
        #------------------------------------------------------------------------------
        #main process
        #------------------------------------------------------------------------------
        if __name__ == "__main__":
            CreateSubsetFile(in_drainage_line='/path/to/drainageline.shp',
                             river_id='LINKNO',
                             out_riv_bas_id_file='/path/to/riv_bas_id.csv',
                             )
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
        hydroid_list.append(drainage_line_feature.GetField(river_id))
        if sort_field:
            hydroseq_list.append(drainage_line_feature.GetField(sort_field))

    hydroid_list = np.array(hydroid_list, dtype=np.int32)
    if hydroseq_list:
        hydroseq_list = np.array(hydroseq_list, dtype=np.int32)
        sort_order = hydroseq_list.argsort()[::-1]
        hydroid_list = hydroid_list[sort_order]
    else:
        hydroid_list = np.sort(hydroid_list)
          

    with open_csv(out_riv_bas_id_file,'w') as csvfile:
        connectwriter = csv_writer(csvfile)
        for hydroid in hydroid_list:
            connectwriter.writerow([hydroid])
