# -*- coding: utf-8 -*-
##
##  muskingum.py
##  RAPIDpy
##
##  Created by Alan D Snow.
##
##  Copyright © 2016 Alan D Snow. All rights reserved.
##  License: BSD 3-Clause

from csv import reader as csv_reader
from csv import writer as csv_writer
import numpy as np

try:
    from osgeo import ogr
except Exception:
    raise Exception("You need the gdal python package to run this tool ...")
   
#local
from ..helper_functions import csv_to_list, open_csv

#in Python 3 xrange is now range
try:
    xrange
except NameError:
    xrange = range
    pass

def CreateMuskingumKfacFile(in_drainage_line,
                            stream_id,
                            length_id,
                            slope_id,
                            celerity,
                            formula_type,
                            in_connectivity_file,
                            out_kfac_file,
                            length_units="km",
                            file_geodatabase=None):
    """
    Creates the Kfac file for calibration
    """
    if file_geodatabase:
        gdb_driver = ogr.GetDriverByName("OpenFileGDB")
        ogr_file_geodatabase = gdb_driver.Open(file_geodatabase)
        ogr_drainage_line_shapefile_lyr = ogr_file_geodatabase.GetLayer(in_drainage_line)
    else:
        ogr_drainage_line_shapefile = ogr.Open(in_drainage_line)
        ogr_drainage_line_shapefile_lyr = ogr_drainage_line_shapefile.GetLayer()

    number_of_features = ogr_drainage_line_shapefile_lyr.GetFeatureCount()
    river_id_list = np.zeros(number_of_features, dtype=np.int32)
    length_list = np.zeros(number_of_features, dtype=np.float32)
    slope_list = np.zeros(number_of_features, dtype=np.float32)
    for feature_idx, drainage_line_feature in enumerate(ogr_drainage_line_shapefile_lyr):
        river_id_list[feature_idx] = drainage_line_feature.GetField(stream_id)
        length = drainage_line_feature.GetField(length_id)
        if length is not None:
            length_list[feature_idx] = length
        slope = drainage_line_feature.GetField(slope_id)
        if slope is not None:
            slope_list[feature_idx] = slope

    if length_units == "m":
        length_list = length_list/1000.0
    elif length_units != "km":
        raise Exception("ERROR: Invalid length units supplied. Supported units are m and km.")
        
    connectivity_table = csv_to_list(in_connectivity_file)
    
    length_slope_array = []
    kfac2_array = []
    if formula_type == 1:
        print("River Length/Celerity")
    elif formula_type == 2:
        print("Eta*River Length/Sqrt(River Slope)")
    elif formula_type == 3:
        print("Eta*River Length/Sqrt(River Slope) [0.05, 0.95]")
    else:
        raise Exception("Invalid formula type. Valid range: 1-3 ...")
    
    with open_csv(out_kfac_file,'w') as kfacfile:
        kfac_writer = csv_writer(kfacfile)
        for row in connectivity_table:
            streamID = int(float(row[0]))
            
            streamIDindex = river_id_list==streamID
            # find the length
            stream_length = length_list[streamIDindex]*1000

            if formula_type >= 2:
                # find the slope
                stream_slope = slope_list[streamIDindex]
                
                if stream_slope <= 0:
                    #if no slope, take average of upstream and downstream to get it
                    nextDownID = int(float(row[1]))
                    next_down_slope = 0
                    try:
                        next_down_index = np.where(river_id_list==nextDownID)[0][0]
                        next_down_slope = slope_list[next_down_index]
                    except IndexError:
                        pass
                        
                    nextUpID = int(float(row[3]))
                    next_up_slope = 0
                    try:
                        next_up_index = np.where(river_id_list==nextUpID)[0][0]
                        next_up_slope = slope_list[next_up_index]
                    except IndexError:
                        pass
                        
                    stream_slope = (next_down_slope+next_up_slope)/2
                    if stream_slope <=0:
                        #if still no slope, set to 0.001
                        stream_slope = 0.001
                
                length_slope_array.append(stream_length/stream_slope**0.5)
                kfac2_array.append(stream_length/celerity)
            else:
                kfac = stream_length/celerity
                kfac_writer.writerow(kfac)
        
        if formula_type >= 2:
            if formula_type == 3:
                print("Filtering Data by 5th and 95th Percentiles ...")
                length_slope_array = np.array(length_slope_array)
                percentile_5 = np.percentile(length_slope_array, 5)
                percentile_95 = np.percentile(length_slope_array, 95)
                
                length_slope_array[length_slope_array<percentile_5] = percentile_5
                length_slope_array[length_slope_array>percentile_95] = percentile_95
            
            eta = np.mean(kfac2_array) / np.mean(length_slope_array)
            print("Kfac2_Avg {0}".format(np.mean(kfac2_array)))
            print("Length_Slope Avg {0}".format( np.mean(length_slope_array)))
            print("Eta {0}".format(eta))
            print("Writing Data ...")
            for len_slope in length_slope_array:
                kfac_writer.writerow(eta*len_slope)

def CreateMuskingumKFile(lambda_k,
                         in_kfac_file,
                         out_k_file):
    """
    Creates muskingum k file from kfac file
    """
    kfac_table = csv_to_list(in_kfac_file)
    
    with open_csv(out_k_file,'w') as kfile:
        k_writer = csv_writer(kfile)
        for row in kfac_table:
             k_writer.writerow([lambda_k*float(row[0])])

def CreateMuskingumXFileFromDranageLine(in_drainage_line,
                                        x_id,
                                        out_x_file,
                                        file_geodatabase=None):
    """
    Create muskingum X file from drainage line
    """
    if file_geodatabase:
        gdb_driver = ogr.GetDriverByName("OpenFileGDB")
        ogr_file_geodatabase = gdb_driver.Open(file_geodatabase)
        ogr_drainage_line_shapefile_lyr = ogr_file_geodatabase.GetLayer(in_drainage_line)
    else:
        ogr_drainage_line_shapefile = ogr.Open(in_drainage_line)
        ogr_drainage_line_shapefile_lyr = ogr_drainage_line_shapefile.GetLayer()

    with open_csv(out_x_file,'w') as kfile:
        x_writer = csv_writer(kfile)
        for drainage_line_feature in ogr_drainage_line_shapefile_lyr:
            x_writer.writerow([drainage_line_feature.GetField(x_id)])    

def CreateConstMuskingumXFile(x_val,
                              rapid_connect_file,
                              out_x_file):
    """
    Create muskingum X file from value that is constant all the way through
    """
    num_rivers = 0
    with open_csv(rapid_connect_file, "r") as csvfile:
        reader = csv_reader(csvfile)
        for row in reader:
            num_rivers+=1

    with open_csv(out_x_file,'w') as kfile:
        x_writer = csv_writer(kfile)
        for idx in xrange(num_rivers):
            x_writer.writerow([x_val])    