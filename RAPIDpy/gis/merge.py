# -*- coding: utf-8 -*-
##
##  merge.py
##  RAPIDpy
##
##  Created by Alan D Snow.
##
##  Copyright © 2016 Alan D Snow. All rights reserved.
##  License: BSD 3-Clause

import csv
import numpy as np

from ..helper_functions import csv_to_list
 
def MergeWeightTables(weight_table_file,
                      connectivity_file,
                      new_weight_table_file):
    """
    This function merges multiple weight tables combined into one file
    with duplicate headers removed
    """
    
    weight_table = csv_to_list(weight_table_file)
    weight_comid_list = np.array([row[0] for row in weight_table[1:]], dtype=int)
    #FEATUREID,area_sqm,lon_index,lat_index,npoints,weight,Lon,Lat
    new_weight_table = weight_table[0:1]
    replacement_row = weight_table[1][1:]
    #set area_sqm to zero
    replacement_row[0] = 0
    #set npoints to one
    replacement_row[3] = 1
 
    print("Looping ...")
    with open(connectivity_file, "rb") as fconnect:
        for row in fconnect:
            connect_rivid = int(float(row.split(",")[0]))
            try:
                #find all occurences
                comid_indicies = np.where(weight_comid_list==connect_rivid)[0]
                weight_indices = [int(d)+1 for d in comid_indicies]
                #if num occurences don't match what table says
                if len(weight_indices) > int(weight_table[weight_indices[0]][4]):
                    #print weight_table[weight_indices[0]]
                    for weight_index in weight_indices:
                        #remove if it has an area of zero
                        if float(weight_table[weight_index][1]) == 0.0:
                            #print "REMOVED:", weight_table[weight_index]
                            weight_indices.remove(weight_index)
     
                if len(weight_indices) != int(weight_table[weight_indices[0]][4]):
                    for weight_index in weight_indices:
                        print("ERROR: {0} {1}".format(weight_index, weight_table[weight_index]))
     
                for weight_index in weight_indices:
                    new_weight_table.append(weight_table[weight_index])
            except IndexError:
                print("{0} not found ...".format(connect_rivid))
                #skip if not found
                continue
 
    print("Writing ...")
    with open(new_weight_table_file, 'wb') as outfile:
        writer = csv.writer(outfile)
        writer.writerows(new_weight_table)
        
def MergeNetworkConnectFiles(old_connectivity_file,
                             new_connectivity_file):

    """
    This function merges multiple rapid_connect files combined into one file
    """
    connectivity_table = csv_to_list(old_connectivity_file)

    max_num_upstream = max([int(float(row[2])) for row in connectivity_table])
    
    print("Maximum number of upstream reaches: {0}".format(max_num_upstream))
    print("Looping ...")
    new_comid_list = np.zeros(len(connectivity_table), dtype=np.int32)
    new_connectivity_table = []
    index = 0
    for row in connectivity_table:
        try:
            comid_index = np.where(new_comid_list==int(float(row[0])))[0][0]
            if int(float(new_connectivity_table[comid_index][2]))<int(float(row[2])):
                #replace with row with more upstreams
                while len(row) < max_num_upstream + 3:
                    row.append(0)
                new_connectivity_table[comid_index] = row[:]
        except IndexError:
            #replace with row with more upstreams
            while len(row) < max_num_upstream + 3:
                row.append(0)
            new_connectivity_table.append(row)
            new_comid_list[index] = row[0]
            index += 1
            
    print("Writing ...")
    with open(new_connectivity_file, 'wb') as outfile:
        writer = csv.writer(outfile)
        writer.writerows(new_connectivity_table)

def MergeMuskingumFiles(old_connectivity_file,
                        new_connectivity_file,
                        old_muskingum_file,
                        new_muskingum_file):
    """
    This function reorganizes your combined and in order muskingum k,kfac,x files
    from the old rapid_connect file to the updated version
    """

    old_connectivity_table = csv_to_list(old_connectivity_file)
    old_comid_list = np.array([row[0] for row in old_connectivity_table], dtype=int)
    old_connectivity_table = None
    new_connectivity_table = csv_to_list(new_connectivity_file)
    old_muskingum_table = csv_to_list(old_muskingum_file)
    new_table = []

    print("Looping ...")
    for row in new_connectivity_table:
        try:
            comid_index = np.where(old_comid_list==int(float(row[0])))[0][0]
            new_table.append(old_muskingum_table[comid_index])
        except IndexError:
            #print "SKIPPED:", row[0]
            continue

    print("Writing ...")
    with open(new_muskingum_file, 'wb') as outfile:
        writer = csv.writer(outfile)
        writer.writerows(new_table)    