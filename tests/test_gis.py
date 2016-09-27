# -*- coding: utf-8 -*-
##
##  test_gis.py
##  RAPIDpy
##
##  Created by Alan D. Snow.
##  Copyright © 2016 Alan D Snow. All rights reserved.
##

from nose.tools import ok_
import os

#local import
from RAPIDpy.gis.workflow import CreateAllStaticECMWFRAPIDFiles
from RAPIDpy.gis.network import CreateNetworkConnectivityNHDPlus
from RAPIDpy.helper_functions import (compare_csv_decimal_files,
                                      remove_files)
#GLOBAL VARIABLES
MAIN_TESTS_FOLDER = os.path.dirname(os.path.abspath(__file__))
COMPARE_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'compare', 'gis')
INPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'data', 'gis')
OUTPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'output')

#------------------------------------------------------------------------------
# MAIN TEST SCRIPTS
#------------------------------------------------------------------------------
def test_gen_static_rapid_input():
    """
    Checks generating static RAPID input
    """
    print("TEST 1: TEST GENERATE STATIC RAPID INPUT DATA")
    CreateAllStaticECMWFRAPIDFiles(in_drainage_line=os.path.join(INPUT_DATA_PATH, 'flowline.shp'),
                                   river_id="COMID",
                                   length_id="LENGTHKM",
                                   slope_id="Slope",
                                   next_down_id="NextDownID",
                                   in_catchment=os.path.join(INPUT_DATA_PATH, 'catchment.shp'),
                                   catchment_river_id="FEATUREID",
                                   rapid_output_folder=OUTPUT_DATA_PATH,
                                   kfac_length_units="km",
                                   )
    
    #CHECK OUTPUT   
    #comid_lat_lon_z
    generated_comid_lat_lon_z_file = os.path.join(OUTPUT_DATA_PATH, 
                                                  "comid_lat_lon_z.csv")
    generated_comid_lat_lon_z_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                           "comid_lat_lon_z.csv")
    ok_(compare_csv_decimal_files(generated_comid_lat_lon_z_file, 
                                  generated_comid_lat_lon_z_file_solution))

    #rapid_connect
    generated_rapid_connect_file = os.path.join(OUTPUT_DATA_PATH, 
                                                "rapid_connect.csv")
    generated_rapid_connect_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                         "rapid_connect.csv")
    ok_(compare_csv_decimal_files(generated_rapid_connect_file, 
                                  generated_rapid_connect_file_solution))

    #riv_bas_id
    generated_riv_bas_id_file = os.path.join(OUTPUT_DATA_PATH, 
                                             "riv_bas_id.csv")
    generated_riv_bas_id_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                      "riv_bas_id.csv")
    ok_(compare_csv_decimal_files(generated_riv_bas_id_file, 
                                  generated_riv_bas_id_file_solution))

    #kfac
    generated_kfac_file = os.path.join(OUTPUT_DATA_PATH, 
                                       "kfac.csv")
    generated_kfac_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                "kfac.csv")
    ok_(compare_csv_decimal_files(generated_kfac_file, 
                                  generated_kfac_file_solution))
    
    #k
    generated_k_file = os.path.join(OUTPUT_DATA_PATH, 
                                    "k.csv")
    generated_k_file_solution = os.path.join(COMPARE_DATA_PATH,
                                             "k.csv")
    ok_(compare_csv_decimal_files(generated_k_file, 
                                  generated_k_file_solution))

    #x
    generated_x_file = os.path.join(OUTPUT_DATA_PATH, 
                                    "x.csv")
    generated_x_file_solution = os.path.join(COMPARE_DATA_PATH,
                                             "x.csv")
    ok_(compare_csv_decimal_files(generated_x_file, 
                                  generated_x_file_solution))

    #weight_ecmwf_t1279
    generated_weight_ecmwf_t1279_file = os.path.join(OUTPUT_DATA_PATH, 
                                                     "weight_ecmwf_t1279.csv")
    generated_weight_ecmwf_t1279_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                              "weight_ecmwf_t1279.csv")
    ok_(compare_csv_decimal_files(generated_weight_ecmwf_t1279_file, 
                                  generated_weight_ecmwf_t1279_file_solution))

    #weight_ecmwf_tco369
    generated_weight_ecmwf_tco639_file = os.path.join(OUTPUT_DATA_PATH, 
                                                      "weight_ecmwf_tco639.csv")
    generated_weight_ecmwf_tco639_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                               "weight_ecmwf_tco639.csv")
    ok_(compare_csv_decimal_files(generated_weight_ecmwf_tco639_file, 
                                  generated_weight_ecmwf_tco639_file_solution))

    #weight_era_t511
    generated_weight_era_t511_file = os.path.join(OUTPUT_DATA_PATH, 
                                                  "weight_era_t511.csv")
    generated_weight_era_t511_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                           "weight_era_t511.csv")
    ok_(compare_csv_decimal_files(generated_weight_era_t511_file, 
                                  generated_weight_era_t511_file_solution))

    remove_files(generated_comid_lat_lon_z_file,
                 generated_rapid_connect_file,
                 generated_riv_bas_id_file,
                 generated_kfac_file,
                 generated_k_file,
                 generated_x_file,
                 generated_weight_ecmwf_t1279_file,
                 generated_weight_ecmwf_tco639_file,
                 generated_weight_era_t511_file)

def test_gen_static_nhd_connect_rapid_input():
    """
    Checks generating static NHDPlus connect RAPID input
    """
    print("TEST 1: TEST GENERATE STATIC NHDPlus CONNECT RAPID INPUT DATA")
    generated_rapid_connect_file = os.path.join(OUTPUT_DATA_PATH, 
                                                "rapid_connect_nhd.csv")
    CreateNetworkConnectivityNHDPlus(in_drainage_line=os.path.join(INPUT_DATA_PATH, 'flowline.shp'),
                                     out_connectivity_file=generated_rapid_connect_file)
    #rapid_connect
    generated_rapid_connect_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                         "rapid_connect.csv")
                                                         
    ok_(compare_csv_decimal_files(generated_rapid_connect_file, 
                                  generated_rapid_connect_file_solution))

    remove_files(generated_rapid_connect_file)

if __name__ == '__main__':
    import nose
    nose.main()