# -*- coding: utf-8 -*-
##
##  test_gis.py
##  RAPIDpy
##
##  Created by Alan D. Snow.
##  Copyright © 2016 Alan D Snow. All rights reserved.
##

from glob import glob
from numpy.testing import assert_almost_equal
from numpy import array
import os
from osgeo import ogr
import pytest

#local import
from RAPIDpy.gis.weight import CreateWeightTableECMWF, CreateWeightTableLDAS
from RAPIDpy.gis.workflow import CreateAllStaticECMWFRAPIDFiles
from RAPIDpy.gis.network import CreateNetworkConnectivityNHDPlus
from RAPIDpy.gis.muskingum import CreateMuskingumKfacFile, CreateMuskingumXFileFromDranageLine
from RAPIDpy.gis.taudem import TauDEM
from RAPIDpy.helper_functions import (compare_csv_decimal_files,
                                      remove_files)
#GLOBAL VARIABLES
MAIN_TESTS_FOLDER = os.path.dirname(os.path.abspath(__file__))
COMPARE_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'compare', 'gis')
GIS_INPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'data', 'gis')
RAPID_INPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'input')
OUTPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'output')
LSM_INPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'data','lsm_grids')
TAUDEM_EXE_PATH = os.path.join(MAIN_TESTS_FOLDER,
                               "..", "..", "TauDEM", "src")


#------------------------------------------------------------------------------
# MAIN TEST SCRIPTS
#------------------------------------------------------------------------------
def test_gen_static_rapid_input():
    """
    Checks generating static RAPID input
    """
    print("TEST 1: TEST GENERATE STATIC RAPID INPUT DATA")
    CreateAllStaticECMWFRAPIDFiles(in_drainage_line=os.path.join(GIS_INPUT_DATA_PATH, 'flowline.shp'),
                                   river_id="COMID",
                                   length_id="LENGTHKM",
                                   slope_id="Slope",
                                   next_down_id="NextDownID",
                                   in_catchment=os.path.join(GIS_INPUT_DATA_PATH, 'catchment.shp'),
                                   catchment_river_id="FEATUREID",
                                   rapid_output_folder=OUTPUT_DATA_PATH,
                                   kfac_length_units="km",
                                   )

    #CHECK OUTPUT
    #comid_lat_lon_z
    generated_comid_lat_lon_z_file = os.path.join(OUTPUT_DATA_PATH,
                                                  "comid_lat_lon_z.csv")
    generated_comid_lat_lon_z_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                           "comid_lat_lon_z.csv")
    assert (compare_csv_decimal_files(generated_comid_lat_lon_z_file,
                                  generated_comid_lat_lon_z_file_solution))

    #rapid_connect
    generated_rapid_connect_file = os.path.join(OUTPUT_DATA_PATH,
                                                "rapid_connect.csv")
    generated_rapid_connect_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                         "rapid_connect.csv")
    assert (compare_csv_decimal_files(generated_rapid_connect_file,
                                  generated_rapid_connect_file_solution))

    #riv_bas_id
    generated_riv_bas_id_file = os.path.join(OUTPUT_DATA_PATH,
                                             "riv_bas_id.csv")
    generated_riv_bas_id_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                      "riv_bas_id.csv")
    assert (compare_csv_decimal_files(generated_riv_bas_id_file,
                                  generated_riv_bas_id_file_solution))

    #kfac
    generated_kfac_file = os.path.join(OUTPUT_DATA_PATH,
                                       "kfac.csv")
    generated_kfac_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                "kfac.csv")
    assert (compare_csv_decimal_files(generated_kfac_file,
                                  generated_kfac_file_solution))

    #k
    generated_k_file = os.path.join(OUTPUT_DATA_PATH,
                                    "k.csv")
    generated_k_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                             "k.csv")
    assert (compare_csv_decimal_files(generated_k_file,
                                  generated_k_file_solution))

    #x
    generated_x_file = os.path.join(OUTPUT_DATA_PATH,
                                    "x.csv")
    generated_x_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                             "x.csv")
    assert (compare_csv_decimal_files(generated_x_file,
                                  generated_x_file_solution))

    #weight_ecmwf_t1279
    generated_weight_ecmwf_t1279_file = os.path.join(OUTPUT_DATA_PATH,
                                                     "weight_ecmwf_t1279.csv")
    generated_weight_ecmwf_t1279_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                              "weight_ecmwf_t1279.csv")
    assert (compare_csv_decimal_files(generated_weight_ecmwf_t1279_file,
                                  generated_weight_ecmwf_t1279_file_solution))

    #weight_ecmwf_tco369
    generated_weight_ecmwf_tco639_file = os.path.join(OUTPUT_DATA_PATH,
                                                      "weight_ecmwf_tco639.csv")
    generated_weight_ecmwf_tco639_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                               "weight_ecmwf_tco639.csv")
    assert (compare_csv_decimal_files(generated_weight_ecmwf_tco639_file,
                                  generated_weight_ecmwf_tco639_file_solution))

    #weight_era_t511
    generated_weight_era_t511_file = os.path.join(OUTPUT_DATA_PATH,
                                                  "weight_era_t511.csv")
    generated_weight_era_t511_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                           "weight_era_t511.csv")
    assert (compare_csv_decimal_files(generated_weight_era_t511_file,
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
    print("TEST 2: TEST GENERATE STATIC NHDPlus CONNECT RAPID INPUT DATA")
    generated_rapid_connect_file = os.path.join(OUTPUT_DATA_PATH,
                                                "rapid_connect_nhd.csv")
    CreateNetworkConnectivityNHDPlus(in_drainage_line=os.path.join(GIS_INPUT_DATA_PATH, 'flowline.shp'),
                                     out_connectivity_file=generated_rapid_connect_file)
    #rapid_connect
    generated_rapid_connect_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                         "rapid_connect.csv")

    assert (compare_csv_decimal_files(generated_rapid_connect_file,
                                  generated_rapid_connect_file_solution))

    remove_files(generated_rapid_connect_file)

def test_gen_weight_table_era20cm():
    """
    Checks generating weight table for ERA 20CM grid
    """
    print("TEST 3: TEST GENERATE WEIGHT TABLE FOR ERA 20CM GRIDS")
    generated_weight_table_file = os.path.join(OUTPUT_DATA_PATH,
                                               "weight_era_t159.csv")
    #rapid_connect
    rapid_connect_file = os.path.join(COMPARE_DATA_PATH, "x-x",
                                      "rapid_connect.csv")

    lsm_grid = os.path.join(LSM_INPUT_DATA_PATH, "era20cm", "era_20cm_runoff_20000129_0.nc")
    CreateWeightTableECMWF(in_ecmwf_nc=lsm_grid,
                           in_catchment_shapefile=os.path.join(GIS_INPUT_DATA_PATH, 'catchment.shp'),
                           river_id="FEATUREID",
                           in_connectivity_file=rapid_connect_file,
                           out_weight_table=generated_weight_table_file)

    generated_weight_table_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                        "weight_era_t159.csv")
    assert (compare_csv_decimal_files(generated_weight_table_file,
                                  generated_weight_table_file_solution))

    remove_files(generated_weight_table_file)

def test_gen_weight_table_era_t255():
    """
    Checks generating weight table for ERA T255 grid
    """
    print("TEST 4: TEST GENERATE WEIGHT TABLE FOR ERA T255 GRIDS")
    generated_weight_table_file = os.path.join(OUTPUT_DATA_PATH,
                                               "weight_era_t255.csv")
    #rapid_connect
    rapid_connect_file = os.path.join(COMPARE_DATA_PATH, "x-x",
                                      "rapid_connect.csv")

    lsm_grid = os.path.join(LSM_INPUT_DATA_PATH, "erai3t255", "era_interim_runoff_20140820.nc")
    CreateWeightTableECMWF(in_ecmwf_nc=lsm_grid,
                           in_catchment_shapefile=os.path.join(GIS_INPUT_DATA_PATH, 'catchment.shp'),
                           river_id="FEATUREID",
                           in_connectivity_file=rapid_connect_file,
                           out_weight_table=generated_weight_table_file)

    generated_weight_table_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                        "weight_era_t255.csv")
    assert (compare_csv_decimal_files(generated_weight_table_file,
                                  generated_weight_table_file_solution))

    remove_files(generated_weight_table_file)

def test_gen_weight_table_era_t511_24hr():
    """
    Checks generating weight table for ERA T511 24hr grid
    """
    print("TEST 5: TEST GENERATE WEIGHT TABLE FOR ERA T511 24hr GRIDS")
    generated_weight_table_file = os.path.join(OUTPUT_DATA_PATH,
                                               "weight_era_t511.csv")
    #rapid_connect
    rapid_connect_file = os.path.join(COMPARE_DATA_PATH,"x-x",
                                      "rapid_connect.csv")

    lsm_grid = os.path.join(LSM_INPUT_DATA_PATH, "erai24", "19990109_erai_runoff.grib.nc")
    CreateWeightTableECMWF(in_ecmwf_nc=lsm_grid,
                           in_catchment_shapefile=os.path.join(GIS_INPUT_DATA_PATH, 'catchment.shp'),
                           river_id="FEATUREID",
                           in_connectivity_file=rapid_connect_file,
                           out_weight_table=generated_weight_table_file)

    generated_weight_table_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                        "weight_era_t511.csv")
    assert (compare_csv_decimal_files(generated_weight_table_file,
                                  generated_weight_table_file_solution))

    remove_files(generated_weight_table_file)

def test_gen_weight_table_gldas2():
    """
    Checks generating weight table for GLDAS V2 grid
    """
    print("TEST 6: TEST GENERATE WEIGHT TABLE FOR GLDAS V2 GRIDS")
    generated_weight_table_file = os.path.join(OUTPUT_DATA_PATH,
                                               "weight_gldas2.csv")
    #rapid_connect
    rapid_connect_file = os.path.join(COMPARE_DATA_PATH, "x-x",
                                      "rapid_connect.csv")

    lsm_grid = os.path.join(LSM_INPUT_DATA_PATH, "gldas2", "GLDAS_NOAH025_3H.A20101231.0000.020.nc4")
    CreateWeightTableLDAS(in_ldas_nc=lsm_grid,
                          in_nc_lon_var="lon",
                          in_nc_lat_var="lat",
                          in_catchment_shapefile=os.path.join(GIS_INPUT_DATA_PATH, 'catchment.shp'),
                          river_id="FEATUREID",
                          in_connectivity_file=rapid_connect_file,
                          out_weight_table=generated_weight_table_file)

    generated_weight_table_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                        "weight_gldas2.csv")
    assert (compare_csv_decimal_files(generated_weight_table_file,
                                      generated_weight_table_file_solution))

    remove_files(generated_weight_table_file)

def test_gen_weight_table_lis():
    """
    Checks generating weight table for LIS grid
    """
    print("TEST 7: TEST GENERATE WEIGHT TABLE FOR LIS GRIDS")
    generated_weight_table_file = os.path.join(OUTPUT_DATA_PATH,
                                               "weight_lis.csv")
    #rapid_connect
    rapid_connect_file = os.path.join(COMPARE_DATA_PATH, "u-k",
                                      "rapid_connect.csv")

    lsm_grid = os.path.join(LSM_INPUT_DATA_PATH, "lis", "LIS_HIST_201101210000.d01.nc")
    CreateWeightTableLDAS(in_ldas_nc=lsm_grid,
                          in_nc_lon_var="lon",
                          in_nc_lat_var="lat",
                          in_catchment_shapefile=os.path.join(GIS_INPUT_DATA_PATH, 'u-k', 'CatchmentSubset.shp'),
                          river_id="DrainLnID",
                          in_connectivity_file=rapid_connect_file,
                          out_weight_table=generated_weight_table_file)

    generated_weight_table_file_solution = os.path.join(COMPARE_DATA_PATH, "u-k",
                                                        "weight_lis.csv")
    assert (compare_csv_decimal_files(generated_weight_table_file,
                                      generated_weight_table_file_solution))

    remove_files(generated_weight_table_file)

def test_gen_weight_table_lis_no_intersect():
    """
    Checks generating weight table for LIS grid with no intersect
    """
    print("TEST 8: TEST GENERATE WEIGHT TABLE FOR LIS GRIDS WITH NO INTERSECT")
    generated_weight_table_file = os.path.join(OUTPUT_DATA_PATH,
                                               "weight_lis_no_intersect.csv")
    #rapid_connect
    rapid_connect_file = os.path.join(GIS_INPUT_DATA_PATH, "uk-no_intersect",
                                      "rapid_connect_45390.csv")

    lsm_grid = os.path.join(LSM_INPUT_DATA_PATH, "lis", "LIS_HIST_201101210000.d01.nc")
    CreateWeightTableLDAS(in_ldas_nc=lsm_grid,
                          in_nc_lon_var="lon",
                          in_nc_lat_var="lat",
                          in_catchment_shapefile=os.path.join(GIS_INPUT_DATA_PATH, 'uk-no_intersect', 'Catchment_thames_drainID45390.shp'),
                          river_id="DrainLnID",
                          in_connectivity_file=rapid_connect_file,
                          out_weight_table=generated_weight_table_file)

    generated_weight_table_file_solution = os.path.join(COMPARE_DATA_PATH, "uk-no_intersect",
                                                        "weight_lis_no_intersect.csv")
    assert (compare_csv_decimal_files(generated_weight_table_file,
                                      generated_weight_table_file_solution))

    remove_files(generated_weight_table_file)

def test_gen_weight_table_joules():
    """
    Checks generating weight table for Joules grid
    """
    print("TEST 9: TEST GENERATE WEIGHT TABLE FOR Joules GRIDS")
    generated_weight_table_file = os.path.join(OUTPUT_DATA_PATH,
                                               "weight_joules.csv")
    #rapid_connect
    rapid_connect_file = os.path.join(COMPARE_DATA_PATH, "u-k",
                                      "rapid_connect.csv")

    lsm_grid = os.path.join(LSM_INPUT_DATA_PATH, "joules", "ukv_test.runoff.20080803_00.nc")
    CreateWeightTableLDAS(in_ldas_nc=lsm_grid,
                          in_nc_lon_var="east_west",
                          in_nc_lat_var="north_south",
                          in_catchment_shapefile=os.path.join(GIS_INPUT_DATA_PATH, 'u-k', 'CatchmentSubset.shp'),
                          river_id="DrainLnID",
                          in_connectivity_file=rapid_connect_file,
                          out_weight_table=generated_weight_table_file)

    generated_weight_table_file_solution = os.path.join(COMPARE_DATA_PATH, "u-k",
                                                        "weight_joules.csv")
    assert (compare_csv_decimal_files(generated_weight_table_file,
                                      generated_weight_table_file_solution))

    remove_files(generated_weight_table_file)

#==============================================================================
# def test_gen_weight_table_wrf():
#     """
#     Checks generating weight table for WRF grid
#     """
#     print("TEST 9: TEST GENERATE WEIGHT TABLE FOR WRF GRIDS")
#     generated_weight_table_file = os.path.join(OUTPUT_DATA_PATH,
#                                                "weight_wrf.csv")
#     #rapid_connect
#     rapid_connect_file = os.path.join(COMPARE_DATA_PATH, "m-s",
#                                       "rapid_connect.csv")
#
# #==============================================================================
# #     drainage_line = os.path.join(GIS_INPUT_DATA_PATH, 'm-s', 'flowline_subset.shp')
# #     CreateAllStaticRAPIDFiles(in_drainage_line=drainage_line,
# #                               river_id="COMID",
# #                               length_id="LENGTHKM",
# #                               slope_id="SLOPE",
# #                               next_down_id="NextDownID",
# #                               rapid_output_folder=os.path.join(GIS_INPUT_DATA_PATH, 'm-s')
# #                               )
# #==============================================================================
#
#
#     lsm_grid = os.path.join(LSM_INPUT_DATA_PATH, "wrf", "diffro_d02_20080601010000.nc")
#     CreateWeightTableLDAS(in_ldas_nc=lsm_grid,
#                           in_nc_lon_var="XLONG",
#                           in_nc_lat_var="XLAT",
#                           in_catchment_shapefile=os.path.join(GIS_INPUT_DATA_PATH, 'm-s', 'catchment_subset.shp'),
#                           river_id="FEATUREID",
#                           in_connectivity_file=rapid_connect_file,
#                           out_weight_table=generated_weight_table_file)
#
#     generated_weight_table_file_solution = os.path.join(COMPARE_DATA_PATH, "m-s",
#                                                         "weight_wrf.csv")
#     assert (compare_csv_decimal_files(generated_weight_table_file,
#                                   generated_weight_table_file_solution))
#
#     remove_files(generated_weight_table_file)
#==============================================================================


def test_extract_sub_network_taudem():
    """
    Checks extracting sub network from larger network
    """
    print("TEST 10: TEST EXTRACTING SUB NETWORK FROM LARGER NETWORK")
    td = TauDEM()

    subset_network_file = os.path.join(OUTPUT_DATA_PATH, "DrainageLineSubset2.shp")
    #to extract a specific network
    td.extractSubNetwork(network_file=os.path.join(GIS_INPUT_DATA_PATH, 'u-k', "DrainageLineSubset.shp"),
                         out_subset_network_file=subset_network_file,
                         outlet_ids=[42911], #list of outlet ids
                         river_id_field="HydroID",
                         next_down_id_field="NextDownID",
                         river_magnitude_field="HydroID",
                         safe_mode=False,
                         )

    #to extract the subset watersheds using subset river network
    subset_watershed_file = os.path.join(OUTPUT_DATA_PATH,"CatchmentSubset2.shp")
    td.extractSubsetFromWatershed(subset_network_file=subset_network_file,
                                  subset_network_river_id_field="HydroID",
                                  watershed_file=os.path.join(GIS_INPUT_DATA_PATH, 'u-k', 'CatchmentSubset.shp'),
                                  watershed_network_river_id_field="DrainLnID",
                                  out_watershed_subset_file=subset_watershed_file)


    #Test results
    subset_network_shapefile = ogr.Open(subset_network_file)
    subset_network_layer = subset_network_shapefile.GetLayer()

    ogr_watershed_shapefile = ogr.Open(subset_watershed_file)
    ogr_watershed_shapefile_lyr = ogr_watershed_shapefile.GetLayer()

    number_of_network_features = subset_network_layer.GetFeatureCount()
    number_of_watershed_features = ogr_watershed_shapefile_lyr.GetFeatureCount()

    #count number of features
    assert (number_of_network_features==7)
    assert (number_of_watershed_features==7)

    #make sure IDs correct
    network_id_list = [42911,42891,42747,42748,42892,42841,42846]
    for feature_idx, network_feature in enumerate(subset_network_layer):
        assert (network_feature.GetField("HydroID") in network_id_list)
    for feature_idx, watershed_feature in enumerate(ogr_watershed_shapefile_lyr):
        assert (watershed_feature.GetField("DrainLnID") in network_id_list)

    #make sure all fields are there

     #TEST WATERSHED
    subset_watershed_layer_defn = ogr_watershed_shapefile_lyr.GetLayerDefn()
    num_watershed_fields = subset_watershed_layer_defn.GetFieldCount()

    watershed_field_names = ['Shape_Leng','Shape_Area','HydroID','GridID','DrainLnID']
    assert (num_watershed_fields==len(watershed_field_names))
    for i in range(num_watershed_fields):
        assert (subset_watershed_layer_defn.GetFieldDefn(i).GetNameRef() in watershed_field_names)

    #TEST NETWORK
    subset_network_layer_defn = subset_network_layer.GetLayerDefn()
    num_network_fields = subset_network_layer_defn.GetFieldCount()

    network_field_names = ['arcid','from_node','to_node','HydroID','GridID',
                           'NextDownID','SLength','Avg_Slope','LENGTHKM',
                           'Shape_Leng','Musk_x','watershed','subbasin']

    assert (num_network_fields==len(network_field_names))

    for i in range(num_network_fields):
        assert (subset_network_layer_defn.GetFieldDefn(i).GetNameRef() in network_field_names)

    #cleanup
    remove_files(*glob(os.path.join(OUTPUT_DATA_PATH,"DrainageLineSubset2.*")))
    remove_files(*glob(os.path.join(OUTPUT_DATA_PATH,"CatchmentSubset2.*")))

def test_add_length_to_network_taudem():
    """
    Checks adding length to network
    """
    print("TEST 11: TEST ADD LENGTH TO NETWORK")
    td = TauDEM()

    subset_network_file = os.path.join(OUTPUT_DATA_PATH, "DrainageLineSubset2.shp")
    #to extract a specific network
    td.extractSubNetwork(network_file=os.path.join(GIS_INPUT_DATA_PATH, 'u-k', "DrainageLineSubset.shp"),
                         out_subset_network_file=subset_network_file,
                         outlet_ids=[42911], #list of outlet ids
                         river_id_field="HydroID",
                         next_down_id_field="NextDownID",
                         river_magnitude_field="HydroID",
                         safe_mode=False,
                         )

    #add length m field
    td.addLengthMeters(subset_network_file)

    #Test results
    subset_network_shapefile = ogr.Open(subset_network_file)
    subset_network_layer = subset_network_shapefile.GetLayer()

    #make sure all fields are there
    subset_network_layer_defn = subset_network_layer.GetLayerDefn()
    num_network_fields = subset_network_layer_defn.GetFieldCount()

    network_field_names = ['arcid','from_node','to_node','HydroID','GridID',
                           'NextDownID', 'SLength', 'Avg_Slope','LENGTHKM',
                           'Shape_Leng','Musk_x','watershed','subbasin', 'LENGTH_M']

    assert (num_network_fields==len(network_field_names))

    for i in range(num_network_fields):
        assert (subset_network_layer_defn.GetFieldDefn(i).GetNameRef() in network_field_names)


    #make sure values are OK
    length_m_list = array([194.440898134, 601.443392962, 1306.53179652, 1501.27444279,
                           3437.46584922, 5579.56507836, 6347.04650903])
    generated_list = []
    for network_feature in subset_network_layer:
        generated_list.append(network_feature.GetField('LENGTH_M'))

    assert_almost_equal(length_m_list, array(sorted(generated_list)), decimal=2)

    #cleanup
    remove_files(*glob(os.path.join(OUTPUT_DATA_PATH,"DrainageLineSubset2.*")))

@pytest.mark.skipif(not os.path.exists(TAUDEM_EXE_PATH), reason='Only run if TauDEM installed')
def test_generate_network_taudem():
    """
    Checks generate TauDEM network
    """
    print("TEST 12: TEST GENERATE TauDEM NETWORK")
    td = TauDEM(TAUDEM_EXE_PATH, use_all_processors=True)

    elevation_dem = os.path.join(GIS_INPUT_DATA_PATH, 'jamaica_dem.tif')

    td.demToStreamNetwork(OUTPUT_DATA_PATH,
                          elevation_dem,
                          threshold=1000)

    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'pit_filled_elevation_grid.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'pit_filled_elevation_grid.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'flow_dir_grid_d8.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'flow_dir_grid_d8.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'contributing_area_grid_d8.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'contributing_area_grid_d8.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'slope_grid_d8.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'slope_grid_d8.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_raster_grid.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_raster_grid.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_order_grid.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_order_grid.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'network_connectivity_tree.txt')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'network_coordinates.txt')))
#    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_reach_file.shp')))
#    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_reach_file.shx')))
#    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_reach_file.dbf')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_reach_file.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'watershed_grid.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'watershed_grid.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'watershed_shapefile.shp')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'watershed_shapefile.shx')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'watershed_shapefile.dbf')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'watershed_shapefile.prj')))
    #cleanup
    remove_files(*[f for f in glob(os.path.join(OUTPUT_DATA_PATH,"*")) if not f.endswith(".gitignore")])

@pytest.mark.skipif(not os.path.exists(TAUDEM_EXE_PATH), reason='Only run if TauDEM installed')
def test_generate_network_taudem_dinf():
    """
    Checks generate TauDEM network dinf
    """
    print("TEST 13: TEST GENERATE TauDEM NETWORK DINF")
    td = TauDEM(TAUDEM_EXE_PATH)

    elevation_dem = os.path.join(GIS_INPUT_DATA_PATH, 'jamaica_dem.tif')

    td.demToStreamNetwork(OUTPUT_DATA_PATH,
                          pit_filled_elevation_grid=elevation_dem,
                          threshold=1000,
                          use_dinf=True)
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'flow_dir_grid_d8.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'flow_dir_grid_d8.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'flow_dir_grid_dinf.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'flow_dir_grid_dinf.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'contributing_area_grid_d8.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'contributing_area_grid_d8.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'contributing_area_grid_dinf.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'contributing_area_grid_dinf.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'slope_grid_d8.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'slope_grid_d8.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'slope_grid_dinf.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'slope_grid_dinf.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_raster_grid.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_raster_grid.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_order_grid.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_order_grid.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'network_connectivity_tree.txt')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'network_coordinates.txt')))
#    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_reach_file.shp')))
#    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_reach_file.shx')))
#    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_reach_file.dbf')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'stream_reach_file.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'watershed_grid.tif')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'watershed_grid.prj')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'watershed_shapefile.shp')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'watershed_shapefile.shx')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'watershed_shapefile.dbf')))
    assert (os.path.exists(os.path.join(OUTPUT_DATA_PATH, 'watershed_shapefile.prj')))
    #cleanup
    remove_files(*[f for f in glob(os.path.join(OUTPUT_DATA_PATH,"*")) if not f.endswith(".gitignore")])

def test_gen_muskingum_kfac2():
    """
    Checks generating Muskingum Kfac option 2
    """
    print("TEST 14: TEST GENERATE MUSKINGUM KFAC OPTION 2")
    generated_kfac_file = os.path.join(OUTPUT_DATA_PATH,
                                       "kfac2.csv")
    #rapid_connect
    rapid_connect_file = os.path.join(COMPARE_DATA_PATH, "x-x",
                                      "rapid_connect.csv")
    CreateMuskingumKfacFile(in_drainage_line=os.path.join(GIS_INPUT_DATA_PATH, 'flowline.shp'),
                            river_id="COMID",
                            length_id="LENGTHKM",
                            slope_id="Slope",
                            celerity=1000.0/3600.0,
                            formula_type=2,
                            in_connectivity_file=rapid_connect_file,
                            out_kfac_file=generated_kfac_file)

    #CHECK OUTPUT
    #kfac
    generated_kfac_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                "kfac2.csv")
    assert (compare_csv_decimal_files(generated_kfac_file,
                                  generated_kfac_file_solution))

    remove_files(generated_kfac_file)

def test_gen_muskingum_kfac1():
    """
    Checks generating Muskingum Kfac option 1
    """
    print("TEST 14: TEST GENERATE MUSKINGUM KFAC OPTION 1")
    generated_kfac_file = os.path.join(OUTPUT_DATA_PATH,
                                       "kfac1.csv")
    #rapid_connect
    rapid_connect_file = os.path.join(COMPARE_DATA_PATH, "x-x",
                                      "rapid_connect.csv")
    CreateMuskingumKfacFile(in_drainage_line=os.path.join(GIS_INPUT_DATA_PATH, 'flowline.shp'),
                            river_id="COMID",
                            length_id="LENGTHKM",
                            slope_id="Slope",
                            celerity=1000.0/3600.0,
                            formula_type=1,
                            in_connectivity_file=rapid_connect_file,
                            out_kfac_file=generated_kfac_file)

    #CHECK OUTPUT
    #kfac
    generated_kfac_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                "kfac1.csv")
    assert (compare_csv_decimal_files(generated_kfac_file,
                                  generated_kfac_file_solution))
    remove_files(generated_kfac_file)

def test_gen_muskingum_x_drainage():
    """
    Checks generating Muskingum X from draiange line
    """
    print("TEST 15: TEST GENERATE MUSKINGUM X FROM DRAINAGE LINE")
    generated_x_file = os.path.join(OUTPUT_DATA_PATH,
                                    "x_drain.csv")

    CreateMuskingumXFileFromDranageLine(in_drainage_line=os.path.join(GIS_INPUT_DATA_PATH, 'u-k', "DrainageLineSubset.shp"),
                                        x_id="Musk_x",
                                        out_x_file=generated_x_file)

    #CHECK OUTPUT
    generated_x_file_solution = os.path.join(COMPARE_DATA_PATH, "u-k",
                                             "x_drain.csv")
    assert (compare_csv_decimal_files(generated_x_file,
                                  generated_x_file_solution))
    remove_files(generated_x_file)


def test_weight_table_with_invalid_polygon():
    """
    Checks generating weight table with invalid polygon
    """
    print("TEST 16: TEST GENERATE WEIGHT TABLE WITH INVALID POLYGON")
    generated_weight_table_file = os.path.join(OUTPUT_DATA_PATH,
                                               "weight_polygons.csv")
    # rapid_connect
    rapid_connect_file = os.path.join(COMPARE_DATA_PATH, "x-x",
                                      "rapid_connect.csv")

    lsm_grid = os.path.join(LSM_INPUT_DATA_PATH, "gldas2", "GLDAS_NOAH025_3H.A20101231.0000.020.nc4")
    CreateWeightTableLDAS(in_ldas_nc=lsm_grid,
                          in_nc_lon_var="lon",
                          in_nc_lat_var="lat",
                          in_catchment_shapefile=os.path.join(GIS_INPUT_DATA_PATH, 'test_catchments.shp'),
                          river_id="DrainLnID",
                          in_connectivity_file=rapid_connect_file,
                          out_weight_table=generated_weight_table_file)

    generated_weight_table_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                        "weight_polygons.csv")
    assert (compare_csv_decimal_files(generated_weight_table_file,
                                  generated_weight_table_file_solution))

    remove_files(generated_weight_table_file)


def test_weight_table_with_area_id():
    """
    Checks generating weight table with area id
    """
    print("TEST 17: TEST GENERATE WEIGHT TABLE WITH INVALID POLYGON")
    generated_weight_table_file = os.path.join(OUTPUT_DATA_PATH,
                                               "weight_area.csv")
    # rapid_connect
    rapid_connect_file = os.path.join(COMPARE_DATA_PATH, "x-x",
                                      "rapid_connect.csv")

    lsm_grid = os.path.join(LSM_INPUT_DATA_PATH, "gldas2", "GLDAS_NOAH025_3H.A20101231.0000.020.nc4")
    CreateWeightTableLDAS(in_ldas_nc=lsm_grid,
                          in_nc_lon_var="lon",
                          in_nc_lat_var="lat",
                          in_catchment_shapefile=os.path.join(GIS_INPUT_DATA_PATH, 'test_catchments.shp'),
                          river_id="DrainLnID",
                          area_id="Shape_Area",
                          in_connectivity_file=rapid_connect_file,
                          out_weight_table=generated_weight_table_file)

    generated_weight_table_file_solution = os.path.join(COMPARE_DATA_PATH, "x-x",
                                                        "weight_area.csv")
    assert (compare_csv_decimal_files(generated_weight_table_file,
                                  generated_weight_table_file_solution))

    remove_files(generated_weight_table_file)
