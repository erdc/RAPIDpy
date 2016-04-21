# -*- coding: utf-8 -*-
##
##  workflow.py
##  RAPIDpy
##
##  Created by Alan D Snow.
##  Based on RAPID_Toolbox for ArcMap
##  Copyright © 2016 Alan D Snow. All rights reserved.
##  License: BSD 3-Clause

import os

from .network import (CreateNetworkConnectivity, CreateNetworkConnectivityTauDEMTree,
                      CreateNetworkConnectivityNHDPlus, CreateSubsetFile)
from .muskingum import (CreateMuskingumKfacFile, CreateMuskingumKFile, 
                        CreateConstMuskingumXFile)
from .weight import CreateWeightTableECMWF
from .centroid import FlowlineToPoint

def CreateAllStaticRAPIDFiles(in_drainage_line,
                              river_id,
                              length_id,
                              slope_id,
                              next_down_river_id,
                              in_catchment,
                              catchment_river_id,
                              rapid_output_folder,
                              kfac_celerity=1000.0/3600.0,
                              kfac_formula_type=3,
                              kfac_length_units="km",
                              lambda_k=0.35,
                              x_value=0.3,
                              nhdplus=False,
                              taudem_network_connectivity_tree_file=None,
                              file_geodatabase=None
                             ):
    """
    Create All Static RAPID Input Files
    """
    #RAPID connect file
    rapid_connect_file = os.path.join(rapid_output_folder, 'rapid_connect.csv')
    if nhdplus:
        CreateNetworkConnectivityNHDPlus(in_drainage_line,
                                         rapid_connect_file,
                                         file_geodatabase)
    elif taudem_network_connectivity_tree_file:
        CreateNetworkConnectivityTauDEMTree(taudem_network_connectivity_tree_file,
                                            rapid_connect_file)
    else:
        CreateNetworkConnectivity(in_drainage_line,
                                  river_id,
                                  next_down_river_id,
                                  rapid_connect_file,
                                  file_geodatabase)
        
    #river basin id file                          
    riv_bas_id_file = os.path.join(rapid_output_folder, 'riv_bas_id.csv')
    CreateSubsetFile(in_drainage_line,
                     river_id, 
                     riv_bas_id_file,
                     file_geodatabase)
    #kfac file                         
    kfac_file = os.path.join(rapid_output_folder, 'kfac.csv')
    CreateMuskingumKfacFile(in_drainage_line,
                            river_id,
                            length_id,
                            slope_id,
                            kfac_celerity,
                            kfac_formula_type,
                            rapid_connect_file,
                            kfac_file,
                            length_units=kfac_length_units,
                            file_geodatabase=file_geodatabase)
    #k file                        
    k_file = os.path.join(rapid_output_folder, 'k.csv')
    CreateMuskingumKFile(lambda_k,
                         kfac_file,
                         k_file)
    #x file
    x_file = os.path.join(rapid_output_folder, 'x.csv')
    CreateConstMuskingumXFile(x_value,
                              rapid_connect_file,
                              x_file)
    #comid lat lon z file
    comid_lat_lon_z_file = os.path.join(rapid_output_folder, 'comid_lat_lon_z.csv')
    FlowlineToPoint(in_drainage_line,
                    river_id,
                    comid_lat_lon_z_file,
                    file_geodatabase)
                           
def CreateAllStaticECMWFFiles(in_catchment,
                              catchment_river_id,
                              rapid_output_folder,
                              rapid_connect_file,
                              file_geodatabase=None
                              ):
    """
    This creates all of the ECMWF grid weight tables
    """
    lsm_grid_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'lsm_grids')    

    #create from ECMWF high reslution grid
    ecmwf_t1279_grid_file = os.path.join(lsm_grid_folder, 'runoff_ecmwf_t1279_grid.nc')
    weight_ecmwf_t1279_file = os.path.join(rapid_output_folder, 'weight_ecmwf_t1279.csv')
    CreateWeightTableECMWF(ecmwf_t1279_grid_file, 
                           in_catchment, 
                           catchment_river_id,
                           rapid_connect_file, 
                           weight_ecmwf_t1279_file,
                           file_geodatabase)

    #create from ECMWF low reslution grid
    ecmwf_tco639_grid_file = os.path.join(lsm_grid_folder, 'runoff_ecmwf_tco639_grid.nc')
    weight_ecmwf_tco639_file = os.path.join(rapid_output_folder, 'weight_ecmwf_tco639.csv')
    CreateWeightTableECMWF(ecmwf_tco639_grid_file, 
                           in_catchment, 
                           catchment_river_id,
                           rapid_connect_file, 
                           weight_ecmwf_tco639_file,
                           file_geodatabase)

    #create from ERA Interim grid
    era_t511_grid_file = os.path.join(lsm_grid_folder, 'runoff_era_t511_grid.nc')
    weight_era_t511_file = os.path.join(rapid_output_folder, 'weight_era_t511.csv')
    CreateWeightTableECMWF(era_t511_grid_file, 
                           in_catchment, 
                           catchment_river_id,
                           rapid_connect_file, 
                           weight_era_t511_file,
                           file_geodatabase)

def CreateAllStaticECMWFRAPIDFiles(in_drainage_line,
                                   river_id,
                                   length_id,
                                   slope_id,
                                   next_down_river_id,
                                   in_catchment,
                                   catchment_river_id,
                                   rapid_output_folder,
                                   kfac_celerity=1000.0/3600.0,
                                   kfac_formula_type=3,
                                   kfac_length_units="km",
                                   lambda_k=0.35,
                                   x_value=0.3,
                                   nhdplus=False,
                                   taudem_network_connectivity_tree_file=None,
                                   file_geodatabase=None
                                   ):
    """
    This creates all of the static rapid files and ECMWF grid weight tables
    """
    #create all RAPID files
    CreateAllStaticRAPIDFiles(in_drainage_line,
                              river_id,
                              length_id,
                              slope_id,
                              next_down_river_id,
                              in_catchment,
                              catchment_river_id,
                              rapid_output_folder,
                              kfac_celerity,
                              kfac_formula_type,
                              kfac_length_units,
                              lambda_k,
                              x_value,
                              nhdplus,
                              taudem_network_connectivity_tree_file,
                              file_geodatabase)
                        

    rapid_connect_file = os.path.join(rapid_output_folder, 'rapid_connect.csv')

    CreateAllStaticECMWFFiles(in_catchment,
                              catchment_river_id,
                              rapid_output_folder,
                              rapid_connect_file,
                              file_geodatabase)