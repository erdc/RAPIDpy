# -*- coding: utf-8 -*-
"""
    workflow.py
    RAPIDpy

    Created by Alan D Snow, 2016.
    Based on RAPID_Toolbox for ArcMap
    License: BSD 3-Clause
"""
import os

from .network import (CreateNetworkConnectivity,
                      CreateNetworkConnectivityTauDEMTree,
                      CreateNetworkConnectivityNHDPlus,
                      CreateSubsetFile)
from .muskingum import (CreateMuskingumKfacFile, CreateMuskingumKFile,
                        CreateConstMuskingumXFile)
from .weight import CreateWeightTableECMWF
from .centroid import FlowlineToPoint


def CreateAllStaticRAPIDFiles(in_drainage_line,
                              river_id,
                              length_id,
                              slope_id,
                              next_down_id,
                              rapid_output_folder,
                              kfac_celerity=1000.0/3600.0,
                              kfac_formula_type=3,
                              kfac_length_units="km",
                              lambda_k=0.35,
                              x_value=0.3,
                              nhdplus=False,
                              taudem_network_connectivity_tree_file=None,
                              file_geodatabase=None):
    """
    To generate the static RAPID files (rapid_connect.csv, riv_bas_id.csv,
    kfac.csv, k.csv, x.csv, comid_lat_lon_z.csv) with default values.

    Parameters
    ----------
    in_drainage_line: str
        Path to the stream network (i.e. Drainage Line) shapefile.
    river_id: str
        The name of the field with the river ID
        (Ex. 'HydroID', 'COMID', or 'LINKNO').
    length_id: str
        The field name containging the length of the river segment
        (Ex. 'LENGTHKM' or 'Length').
    slope_id: str
        The field name containging the slope of the river segment
        (Ex. 'Avg_Slope' or 'Slope').
    next_down_id: str
        The name of the field with the river ID of the next downstream river
        segment (Ex. 'NextDownID' or 'DSLINKNO').
    rapid_output_folder: str
        The path to the folder where all of the RAPID output will be generated.
    kfac_celerity: float, optional
        The flow wave celerity for the watershed in meters per second.
        1 km/hr or 1000.0/3600.0 m/s is a reasonable value if unknown.
    kfac_formula_type: int, optional
        An integer representing the formula type to use when calculating kfac.
        Default is 3.
    kfac_length_units: str, optional
        The units for the length_id field. Supported types are "m" for meters
        and "km" for kilometers. Default is "km".
    lambda_k: float, optional
        The value for lambda given from RAPID after the calibration process.
        Default is 0.35.
    x_value: float, optional
        Value for the muskingum X parameter [0-0.5]. Default is 0.3.
    nhdplus: bool, optional
        If True, the drainage line is from the NHDPlus dataset with the VAA
        fields COMID, FROMNODE, TONODE, and DIVERGENCE. Default is False.
    taudem_network_connectivity_tree_file: str, optional
        If set, the connectivity file will be generated from the TauDEM
        connectivity tree file.
    file_geodatabase: str, optional
        Path to the file geodatabase. If you use this option,
        in_drainage_line is the name of the stream network feature class.
        (WARNING: Not always stable with GDAL.)


    Example::

        from RAPIDpy.gis.workflow import CreateAllStaticRAPIDFiles

        CreateAllStaticRAPIDFiles(
            in_drainage_line="/path/to/drainage_line.shp",
            river_id="HydroID",
            length_id="LENGTHKM",
            slope_id="SLOPE",
            next_down_river_id="NextDownID",
            rapid_output_folder="/path/to/rapid/output",
        )
    """
    # RAPID connect file
    rapid_connect_file = os.path.join(rapid_output_folder, 'rapid_connect.csv')
    if nhdplus:
        CreateNetworkConnectivityNHDPlus(in_drainage_line,
                                         rapid_connect_file,
                                         file_geodatabase)
    elif taudem_network_connectivity_tree_file:
        CreateNetworkConnectivityTauDEMTree(
            taudem_network_connectivity_tree_file,
            rapid_connect_file)
    else:
        CreateNetworkConnectivity(in_drainage_line,
                                  river_id,
                                  next_down_id,
                                  rapid_connect_file,
                                  file_geodatabase)

    # river basin id file
    riv_bas_id_file = os.path.join(rapid_output_folder, 'riv_bas_id.csv')
    CreateSubsetFile(in_drainage_line,
                     river_id,
                     riv_bas_id_file,
                     file_geodatabase)
    # kfac file
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
    # k file
    k_file = os.path.join(rapid_output_folder, 'k.csv')
    CreateMuskingumKFile(lambda_k,
                         kfac_file,
                         k_file)
    # x file
    x_file = os.path.join(rapid_output_folder, 'x.csv')
    CreateConstMuskingumXFile(x_value,
                              rapid_connect_file,
                              x_file)
    # comid lat lon z file
    comid_lat_lon_z_file = \
        os.path.join(rapid_output_folder, 'comid_lat_lon_z.csv')
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
    This creates all of the ECMWF grid weight tables using an area
    weighted method based on Esri's RAPID_Toolbox.

    Parameters
    ----------
    in_catchment: str
        Path to the Catchment shapefile.
    catchment_river_id: str
        The name of the field with the river ID (Ex. 'DrainLnID' or 'LINKNO').
    rapid_output_folder: str
        The path to the folder where all of the RAPID output will be generated.
    rapid_connect_file: str
        The path to the RAPID connectivity file.
    file_geodatabase: str, optional
        Path to the file geodatabase. If you use this option,
        in_drainage_line is the name of the stream network feature class.
        (WARNING: Not always stable with GDAL.)


    Example::

        from RAPIDpy.gis.workflow import CreateAllStaticECMWFFiles

        CreateAllStaticECMWFFiles(
            in_catchment="/path/to/catchment.shp",
            catchment_river_id="DrainLnID",
            rapid_output_folder="/path/to/rapid/output",
            rapid_connect_file="/path/to/rapid_connect.csv",
        )

    """
    lsm_grid_folder = \
        os.path.join(os.path.dirname(os.path.realpath(__file__)), 'lsm_grids')

    # create from ECMWF high reslution grid
    ecmwf_t1279_grid_file = \
        os.path.join(lsm_grid_folder, 'runoff_ecmwf_t1279_grid.nc')
    weight_ecmwf_t1279_file = \
        os.path.join(rapid_output_folder, 'weight_ecmwf_t1279.csv')
    CreateWeightTableECMWF(ecmwf_t1279_grid_file,
                           in_catchment,
                           catchment_river_id,
                           rapid_connect_file,
                           weight_ecmwf_t1279_file,
                           file_geodatabase=file_geodatabase)

    # create from ECMWF low reslution grid
    ecmwf_tco639_grid_file = \
        os.path.join(lsm_grid_folder, 'runoff_ecmwf_tco639_grid.nc')
    weight_ecmwf_tco639_file = \
        os.path.join(rapid_output_folder, 'weight_ecmwf_tco639.csv')
    CreateWeightTableECMWF(ecmwf_tco639_grid_file,
                           in_catchment,
                           catchment_river_id,
                           rapid_connect_file,
                           weight_ecmwf_tco639_file,
                           file_geodatabase=file_geodatabase)

    # create from ERA Interim grid
    era_t511_grid_file = \
        os.path.join(lsm_grid_folder, 'runoff_era_t511_grid.nc')
    weight_era_t511_file = \
        os.path.join(rapid_output_folder, 'weight_era_t511.csv')
    CreateWeightTableECMWF(era_t511_grid_file,
                           in_catchment,
                           catchment_river_id,
                           rapid_connect_file,
                           weight_era_t511_file,
                           file_geodatabase=file_geodatabase)


def CreateAllStaticECMWFRAPIDFiles(in_drainage_line,
                                   river_id,
                                   length_id,
                                   slope_id,
                                   next_down_id,
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
                                   file_geodatabase=None):
    """
    This creates all of the static RAPID files and ECMWF grid weight tables.

    Parameters
    ----------
    in_drainage_line: str
        Path to the stream network (i.e. Drainage Line) shapefile.
    river_id: str
        The name of the field with the river ID
        (Ex. 'HydroID', 'COMID', or 'LINKNO').
    length_id: str
        The field name containging the length of the river segment
        (Ex. 'LENGTHKM' or 'Length').
    slope_id: str
        The field name containging the slope of the river segment
        (Ex. 'Avg_Slope' or 'Slope').
    next_down_id: str
        The name of the field with the river ID of the next downstream
        river segment (Ex. 'NextDownID' or 'DSLINKNO').
    in_catchment: str
        Path to the Catchment shapefile.
    catchment_river_id: str
        The name of the field with the river ID (Ex. 'DrainLnID' or 'LINKNO').
    rapid_output_folder: str
        The path to the folder where all of the RAPID output will be generated.
    kfac_celerity: float, optional
        The flow wave celerity for the watershed in meters per second.
        1 km/hr or 1000.0/3600.0 m/s is a reasonable value if unknown.
    kfac_formula_type: int, optional
        An integer representing the formula type to use when calculating kfac.
        Default is 3.
    kfac_length_units: str, optional
        The units for the length_id field. Supported types are "m" for meters
        and "km" for kilometers. Default is "km".
    lambda_k: float, optional
        The value for lambda given from RAPID after the calibration process.
        Default is 0.35.
    x_value: float, optional
        Value for the muskingum X parameter [0-0.5].Default is 0.3.
    nhdplus: bool, optional
        If True, the drainage line is from the NHDPlus dataset with the
        VAA fields COMID, FROMNODE, TONODE, and DIVERGENCE. Default is False.
    taudem_network_connectivity_tree_file: str, optional
        If set, the connectivity file will be generated from the
        TauDEM connectivity tree file.
    file_geodatabase: str, optional
        Path to the file geodatabase. If you use this option,
        in_drainage_line is the name of the stream network feature class
        (WARNING: Not always stable with GDAL).


    Example::

        from RAPIDpy.gis.workflow import CreateAllStaticECMWFRAPIDFiles

        CreateAllStaticECMWFRAPIDFiles(
            in_drainage_line="/path/to/drainage_line.shp",
            river_id="HydroID",
            length_id="LENGTHKM",
            slope_id="SLOPE",
            next_down_id="NextDownID",
            in_catchment="/path/to/catchment.shp",
            catchment_river_id="DrainLnID",
            rapid_output_folder="/path/to/rapid/output",
        )
    """
    CreateAllStaticRAPIDFiles(in_drainage_line,
                              river_id,
                              length_id,
                              slope_id,
                              next_down_id,
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
