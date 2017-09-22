# -*- coding: utf-8 -*-
"""
    weight.py
    RAPIDpy

    Created by Alan D Snow, 2016.
    Based on RAPID_Toolbox for ArcMap
    License: BSD 3-Clause
"""
import csv
from datetime import datetime
from functools import partial

from netCDF4 import Dataset
import numpy as np
from pyproj import Proj, transform
from shapely.wkb import loads as shapely_loads
from shapely.ops import transform as shapely_transform
from shapely.geos import TopologicalError
import rtree  # http://toblerity.org/rtree/install.html
from osgeo import gdal, ogr, osr

# local
from .voronoi import pointsToVoronoiGridArray
from ..helper_functions import log, open_csv

gdal.UseExceptions()


def get_poly_area_geo(poly):
    """
    Calculates the area in meters squared of the individual polygon
    """
    minx, miny, maxx, maxy = poly.bounds
    # reproject polygon to get area
    reprojected_for_area = Proj("+proj=aea +lat_1={0} +lat_1={1} "
                                "+lat_0={2} +lon_0={3}"
                                .format(miny,
                                        maxy,
                                        (miny + maxy) / 2.0,
                                        (minx + maxx) / 2.0))
    geographic_proj = Proj(init='epsg:4326')
    project_func = partial(transform,
                           geographic_proj,
                           reprojected_for_area)
    reprojected_poly = shapely_transform(project_func, poly)
    return reprojected_poly.area


def _get_lat_lon_indices(lsm_lat_array, lsm_lon_array, lat, lon):
    """
    Determines the index in the array (1D or 2D) where the
    lat/lon point is
    """
    if lsm_lat_array.ndim == 2 and lsm_lon_array.ndim == 2:
        lsm_lat_indices_from_lat, lsm_lon_indices_from_lat = \
            np.where((lsm_lat_array == lat))
        lsm_lat_indices_from_lon, lsm_lon_indices_from_lon = \
            np.where((lsm_lon_array == lon))

        index_lsm_grid_lat = np.intersect1d(lsm_lat_indices_from_lat,
                                            lsm_lat_indices_from_lon)[0]
        index_lsm_grid_lon = np.intersect1d(lsm_lon_indices_from_lat,
                                            lsm_lon_indices_from_lon)[0]

    elif lsm_lat_array.ndim == 1 and lsm_lon_array.ndim == 1:
        index_lsm_grid_lon = np.where(lsm_lon_array == lon)[0][0]
        index_lsm_grid_lat = np.where(lsm_lat_array == lat)[0][0]
    else:
        raise IndexError("Lat/Lon lists have invalid dimensions. "
                         "Only 1D or 2D arrays allowed ...")

    return index_lsm_grid_lat, index_lsm_grid_lon


def find_nearest(array, value):
    """
    Get the nearest index to value searching for
    """
    return (np.abs(array-value)).argmin()


def rtree_create_weight_table(lsm_grid_lat, lsm_grid_lon,
                              in_catchment_shapefile, river_id,
                              in_rapid_connect, out_weight_table,
                              file_geodatabase=None, area_id=None):
    """
    Create Weight Table for Land Surface Model Grids
    """
    time_start_all = datetime.utcnow()

    if lsm_grid_lat.ndim == 3 and lsm_grid_lon.ndim == 3:
        # assume first dimension is time
        lsm_grid_lat = lsm_grid_lat[0]
        lsm_grid_lon = lsm_grid_lon[0]

    log("Generating LSM Grid Thiessen Array ...")
    if file_geodatabase:
        gdb_driver = ogr.GetDriverByName("OpenFileGDB")
        ogr_file_geodatabase = gdb_driver.Open(file_geodatabase, 0)
        ogr_catchment_shapefile_lyr = \
            ogr_file_geodatabase.GetLayer(in_catchment_shapefile)
    else:
        ogr_catchment_shapefile = ogr.Open(in_catchment_shapefile)
        ogr_catchment_shapefile_lyr = ogr_catchment_shapefile.GetLayer()

    ogr_catchment_shapefile_lyr_proj = \
        ogr_catchment_shapefile_lyr.GetSpatialRef()
    original_catchment_proj = \
        Proj(ogr_catchment_shapefile_lyr_proj.ExportToProj4())
    geographic_proj = Proj(init='EPSG:4326')
    extent = ogr_catchment_shapefile_lyr.GetExtent()
    if original_catchment_proj != geographic_proj:
        x, y = transform(original_catchment_proj,
                         geographic_proj,
                         [extent[0], extent[1]],
                         [extent[2], extent[3]])
        extent = [min(x), max(x), min(y), max(y)]

    lsm_grid_feature_list = \
        pointsToVoronoiGridArray(lsm_grid_lat, lsm_grid_lon, extent)

#    ##COMMENTED LINES FOR TESTING
#    import os
#    from .voronoi import pointsToVoronoiGridShapefile
#    vor_shp_path = \
#        os.path.join(os.path.dirname(in_catchment_shapefile), "test_grid.shp")
#    pointsToVoronoiGridShapefile(lsm_grid_lat, lsm_grid_lon,
#                                 vor_shp_path, extent)

    time_end_lsm_grid_thiessen = datetime.utcnow()
    log(time_end_lsm_grid_thiessen - time_start_all)

    log("Generating LSM Grid Rtree ...")
    rtree_idx = rtree.index.Index()
    # Populate R-tree index with bounds of ECMWF grid cells
    for lsm_grid_pos, lsm_grid_feature in enumerate(lsm_grid_feature_list):
        rtree_idx.insert(lsm_grid_pos, lsm_grid_feature['polygon'].bounds)

    time_end_lsm_grid_rtree = datetime.utcnow()
    log(time_end_lsm_grid_rtree - time_end_lsm_grid_thiessen)

    log("Retrieving catchment river id list ...")
    number_of_catchment_features = \
        ogr_catchment_shapefile_lyr.GetFeatureCount()
    catchment_rivid_list = \
        np.zeros(number_of_catchment_features, dtype=np.int32)
    for feature_idx, catchment_feature in \
            enumerate(ogr_catchment_shapefile_lyr):
        catchment_rivid_list[feature_idx] = \
            catchment_feature.GetField(river_id)

    log("Reading in RAPID connect file ...")
    rapid_connect_rivid_list = np.loadtxt(in_rapid_connect,
                                          delimiter=",",
                                          usecols=(0,),
                                          ndmin=1,
                                          dtype=int)
    log("Find LSM grid cells that intersect with each catchment")
    log("and write out weight table ...")

    dummy_lat_index, dummy_lon_index = \
        _get_lat_lon_indices(lsm_grid_lat,
                             lsm_grid_lon,
                             lsm_grid_feature_list[0]['lat'],
                             lsm_grid_feature_list[0]['lon'])
    dummy_row_end = [
        0,
        dummy_lon_index,
        dummy_lat_index,
        1,
        lsm_grid_feature_list[0]['lon'],
        lsm_grid_feature_list[0]['lat']
    ]

    with open_csv(out_weight_table, 'w') as csvfile:
        connectwriter = csv.writer(csvfile)
        connectwriter.writerow(['rivid', 'area_sqm', 'lon_index', 'lat_index',
                                'npoints', 'lsm_grid_lon', 'lsm_grid_lat'])
        geographic_proj = Proj(init='EPSG:4326')
        osr_geographic_proj = osr.SpatialReference()
        osr_geographic_proj.ImportFromEPSG(4326)
        proj_transform = None
        if original_catchment_proj != geographic_proj:
            proj_transform = \
                osr.CoordinateTransformation(ogr_catchment_shapefile_lyr_proj,
                                             osr_geographic_proj)

        for rapid_connect_rivid in rapid_connect_rivid_list:
            intersect_grid_info_list = []
            try:
                catchment_pos = \
                    np.where(catchment_rivid_list == rapid_connect_rivid)[0][0]
            except IndexError:
                # if it is not in the catchment, add dummy row in its place
                connectwriter.writerow([rapid_connect_rivid] + dummy_row_end)
                continue

            get_catchment_feature = \
                ogr_catchment_shapefile_lyr.GetFeature(catchment_pos)
            feat_geom = get_catchment_feature.GetGeometryRef()
            # make sure coordinates are geographic
            if proj_transform:
                feat_geom.Transform(proj_transform)
            catchment_polygon = shapely_loads(feat_geom.ExportToWkb())

            for sub_lsm_grid_pos in \
                    rtree_idx.intersection(catchment_polygon.bounds):
                lsm_grid_polygon = \
                    lsm_grid_feature_list[sub_lsm_grid_pos]['polygon']
                if catchment_polygon.intersects(lsm_grid_polygon):
                    try:
                        intersect_poly = \
                            catchment_polygon.intersection(lsm_grid_polygon)
                    except TopologicalError:
                        log('The catchment polygon with id {0} was '
                            'invalid. Attempting to self clean...'
                            .format(rapid_connect_rivid))
                        original_area = catchment_polygon.area
                        catchment_polygon = catchment_polygon.buffer(0)
                        area_ratio = original_area/catchment_polygon.area
                        log('AREA_RATIO: {0}'.format(area_ratio))
                        msg_level = "INFO"
                        if round(area_ratio, 5) != 1:
                            msg_level = "WARNING"
                        log('The cleaned catchment polygon area '
                            'differs from the original area by {1}%.'
                            .format(abs(area_ratio - 1)), severity=msg_level)
                        intersect_poly = \
                            catchment_polygon.intersection(lsm_grid_polygon)
                    if not area_id:
                        # attempt to calculate AREA
                        poly_area = get_poly_area_geo(intersect_poly)
                    else:
                        poly_area = \
                            float(get_catchment_feature.GetField(area_id)) * \
                            intersect_poly.area/catchment_polygon.area

                    index_lsm_grid_lat, index_lsm_grid_lon = \
                        _get_lat_lon_indices(
                            lsm_grid_lat,
                            lsm_grid_lon,
                            lsm_grid_feature_list[sub_lsm_grid_pos]['lat'],
                            lsm_grid_feature_list[sub_lsm_grid_pos]['lon'])
                    intersect_grid_info_list.append({
                        'rivid': rapid_connect_rivid,
                        'area': poly_area,
                        'lsm_grid_lat':
                            lsm_grid_feature_list[sub_lsm_grid_pos]['lat'],
                        'lsm_grid_lon':
                            lsm_grid_feature_list[sub_lsm_grid_pos]['lon'],
                        'index_lsm_grid_lon': index_lsm_grid_lon,
                        'index_lsm_grid_lat': index_lsm_grid_lat
                    })

            npoints = len(intersect_grid_info_list)
            # If no intersection found, add dummy row
            if npoints <= 0:
                connectwriter.writerow([rapid_connect_rivid] + dummy_row_end)

            for intersect_grid_info in intersect_grid_info_list:
                connectwriter.writerow([
                    intersect_grid_info['rivid'],
                    intersect_grid_info['area'],
                    intersect_grid_info['index_lsm_grid_lon'],
                    intersect_grid_info['index_lsm_grid_lat'],
                    npoints,
                    intersect_grid_info['lsm_grid_lon'],
                    intersect_grid_info['lsm_grid_lat']
                ])

    time_end_all = datetime.utcnow()
    log(time_end_all - time_end_lsm_grid_rtree)
    log("TOTAL TIME: {0}".format(time_end_all - time_start_all))


def CreateWeightTableECMWF(in_ecmwf_nc,
                           in_catchment_shapefile,
                           river_id,
                           in_connectivity_file,
                           out_weight_table,
                           area_id=None,
                           file_geodatabase=None):
    """
    Create Weight Table for ECMWF Grids

    .. note:: The grids are in the RAPIDpy package under
              the gis/lsm_grids folder.

    Parameters
    ----------
    in_ecmwf_nc: str
        Path to the ECMWF NetCDF grid.
    in_catchment_shapefile: str
        Path to the Catchment shapefile.
    river_id: str
        The name of the field with the river ID (Ex. 'DrainLnID' or 'LINKNO').
    in_connectivity_file: str
        The path to the RAPID connectivity file.
    out_weight_table: str
        The path to the output weight table file.
    area_id: str, optional
        The name of the field with the area of each catchment stored in meters
        squared. Default is it calculate the area.
    file_geodatabase: str, optional
        Path to the file geodatabase. If you use this option, in_drainage_line
        is the name of the stream network feature class.
        (WARNING: Not always stable with GDAL.)


    Example:

    .. code:: python

        from RAPIDpy.gis.weight import CreateWeightTableECMWF

        CreateWeightTableECMWF(
            in_ecmwf_nc='/path/to/runoff_ecmwf_grid.nc'
            in_catchment_shapefile='/path/to/catchment.shp',
            river_id='LINKNO',
            in_connectivity_file='/path/to/rapid_connect.csv',
            out_weight_table='/path/to/ecmwf_weight.csv',
        )

    """
    # extract ECMWF GRID
    data_ecmwf_nc = Dataset(in_ecmwf_nc)
    variables_list = data_ecmwf_nc.variables.keys()
    in_ecmwf_lat_var = 'lat'
    if 'latitude' in variables_list:
        in_ecmwf_lat_var = 'latitude'
    in_ecmwf_lon_var = 'lon'
    if 'longitude' in variables_list:
        in_ecmwf_lon_var = 'longitude'

    # convert [0, 360] to [-180, 180]
    ecmwf_lon = \
        (data_ecmwf_nc.variables[in_ecmwf_lon_var][:] + 180) % 360 - 180
    # assume [-90, 90]
    ecmwf_lat = data_ecmwf_nc.variables[in_ecmwf_lat_var][:]
    data_ecmwf_nc.close()

    rtree_create_weight_table(ecmwf_lat, ecmwf_lon,
                              in_catchment_shapefile, river_id,
                              in_connectivity_file, out_weight_table,
                              file_geodatabase, area_id)


def CreateWeightTableLDAS(in_ldas_nc,
                          in_nc_lon_var,
                          in_nc_lat_var,
                          in_catchment_shapefile,
                          river_id,
                          in_connectivity_file,
                          out_weight_table,
                          area_id=None,
                          file_geodatabase=None):
    """
    Create Weight Table for NLDAS, GLDAS grids as well as for 2D Joules,
    or LIS Grids

    Parameters
    ----------
    in_ldas_nc: str
        Path to the land surface model NetCDF grid.
    in_nc_lon_var: str
        The variable name in the NetCDF file for the longitude.
    in_nc_lat_var: str
        The variable name in the NetCDF file for the latitude.
    in_catchment_shapefile: str
        Path to the Catchment shapefile.
    river_id: str
        The name of the field with the river ID (Ex. 'DrainLnID' or 'LINKNO').
    in_connectivity_file: str
        The path to the RAPID connectivity file.
    out_weight_table: str
        The path to the output weight table file.
    area_id: str, optional
        The name of the field with the area of each catchment stored in meters
        squared. Default is it calculate the area.
    file_geodatabase: str, optional
        Path to the file geodatabase. If you use this option, in_drainage_line
        is the name of the stream network feature class.
        (WARNING: Not always stable with GDAL.)


    Example:

    .. code:: python

        from RAPIDpy.gis.weight import CreateWeightTableLDAS

        CreateWeightTableLDAS(
            in_ldas_nc='/path/to/runoff_grid.nc',
            in_nc_lon_var="lon_110",
            in_nc_lat_var="lat_110",
            in_catchment_shapefile='/path/to/catchment.shp',
            river_id='LINKNO',
            in_connectivity_file='/path/to/rapid_connect.csv',
            out_weight_table='/path/to/ldas_weight.csv',
        )
    """
    # extract LDAS GRID
    data_ldas_nc = Dataset(in_ldas_nc)
    variables_list = data_ldas_nc.variables.keys()
    if in_nc_lon_var not in variables_list:
        raise Exception("Invalid longitude variable. Choose from: {0}"
                        .format(variables_list))
    if in_nc_lat_var not in variables_list:
        raise Exception("Invalid latitude variable. Choose from: {0}"
                        .format(variables_list))
    ldas_lon = data_ldas_nc.variables[in_nc_lon_var][:]  # assume [-180, 180]
    ldas_lat = data_ldas_nc.variables[in_nc_lat_var][:]  # assume [-90,90]
    data_ldas_nc.close()

    rtree_create_weight_table(ldas_lat, ldas_lon,
                              in_catchment_shapefile, river_id,
                              in_connectivity_file, out_weight_table,
                              file_geodatabase, area_id)
