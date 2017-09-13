# -*- coding: utf-8 -*-
##
##  weight.py
##  RAPIDpy
##
##  Created by Alan D Snow.
##  Based on RAPID_Toolbox for ArcMap
##  Copyright © 2016 Alan D Snow. All rights reserved.
##  License: BSD 3-Clause

import csv
from datetime import datetime
from functools import partial
from netCDF4 import Dataset
import pangaea as pa
import numpy as np

try:
    from pyproj import Proj, transform
    from shapely.wkb import loads as shapely_loads
    from shapely.ops import transform as shapely_transform
    from shapely.geos import TopologicalError
    from shapely.geometry import Polygon
    import rtree #http://toblerity.org/rtree/install.html
    from osgeo import gdal, ogr, osr
except Exception:
    raise Exception("You need the gdal, pyproj, shapely, and rtree python package to run these tools ...")

#local
from .voronoi import pointsToVoronoiGridArray, pointsArrayToPolygonGridList
from ..helper_functions import open_csv

GRIDDING_METHODS = {
    'array_math': pointsArrayToPolygonGridList,
    'voronoi': pointsToVoronoiGridArray
}

gdal.UseExceptions()


def gdal_error_handler(err_class, err_num, err_msg):
    errtype = {
            gdal.CE_None: 'None',
            gdal.CE_Debug: 'Debug',
            gdal.CE_Warning: 'Warning',
            gdal.CE_Failure: 'Failure',
            gdal.CE_Fatal: 'Fatal'
    }
    err_msg = err_msg.replace('\n', ' ')
    err_class = errtype.get(err_class, 'None')
    print('Error Number: %s' % (err_num,))
    print('Error Type: %s' % (err_class,))
    print('Error Message: %s' % (err_msg,))

# install error handler
gdal.PushErrorHandler(gdal_error_handler)


def get_poly_area_geo(poly):
    """
    Calculates the area in meters squared of the individual polygon
    """
    minx, miny, maxx, maxy = poly.bounds
    #reproject polygon to get area
    reprojected_for_area = Proj("+proj=aea +lat_1={0} +lat_1={1} +lat_0={2} +lon_0={3}".format(miny,
                                                                                               maxy,
                                                                                               (miny+maxy)/2.0,
                                                                                               (minx+maxx)/2.0))
    geographic_proj = Proj(init='epsg:4326')
    project_func = partial(transform,
                           geographic_proj,
                           reprojected_for_area)
    reprojected_poly = shapely_transform(project_func, poly)
    return reprojected_poly.area
   

def find_nearest(array, value):
    """
    Get the nearest index to value searching for
    """
    return (np.abs(array-value)).argmin()


def RTreeCreateWeightTable(lsm_grid_lat, lsm_grid_lon, 
                           in_catchment_shapefile, river_id,
                           in_rapid_connect, out_weight_table,
                           file_geodatabase=None, area_id=None,
                           grid_proj=None, no_data_value=-9999,
                           gridding_method='array_math'
                           ):
                                      
    """
    Create Weight Table for Land Surface Model Grids
    """
    time_start_all = datetime.utcnow()    

    if lsm_grid_lat.ndim == 3 and lsm_grid_lon.ndim == 3:
        #assume first dimension is time
        lsm_grid_lat = lsm_grid_lat[0]
        lsm_grid_lon = lsm_grid_lon[0]
    
    print("Generating LSM Grid Thiessen Array ...")
    if file_geodatabase:
        gdb_driver = ogr.GetDriverByName("OpenFileGDB")
        ogr_file_geodatabase = gdb_driver.Open(file_geodatabase, 0)
        ogr_catchment_shapefile_lyr = ogr_file_geodatabase.GetLayer(in_catchment_shapefile)
    else:
        ogr_catchment_shapefile = ogr.Open(in_catchment_shapefile)
        ogr_catchment_shapefile_lyr = ogr_catchment_shapefile.GetLayer()

    osr_geographic_proj = osr.SpatialReference()
    osr_geographic_proj.ImportFromEPSG(4326)

    ogr_catchment_proj = ogr_catchment_shapefile_lyr.GetSpatialRef()

    extent = ogr_catchment_shapefile_lyr.GetExtent()

    catchment_transform = None
    if ogr_catchment_proj != osr_geographic_proj:
        catchment_transform = osr.CoordinateTransformation(ogr_catchment_proj, osr_geographic_proj)

        point1, point2 = catchment_transform.TransformPoints(((extent[0], extent[2]), (extent[1], extent[3])))
        x = point1[0], point2[0]
        y = point1[1], point2[1]

        extent = [min(x), max(x), min(y), max(y)]

        poly = Polygon([(extent[0], extent[2]),
                        (extent[1], extent[2]),
                        (extent[1], extent[3]),
                        (extent[0], extent[3]),
                        ])

        ogr_poly = ogr.CreateGeometryFromWkb(poly.wkb)
        ogr_poly.Transform(catchment_transform)

        extent = ogr_poly.GetEnvelope()

    grid_transform_in = None
    grid_transform_out = None

    if grid_proj is not None:
        grid_transform_in = osr.CoordinateTransformation(osr_geographic_proj, grid_proj)
        grid_transform_out = osr.CoordinateTransformation(grid_proj, osr_geographic_proj)

        extent = None
        # TODO figure out how to define the extents if grid needs to be reprojected

    generate_grid_function = GRIDDING_METHODS[gridding_method]
    lsm_grid_feature_list = generate_grid_function(lsm_grid_lat, lsm_grid_lon,
                                                     extent, grid_transform_in, grid_transform_out,
                                                     no_data_value)
    
    ##COMMENTED LINES FOR TESTING
#    import os
#    from .voronoi import pointsToVoronoiGridShapefile
#    vor_shp_path = os.path.join(os.path.dirname(in_catchment_shapefile), "test_grid.shp")
#    pointsToVoronoiGridShapefile(lsm_grid_lat, lsm_grid_lon, vor_shp_path, extent)
    
    time_end_lsm_grid_thiessen = datetime.utcnow()
    print(time_end_lsm_grid_thiessen - time_start_all)
    
    print("Generating LSM Grid Rtree ...")
    rtree_idx = rtree.index.Index()
    # Populate R-tree index with bounds of ECMWF grid cells
    for lsm_grid_pos, lsm_grid_feature in enumerate(lsm_grid_feature_list):
        rtree_idx.insert(lsm_grid_pos, lsm_grid_feature['polygon'].bounds)
     
    time_end_lsm_grid_rtree = datetime.utcnow()
    print(time_end_lsm_grid_rtree - time_end_lsm_grid_thiessen)
    
    print("Retrieving catchment river id list ...")
    number_of_catchment_features = ogr_catchment_shapefile_lyr.GetFeatureCount()
    catchment_rivid_list = np.zeros(number_of_catchment_features, dtype=np.int64)
    for feature_idx, catchment_feature in enumerate(ogr_catchment_shapefile_lyr):
        catchment_rivid_list[feature_idx] = catchment_feature.GetField(river_id)
        
    print("Reading in RAPID connect file ...")
    rapid_connect_rivid_list = np.loadtxt(in_rapid_connect, 
                                          delimiter=",", 
                                          usecols=(0,),
                                          ndmin=1,
                                          dtype=int)
    print("Find LSM grid cells that intersect with each catchment")
    print("and write out weight table ...")

    dummy_grid_feature = lsm_grid_feature_list[0]
    dummy_lat_index = dummy_grid_feature['lat_index']
    dummy_lon_index = dummy_grid_feature['lon_index']
    dummy_row_end = [0,
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
            
        for rapid_connect_rivid in rapid_connect_rivid_list:
            intersect_grid_info_list = []
            try:
                catchment_pos = np.where(catchment_rivid_list==rapid_connect_rivid)[0][0]
            except IndexError:
                #if it is not in the catchment, add dummy row in its place
                connectwriter.writerow([rapid_connect_rivid] + dummy_row_end)
                continue
                pass
            get_catchment_feature = ogr_catchment_shapefile_lyr.GetFeature(catchment_pos)
            feat_geom = get_catchment_feature.GetGeometryRef()
            #make sure coordinates are geographic
            if catchment_transform:
                feat_geom.Transform(catchment_transform)
            catchment_polygon = shapely_loads(feat_geom.ExportToWkb())

            for sub_lsm_grid_pos in rtree_idx.intersection(catchment_polygon.bounds):
                lsm_grid_polygon = lsm_grid_feature_list[sub_lsm_grid_pos]['polygon']
                if catchment_polygon.intersects(lsm_grid_polygon):
                    try:
                        intersect_poly = catchment_polygon.intersection(lsm_grid_polygon)
                    except TopologicalError:
                        print('INFO: The catchment polygon with id {0} was invalid. Attempting to self clean...'.format(rapid_connect_rivid))
                        original_area = catchment_polygon.area
                        catchment_polygon = catchment_polygon.buffer(0)
                        area_ratio = original_area/catchment_polygon.area
                        print('AREA_RATIO', area_ratio)
                        msg_level = "INFO"
                        if round(area_ratio, 5) != 1:
                            msg_level = "WARNING"
                        print('{0}: The cleaned catchment polygon area differs from the'
                              ' original area by {1}%.'.format(msg_level, abs(area_ratio - 1)))
                        intersect_poly = catchment_polygon.intersection(lsm_grid_polygon)
                    if not area_id:
                        #attempt to calculate AREA
                        poly_area = get_poly_area_geo(intersect_poly)
                    else:
                        poly_area = float(get_catchment_feature.GetField(area_id))*intersect_poly.area/catchment_polygon.area

                    index_lsm_grid_lat = lsm_grid_feature_list[sub_lsm_grid_pos]['lat_index']
                    index_lsm_grid_lon = lsm_grid_feature_list[sub_lsm_grid_pos]['lon_index']
                    intersect_grid_info_list.append({'rivid': rapid_connect_rivid,
                                                     'area': poly_area,
                                                     'lsm_grid_lat': lsm_grid_feature_list[sub_lsm_grid_pos]['lat'],
                                                     'lsm_grid_lon': lsm_grid_feature_list[sub_lsm_grid_pos]['lon'],
                                                     'index_lsm_grid_lon': index_lsm_grid_lon,
                                                     'index_lsm_grid_lat': index_lsm_grid_lat})

            npoints = len(intersect_grid_info_list)
            #If no intersection found, add dummy row
            if(npoints <=0):
                connectwriter.writerow([rapid_connect_rivid] + dummy_row_end)
                
            for intersect_grid_info in intersect_grid_info_list:
                connectwriter.writerow([intersect_grid_info['rivid'],
                                        intersect_grid_info['area'],
                                        intersect_grid_info['index_lsm_grid_lon'],
                                        intersect_grid_info['index_lsm_grid_lat'],
                                        npoints,
                                        intersect_grid_info['lsm_grid_lon'],
                                        intersect_grid_info['lsm_grid_lat']])
            
    time_end_all = datetime.utcnow()                                        
    print(time_end_all - time_end_lsm_grid_rtree)
    print("TOTAL TIME: {0}".format(time_end_all - time_start_all))  
    

def CreateWeightTableECMWF(in_ecmwf_nc, 
                           in_catchment_shapefile, 
                           river_id,
                           in_connectivity_file, 
                           out_weight_table,
                           area_id=None, 
                           file_geodatabase=None,
                           ):
                                      
    """
    Create Weight Table for ECMWF Grids
    
    .. note:: The grids are in the RAPIDpy package under the gis/lsm_grids folder.

    Args:
        in_ecmwf_nc(str): Path to the ECMWF NetCDF grid.
        in_catchment_shapefile(str): Path to the Catchment shapefile.
        river_id(str): The name of the field with the river ID (Ex. 'DrainLnID' or 'LINKNO').
        in_connectivity_file(str): The path to the RAPID connectivity file.
        out_weight_table(str): The path to the output weight table file.
        area_id(Optional[str]): The name of the field with the area of each catchment stored in meters squared. Default is it calculate the area.
        file_geodatabase(Optional[str]): Path to the file geodatabase. If you use this option, in_drainage_line is the name of the stream network feature class. (WARNING: Not always stable with GDAL.)
    
    Example:
    
    .. code:: python
    
        from RAPIDpy.gis.weight import CreateWeightTableECMWF

        CreateWeightTableECMWF(in_ecmwf_nc='/path/to/runoff_ecmwf_grid.nc'
                               in_catchment_shapefile='/path/to/catchment.shp',
                               river_id='LINKNO',
                               in_connectivity_file='/path/to/rapid_connect.csv',
                               out_weight_table='/path/to/ecmwf_weight.csv',
                               )    
    """
    #extract ECMWF GRID
    data_ecmwf_nc = Dataset(in_ecmwf_nc)
    variables_list = data_ecmwf_nc.variables.keys()
    in_ecmwf_lat_var = 'lat'
    if 'latitude' in variables_list:
        in_ecmwf_lat_var = 'latitude'
    in_ecmwf_lon_var = 'lon'
    if 'longitude' in variables_list:
        in_ecmwf_lon_var = 'longitude'

    ecmwf_lon = (data_ecmwf_nc.variables[in_ecmwf_lon_var][:] + 180) % 360 - 180 # convert [0, 360] to [-180, 180]
    ecmwf_lat = data_ecmwf_nc.variables[in_ecmwf_lat_var][:] #assume [-90,90]
    data_ecmwf_nc.close()
    
    RTreeCreateWeightTable(ecmwf_lat, ecmwf_lon, 
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
    Create Weight Table for NLDAS, GLDAS grids as well as for 2D Joules, or LIS Grids

    Args:
        in_ldas_nc(str): Path to the land surface model NetCDF grid.
        in_nc_lon_var(str): The variable name in the NetCDF file for the longitude.
        in_nc_lat_var(str): The variable name in the NetCDF file for the latitude.
        in_catchment_shapefile(str): Path to the Catchment shapefile.
        river_id(str): The name of the field with the river ID (Ex. 'DrainLnID' or 'LINKNO').
        in_connectivity_file(str): The path to the RAPID connectivity file.
        out_weight_table(str): The path to the output weight table file.
        area_id(Optional[str]): The name of the field with the area of each catchment stored in meters squared. Default is it calculate the area.
        file_geodatabase(Optional[str]): Path to the file geodatabase. If you use this option, in_drainage_line is the name of the stream network feature class. (WARNING: Not always stable with GDAL.)
    
    Example:
    
    .. code:: python
    
        from RAPIDpy.gis.weight import CreateWeightTableLDAS

        CreateWeightTableLDAS(in_ldas_nc='/path/to/runoff_grid.nc',
                              in_nc_lon_var="lon_110",
                              in_nc_lat_var="lat_110",
                              in_catchment_shapefile='/path/to/catchment.shp',
                              river_id='LINKNO',
                              in_connectivity_file='/path/to/rapid_connect.csv',
                              out_weight_table='/path/to/ldas_weight.csv',
                              )    
    """
    #extract ECMWF GRID
    data_ldas_nc = Dataset(in_ldas_nc)
    variables_list = data_ldas_nc.variables.keys()
    if in_nc_lon_var not in variables_list:
        raise Exception("Invalid longitude variable. Choose from: {0}".format(variables_list))
    if in_nc_lat_var not in variables_list:
        raise Exception("Invalid latitude variable. Choose from: {0}".format(variables_list))
    ldas_lon = data_ldas_nc.variables[in_nc_lon_var][:] #assume [-180, 180]
    ldas_lat = data_ldas_nc.variables[in_nc_lat_var][:] #assume [-90,90]
    data_ldas_nc.close()
    
    RTreeCreateWeightTable(ldas_lat, ldas_lon, 
                           in_catchment_shapefile, river_id,
                           in_connectivity_file, out_weight_table, 
                           file_geodatabase, area_id)


def CreateWeightTable(lsm_grid,
                      lat_var,
                      lon_var,
                      time_var,
                      lat_dim,
                      lon_dim,
                      time_dim,
                      in_catchment_shapefile,
                      river_id,
                      in_rapid_connect,
                      out_weight_table,
                      area_id=None,
                      file_geodatabase=None
                      ):

    with pa.open_mfdataset(lsm_grid,
                           lat_var,
                           lon_var,
                           time_var,
                           lat_dim,
                           lon_dim,
                           time_dim) as xds:
        proj = xds.lsm.projection
        lsm_grid_lat, lsm_grid_lon = xds.lsm.latlon

    RTreeCreateWeightTable(lsm_grid_lat,
                           lsm_grid_lon,
                           in_catchment_shapefile,
                           river_id,
                           in_rapid_connect,
                           out_weight_table,
                           file_geodatabase=None,
                           area_id=None,
                           grid_proj=proj,
                           no_data_value=-9999,
                           )