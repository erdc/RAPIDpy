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
from json import loads as json_loads
from netCDF4 import Dataset
import numpy as np
import os
from subprocess import PIPE, Popen

try:
    from osgeo import gdal, ogr, osr
    from pyproj import Proj, transform
    from shapely.geometry import MultiPolygon, Polygon, shape
    import rtree #http://toblerity.org/rtree/install.html
except Exception:
    raise Exception("You need the gdal, pyproj, shapely, and rtree python package to run these tools ...")

#local
from voronoi import pointsToVoronoiGridArray, pointsToVoronoiGridShapefile

gdal.UseExceptions()

def get_poly_area_geo(poly):
    """
    Calculates the area in meters squared of the individual polygon
    """
    lon_list, lat_list = poly.exterior.coords.xy
    reprojected_for_area = Proj("+proj=aea +lat_1={0} +lat_1={1} +lat_0={2} +lon_0={3}".format(np.min(lat_list),
                                                                                               np.max(lat_list),
                                                                                               np.mean(lat_list),
                                                                                               np.mean(lon_list)))
    #reproject polygon to get area
    lon_list, lat_list = reprojected_for_area(lon_list, lat_list)
    poly = Polygon(zip(lon_list, lat_list))
    return poly.area

def find_nearest(array, value):
    """
    Get the nearest index to value searching for
    """
    return (np.abs(array-value)).argmin()
    
def RTreeCreateWeightTable(lsm_grid_lat, lsm_grid_lon, 
                           in_catchment_shapefile, river_id,
                           in_rapid_connect, out_weight_table,
                           file_geodatabase=None, area_id=None):
                                      
    """
    Create Weight Table for Land Surface Model Grids
    """
    time_start_all = datetime.utcnow()    
    
    print "Generating LSM Grid Thiessen Array ..."
    if file_geodatabase:
        gdb_driver = ogr.GetDriverByName("OpenFileGDB")
        ogr_file_geodatabase = gdb_driver.Open(file_geodatabase, 0)
        ogr_catchment_shapefile_lyr = ogr_file_geodatabase.GetLayer(in_catchment_shapefile)
    else:
        ogr_catchment_shapefile = ogr.Open(in_catchment_shapefile)
        ogr_catchment_shapefile_lyr = ogr_catchment_shapefile.GetLayer()
        
    ogr_catchment_shapefile_lyr_proj = ogr_catchment_shapefile_lyr.GetSpatialRef()
    original_catchment_proj = Proj(ogr_catchment_shapefile_lyr_proj.ExportToProj4())
    geographic_proj =  Proj(init='EPSG:4326') #geographic
    extent = ogr_catchment_shapefile_lyr.GetExtent()
    if original_catchment_proj != geographic_proj:
        x, y = transform(original_catchment_proj,
                         geographic_proj,
                         [extent[0], extent[1]], 
                         [extent[2], extent[3]])
        extent = [min(x), max(x), min(y), max(y)]
            
    lsm_grid_feature_list = pointsToVoronoiGridArray(lsm_grid_lat, lsm_grid_lon, extent)
    time_end_lsm_grid_thiessen = datetime.utcnow()
    print time_end_lsm_grid_thiessen - time_start_all
    
    print "Generating LSM Grid Rtree ..."
    rtree_idx = rtree.index.Index()
    # Populate R-tree index with bounds of ECMWF grid cells
    for lsm_grid_pos, lsm_grid_feature in enumerate(lsm_grid_feature_list):
        rtree_idx.insert(lsm_grid_pos, lsm_grid_feature['polygon'].bounds)
     
    time_end_lsm_grid_rtree = datetime.utcnow()
    print time_end_lsm_grid_rtree - time_end_lsm_grid_thiessen
    
    print "Retrieving catchment river id list ..."
    catchment_rivid_list = []
    for catchment_feature in ogr_catchment_shapefile_lyr:
        catchment_rivid_list.append(catchment_feature.GetField(river_id))
     
    catchment_rivid_list = np.array(catchment_rivid_list, dtype=np.int32)
        
    print "Reading in RAPID connect file ..."
    rapid_connect_rivid_list = []
    with open(in_rapid_connect, "rb") as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            rapid_connect_rivid_list.append(row[0])
    
    rapid_connect_rivid_list = np.array(rapid_connect_rivid_list, dtype=np.int32)

    print "Find LSM grid cells that intersect with each catchment"
    print "and write out weight table ..."
    
    
    dummy_row_end = [0,
                    np.where(lsm_grid_lon == lsm_grid_feature_list[0]['lon'])[0][0],
                    np.where(lsm_grid_lat == lsm_grid_feature_list[0]['lat'])[0][0],
                    1,
                    lsm_grid_feature_list[0]['lon'],
                    lsm_grid_feature_list[0]['lat']
                    ]
                    
    with open(out_weight_table, 'wb') as csvfile:
        connectwriter = csv.writer(csvfile)
        connectwriter.writerow(['rivid', 'area', 'index_lsm_grid_lon', 'index_lsm_grid_lat', 
                                'npoints', 'lsm_grid_lon', 'lsm_grid_lat'])
        geographic_proj =  Proj(init='EPSG:4326') #geographic
        osr_geographic_proj = osr.SpatialReference()
        osr_geographic_proj.ImportFromEPSG(4326)
        proj_transform = None
        if original_catchment_proj != geographic_proj:
            proj_transform = osr.CoordinateTransformation(ogr_catchment_shapefile_lyr_proj, osr_geographic_proj)
            
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
            feat_geom_count = feat_geom.GetGeometryCount()
            #make sure coordinates are geographic
            if proj_transform:
                feat_geom.Transform(proj_transform)
            main_polygon_arr = []
            for ring_idx in xrange(feat_geom_count):
                feat_ring = feat_geom.GetGeometryRef(ring_idx)
                num_ring_points = feat_ring.GetPointCount()
                
                if num_ring_points <= 0:
                    feature_info = json_loads(feat_ring.ExportToJson())
                    main_polygon_arr.append(Polygon(shape(feature_info)))
                else:
                    poly_pt_arr = []
                    for pt_idx in xrange(num_ring_points):
                        x, y, z = feat_ring.GetPoint(pt_idx)
                        poly_pt_arr.append((x,y))
                    poly_pt_arr.append((poly_pt_arr[0])) #close the loop
                    main_polygon_arr.append(Polygon(poly_pt_arr))
                    
            if len(main_polygon_arr) == 1:
                catchment_polygon = main_polygon_arr[0]
            else:
                catchment_polygon = MultiPolygon(main_polygon_arr)
                
            for sub_lsm_grid_pos in rtree_idx.intersection(catchment_polygon.bounds):
                if sub_lsm_grid_pos != catchment_pos:
                    if catchment_polygon.intersects(lsm_grid_feature_list[sub_lsm_grid_pos]['polygon']):
                        intersect_poly = catchment_polygon.intersection(lsm_grid_feature_list[sub_lsm_grid_pos]['polygon'])
                        if not area_id:
                            #attempt to calculate AREA
                            try:
                                poly_area = get_poly_area_geo(intersect_poly)
                            except AttributeError:
                                #multipolygon
                                poly_area = 0
                                for subpoly in intersect_poly:
                                    poly_area += get_poly_area_geo(subpoly)
                        else:
                            poly_area = float(catchment_polygon.GetFeature(area_id))*intersect_poly.area/catchment_polygon.area
                            
                        index_lsm_grid_lon = np.where(lsm_grid_lon == lsm_grid_feature_list[sub_lsm_grid_pos]['lon'])[0][0]
                        index_lsm_grid_lat = np.where(lsm_grid_lat == lsm_grid_feature_list[sub_lsm_grid_pos]['lat'])[0][0]            
                        intersect_grid_info_list.append({'rivid' : rapid_connect_rivid,
                                                         'area' : poly_area,
                                                         'lsm_grid_lat': lsm_grid_feature_list[sub_lsm_grid_pos]['lat'],
                                                         'lsm_grid_lon': lsm_grid_feature_list[sub_lsm_grid_pos]['lon'],
                                                         'index_lsm_grid_lon': index_lsm_grid_lon,
                                                         'index_lsm_grid_lat': index_lsm_grid_lat})
            npoints = len(intersect_grid_info_list)
            for intersect_grid_info in intersect_grid_info_list:
                connectwriter.writerow([intersect_grid_info['rivid'],
                                        intersect_grid_info['area'],
                                        intersect_grid_info['index_lsm_grid_lon'],
                                        intersect_grid_info['index_lsm_grid_lat'],
                                        npoints,
                                        intersect_grid_info['lsm_grid_lon'],
                                        intersect_grid_info['lsm_grid_lat']])
            
    time_end_all = datetime.utcnow()                                        
    print time_end_all - time_end_lsm_grid_rtree
    print "TOTAL TIME:",   time_end_all - time_start_all  
    
def GDALCreateWeightTable(lsm_grid_lat, lsm_grid_lon,
                          in_catchment_shapefile, river_id,
                          in_rapid_connect, out_weight_table):
                              
    time_start_all = datetime.utcnow()    
    
    print "Generating LSM Thiessen Shapefile ..."
    time_start_lsm_grid_thiessen = datetime.utcnow()

    shapefile_working_directory = os.path.dirname(in_catchment_shapefile)

    ogr_catchment_shapefile = ogr.Open(in_catchment_shapefile)
    ogr_catchment_shapefile_lyr = ogr_catchment_shapefile.GetLayer()
    ogr_catchment_shapefile_lyr_proj = ogr_catchment_shapefile_lyr.GetSpatialRef()
    osr_geographic_proj = osr.SpatialReference()
    osr_geographic_proj.ImportFromEPSG(4326)
    #reproject if needed
    if ogr_catchment_shapefile_lyr_proj != osr_geographic_proj:
        print "Reprojecting Catchment to EPSG:4326 ..."
        reprojected_shapefile_name = "{0}_epsg4326.shp".format(os.path.splitext(in_catchment_shapefile)[0])
        process = Popen(["ogr2ogr", "-f", "ESRI Shapefile",
                         reprojected_shapefile_name,
                         in_catchment_shapefile,
                         "-t_srs", "EPSG:4326", 
                         "-s_srs", ogr_catchment_shapefile_lyr_proj.ExportToProj4(),
                         ],
                         stdout=PIPE, stderr=PIPE, shell=False)
        out, err = process.communicate()
        if err:
            print "WARNINGS Reprojecting ..."
                
        ogr_catchment_shapefile = ogr_catchment_shapefile_lyr = None
        ogr_catchment_shapefile_lyr_proj = None
        in_catchment_shapefile = reprojected_shapefile_name
        
    print "Generating LSM thiesen polygon ..."
    ogr_catchment_shapefile = ogr.Open(in_catchment_shapefile)
    ogr_catchment_shapefile_lyr = ogr_catchment_shapefile.GetLayer()
    extent = ogr_catchment_shapefile_lyr.GetExtent()
    ogr_catchment_shapefile = ogr_catchment_shapefile_lyr = None
            
    lsm_grid_thiessen_shapefile = os.path.join(shapefile_working_directory, "lsm_grid_thiesen_poly.shp")
    pointsToVoronoiGridShapefile(lsm_grid_lat, lsm_grid_lon, lsm_grid_thiessen_shapefile, extent)
    time_end_lsm_grid_thiessen = datetime.utcnow()
    print time_end_lsm_grid_thiessen - time_start_lsm_grid_thiessen
    
    print "Intersecting LSM thiesen shapefile with catchment ..."
    ogr_ds = ogr.Open(shapefile_working_directory, True)
    SQL = """\
          SELECT ST_Intersection(A.geometry, B.geometry) AS geometry, A.*, B.*
          FROM {0} A, {1} B
          WHERE ST_Intersects(A.geometry, B.geometry)
          """.format(os.path.basename(os.path.splitext(in_catchment_shapefile)[0]),
                     os.path.basename(os.path.splitext(lsm_grid_thiessen_shapefile)[0]))
    intersected_layer = ogr_ds.ExecuteSQL(SQL, dialect='SQLITE')
    
    time_end_intersect = datetime.utcnow()
    print time_end_intersect - time_end_lsm_grid_thiessen
    print "Retrieving intersected layer river id list ..."
    intersect_rivid_list = []
    for intersect_feature in intersected_layer:
        intersect_rivid_list.append(intersect_feature.GetField(river_id))
     
    intersect_rivid_list = np.array(intersect_rivid_list, dtype=np.int32)
        
    print "Reading in RAPID connect file ..."
    rapid_connect_rivid_list = []
    with open(in_rapid_connect, "rb") as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            rapid_connect_rivid_list.append(row[0])
    
    rapid_connect_rivid_list = np.array(rapid_connect_rivid_list, dtype=np.int32)

    print "Find ECMWF grid cells that intersect with each catchment"
    print "and write out weight table ..."
    dummy_intersect_feature = intersected_layer.GetFeature(0)
    dummy_intersect_lsm_grid_lon = dummy_intersect_feature.GetField('GRID_LON')
    dummy_intersect_lsm_grid_lat = dummy_intersect_feature.GetField('GRID_LAT')
    dummy_intersect_feature = None
    
    dummy_row_end = [0,
                    np.where(lsm_grid_lon == dummy_intersect_lsm_grid_lon)[0][0],
                    np.where(lsm_grid_lat == dummy_intersect_lsm_grid_lat)[0][0],
                    1,
                    dummy_intersect_lsm_grid_lon,
                    dummy_intersect_lsm_grid_lat
                    ]
                    
    with open(out_weight_table, 'wb') as csvfile:
        connectwriter = csv.writer(csvfile)
        connectwriter.writerow(['rivid', 'area', 'index_lsm_grid_lon', 'index_lsm_grid_lat', 
                                'npoints', 'lsm_grid_lon', 'lsm_grid_lat'])
            
        for rapid_connect_rivid in rapid_connect_rivid_list:
            intersect_positions = np.where(intersect_rivid_list==rapid_connect_rivid)[0]
            num_ind_points = len(intersect_positions)
            if num_ind_points <= 0:
                # if point not in array, append dummy data for one point of data
                # streamID, area_sqm, lon_index, lat_index, npoints
                connectwriter.writerow([rapid_connect_rivid] + dummy_row_end)
                continue
            
            for intersect_index in intersect_positions:
                intersect_feature = intersected_layer.GetFeature(intersect_index)
                
                feat_geom = intersect_feature.GetGeometryRef()
                feat_geom_count = feat_geom.GetGeometryCount()
                main_polygon_arr = []
                for ring_idx in xrange(feat_geom_count):
                    feat_ring = feat_geom.GetGeometryRef(ring_idx)
                    num_ring_points = feat_ring.GetPointCount()
                    
                    if num_ring_points <= 0:
                        feature_info = json_loads(feat_ring.ExportToJson())
                        main_polygon_arr.append(Polygon(shape(feature_info)))
                    else:
                        poly_pt_arr = []
                        for pt_idx in xrange(num_ring_points):
                            x, y, z = feat_ring.GetPoint(pt_idx)
                            poly_pt_arr.append((x,y))
                        poly_pt_arr.append((poly_pt_arr[0])) #close the loop
                        main_polygon_arr.append(Polygon(poly_pt_arr))
                if len(main_polygon_arr) == 1:
                    intersect_polygon = main_polygon_arr[0]
                else:
                    intersect_polygon = MultiPolygon(main_polygon_arr)

                #attempt to calculate AREA
                try:
                    poly_area = get_poly_area_geo(intersect_polygon)
                except AttributeError:
                    #multipolygon
                    poly_area = 0
                    for subpoly in intersect_polygon:
                        poly_area += get_poly_area_geo(subpoly)
                    
                intersect_lsm_grid_lon = intersect_feature.GetField('GRID_LON')
                intersect_lsm_grid_lat = intersect_feature.GetField('GRID_LAT')

                index_lsm_grid_lon = np.where(lsm_grid_lon == intersect_lsm_grid_lon)[0][0]
                index_lsm_grid_lat = np.where(lsm_grid_lat == intersect_lsm_grid_lat)[0][0]
                
                connectwriter.writerow([rapid_connect_rivid,
                                        poly_area,
                                        index_lsm_grid_lon,
                                        index_lsm_grid_lat,
                                        num_ind_points,
                                        intersect_lsm_grid_lon,
                                        intersect_lsm_grid_lat])
            

    time_end_all = datetime.utcnow()                                        
    print time_end_all - time_end_intersect

    print "TOTAL TIME:",   time_end_all - time_start_all
           
def CreateWeightTableECMWF(in_ecmwf_nc, 
                           in_catchment_shapefile, 
                           river_id,
                           in_rapid_connect, 
                           out_weight_table,
                           file_geodatabase=None,
                           area_id=None, 
                           method="rtree"):
                                      
    """
    Create Weight Table for ECMWF Grids
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
    
    if method.lower() == "rtree":
        RTreeCreateWeightTable(ecmwf_lat, ecmwf_lon, 
                               in_catchment_shapefile, river_id,
                               in_rapid_connect, out_weight_table, 
                               file_geodatabase, area_id)
    elif method.lower() == "gdal" and not file_geodatabase:
        GDALCreateWeightTable(ecmwf_lat, ecmwf_lon, 
                              in_catchment_shapefile, river_id,
                              in_rapid_connect, out_weight_table)
    else:
        raise Exception("ERROR: Invalid run method. Valid run methods are rtree and gdal (no File Geodatabase support for GDAL method).")

def CreateWeightTableLDAS(in_ldas_nc,
                          in_nc_lon_var,
                          in_nc_lat_var,
                          in_catchment_shapefile, 
                          river_id,
                          in_rapid_connect, 
                          out_weight_table,
                          file_geodatabase=None,
                          area_id=None, 
                          method="rtree"):
                                      
    """
    Create Weight Table for LDAS Grids
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
    
    if method.lower() == "rtree":
        RTreeCreateWeightTable(ldas_lat, ldas_lon, 
                               in_catchment_shapefile, river_id,
                               in_rapid_connect, out_weight_table, 
                               file_geodatabase, area_id)
    elif method.lower() == "gdal" and not file_geodatabase:
        GDALCreateWeightTable(ldas_lat, ldas_lon, 
                              in_catchment_shapefile, river_id,
                              in_rapid_connect, out_weight_table)
    else:
        raise Exception("ERROR: Invalid run method. Valid run methods are rtree and gdal (no File Geodatabase support for GDAL method).")

def CreateWeightTableLIS(in_lis_nc,
                         in_nc_lon_var,
                         in_nc_lat_var,
                         in_nc_lon_dim,
                         in_nc_lat_dim,
                         in_catchment_shapefile, 
                         river_id,
                         in_rapid_connect, 
                         out_weight_table,
                         file_geodatabase=None,
                         area_id=None, 
                         method="rtree"):
                                      
    """
    Create Weight Table for LDAS Grids
    """
    data_nc = Dataset(in_lis_nc)

    # Obtain geographic coordinates
    lon = np.unique(np.concatenate(data_nc.variables[in_nc_lon_var][:])) #assume [-180,180]
    lat = np.unique(np.concatenate(data_nc.variables[in_nc_lat_var][:])) #assume [-90,90]
    size_xdim = len(data_nc.dimensions[in_nc_lat_dim])
    #remove missing lon values
    if 'missing_value' in data_nc.variables[in_nc_lon_var].ncattrs():
        lon = lon[lon!=data_nc.variables[in_nc_lon_var].getncattr('missing_value')]
    if size_xdim != len(lon):
        raise Exception("Latitude dimension in data does not match netCDF file dimension")
    size_ydim = len(data_nc.dimensions[in_nc_lon_dim])
    #remove missing lat values
    if 'missing_value' in data_nc.variables[in_nc_lat_var].ncattrs():
        lat = lat[lat!=data_nc.variables[in_nc_lat_var].getncattr('missing_value')]
    if size_ydim != len(lat):
        raise Exception("Latitude dimension in data does not match netCDF file dimension")
    data_nc.close()

    
    if method.lower() == "rtree":
        RTreeCreateWeightTable(lat, lon, 
                               in_catchment_shapefile, river_id,
                               in_rapid_connect, out_weight_table, 
                               file_geodatabase, area_id)
    elif method.lower() == "gdal" and not file_geodatabase:
        GDALCreateWeightTable(lat, lon, 
                              in_catchment_shapefile, river_id,
                              in_rapid_connect, out_weight_table)
    else:
        raise Exception("ERROR: Invalid run method. Valid run methods are rtree and gdal (no File Geodatabase support for GDAL method).")