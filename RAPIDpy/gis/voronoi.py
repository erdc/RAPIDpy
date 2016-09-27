# -*- coding: utf-8 -*-
##
##  voronoi.py
##  RAPIDpy
##
##  Created by Alan D Snow.
##
##  Copyright © 2016 Alan D Snow. All rights reserved.
##  License: BSD 3-Clause

import numpy as np
import os
try:
    from osgeo import ogr, osr
    from scipy.spatial import Voronoi
    from shapely.geometry import Polygon
except Exception:
    raise Exception("You need scipy, gdal, and shapely python packages to run these tools ...")

def _get_voronoi_centroid_array(lsm_lat_array, lsm_lon_array, extent):
    """
    This function generates a voronoi centroid point
    list from arrays of latitude and longitude
    """
    if extent:
        YMin = extent[2]
        YMax = extent[3]
        XMin = extent[0]
        XMax = extent[1]
    else:
        lsm_lat_list = lsm_lat_array
        lsm_lon_list = lsm_lon_array

    ptList = []

    if lsm_lat_array.ndim == 2 and lsm_lon_array.ndim == 2:
        # generate point list with 2D lat lon lists
        if extent:
            #exctract subset within extent
            lsm_dx = np.max(np.absolute(np.diff(lsm_lon_array)))
            lsm_dy = np.max(np.absolute(np.diff(lsm_lat_array, axis=0)))
            
            lsm_lat_indices_from_lat, lsm_lon_indices_from_lat = np.where((lsm_lat_array >= (YMin - 2*lsm_dy)) & (lsm_lat_array <= (YMax + 2*lsm_dy)))
            lsm_lat_indices_from_lon, lsm_lon_indices_from_lon = np.where((lsm_lon_array >= (XMin - 2*lsm_dx)) & (lsm_lon_array <= (XMax + 2*lsm_dx)))

            lsm_lat_indices = np.intersect1d(lsm_lat_indices_from_lat, lsm_lat_indices_from_lon)
            lsm_lon_indices = np.intersect1d(lsm_lon_indices_from_lat, lsm_lon_indices_from_lon)

            lsm_lat_list = lsm_lat_array[lsm_lat_indices,:][:,lsm_lon_indices]
            lsm_lon_list = lsm_lon_array[lsm_lat_indices,:][:,lsm_lon_indices]

        # Create a list of geographic coordinate pairs
        for i in range(len(lsm_lat_indices)):
            for j in range(len(lsm_lon_indices)):
                ptList.append([lsm_lon_list[i][j], lsm_lat_list[i][j]])
    if lsm_lat_array.ndim == 1 and lsm_lon_array.ndim == 1:
        #generate point list with 1D lat lon lists
        if extent:
            Ybuffer = 2 * abs(lsm_lat_array[0]-lsm_lat_array[1])
            Xbuffer = 2 * abs(lsm_lon_array[0]-lsm_lon_array[1])
            # Extract the lat and lon within buffered extent (buffer with 2* interval degree)
            lsm_lat_list = lsm_lat_array[(lsm_lat_array >= (YMin - Ybuffer)) & (lsm_lat_array <= (YMax + Ybuffer))]
            lsm_lon_list = lsm_lon_array[(lsm_lon_array >= (XMin - Xbuffer)) & (lsm_lon_array <= (XMax + Xbuffer))]
            
        # Create a list of geographic coordinate pairs
        for ptX in lsm_lon_list:
            for ptY in lsm_lat_list:
                ptList.append([ptX, ptY])
    else:
        raise IndexError("Lat/Lon lists have invalid dimensions. Only 1D or 2D arrays allowed ...")

    return np.array(ptList) # set-up for input to Delaunay    

def pointsToVoronoiGridShapefile(lat, lon, vor_shp_path, extent=None):
    """
    Converts points to grid via voronoi
    """
    voronoi_centroids = _get_voronoi_centroid_array(lat, lon, extent)

    # set-up output polygon shp
    print("Creating output polygon shp {0}".format(os.path.basename(vor_shp_path)))
    if os.path.exists(vor_shp_path): os.remove(vor_shp_path)
    drv = ogr.GetDriverByName('ESRI Shapefile')
    outShp = drv.CreateDataSource(vor_shp_path)
    osr_geographic_proj = osr.SpatialReference()
    osr_geographic_proj.ImportFromEPSG(4326)
    layer = outShp.CreateLayer('', osr_geographic_proj, ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn('GRID_LAT', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('GRID_LON', ogr.OFTReal))
    layerDefn = layer.GetLayerDefn()
    
    # find nodes surrounding polygon centroid
    # sort nodes in counterclockwise order
    # create polygon perimeter through nodes
    print("Building Voronoi polygons...")
    #compute voronoi
    voronoi_manager  = Voronoi(voronoi_centroids)
    voronoi_verticies = voronoi_manager.vertices
    voronoi_regions = voronoi_manager.regions
    point_id = 0
    for point_index in voronoi_manager.point_region:
        vert_index_list = np.array(voronoi_regions[point_index])
        voronoi_poly_points = []
        if -1 not in vert_index_list and len(vert_index_list) > 3:
            voronoi_poly_points = voronoi_verticies[vert_index_list]
        elif vert_index_list.size>0:
            #ASSUME RECTANGLE
            vert_index_list = vert_index_list[vert_index_list>=0]
            voronoi_poly_points = voronoi_verticies[vert_index_list]
            #CASE 1: 2 valid voronoi vertices
            if vert_index_list.size == 2:
                center_lon = voronoi_centroids[point_id][0]
                center_lat = voronoi_centroids[point_id][1]
                corner_lon1 = voronoi_poly_points[0][0]
                corner_lat1 = voronoi_poly_points[0][1]
                corner_lon2 = voronoi_poly_points[1][0]
                corner_lat2 = voronoi_poly_points[1][1]
                
                #check if need to add points in lon or lat
                if abs(corner_lon1-corner_lon2) > abs(corner_lat1-corner_lat2):
                    dLat = center_lat - corner_lat1
                    #append the corners in order
                    voronoi_poly_points = np.array([[corner_lon1, corner_lat1],
                                                    [corner_lon2, corner_lat2],
                                                    [corner_lon2, center_lat + dLat],
                                                    [corner_lon1, center_lat + dLat]])
                else:
                    dLon = center_lon - corner_lon1
                    #append the corners in order
                    voronoi_poly_points = np.array([[corner_lon1, corner_lat1],
                                                    [corner_lon2, corner_lat2],
                                                    [center_lon+dLon, corner_lat2],
                                                    [center_lon+dLon, corner_lat1]])
            #CADE 2: 1 valid voronoi vertex
            elif vert_index_list.size == 1:
                center_lon = voronoi_centroids[point_id][0]
                center_lat = voronoi_centroids[point_id][1]
                corner_lon = voronoi_poly_points[0][0]
                corner_lat = voronoi_poly_points[0][1]
                dLat = center_lat - corner_lat
                dLon = center_lon - corner_lon
                #append the corners in order
                voronoi_poly_points = np.array([[corner_lon, corner_lat],
                                                [center_lon + dLon, corner_lat],
                                                [center_lon + dLon, center_lat + dLat],
                                                [corner_lon, center_lat + dLat]])

        if len(voronoi_poly_points) == 4:
            poly = ogr.Geometry(ogr.wkbPolygon)
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for i, node in enumerate(voronoi_poly_points):
                if i==0:
                    loopLon = node[0] # grab first node to close ring
                    loopLat = node[1]
                ring.AddPoint(node[0],node[1])
    
            ring.AddPoint(loopLon,loopLat)
            poly.AddGeometry(ring)
            feat = ogr.Feature(layerDefn)
            feat.SetField('GRID_LON', voronoi_centroids[point_id][0])
            feat.SetField('GRID_LAT', voronoi_centroids[point_id][1])
            feat.SetGeometry(poly)  
            layer.CreateFeature(feat)
            feat = poly = ring = None
                
        point_id+=1

def pointsToVoronoiGridArray(lat, lon, extent=None):
    """
    Converts points to grid via voronoi
    """
    voronoi_centroids = _get_voronoi_centroid_array(lat, lon, extent)
    
    # find nodes surrounding polygon centroid
    # sort nodes in counterclockwise order
    # create polygon perimeter through nodes
    print("Building Voronoi polygons...")
    #compute voronoi
    voronoi_manager  = Voronoi(voronoi_centroids)
    voronoi_verticies = voronoi_manager.vertices
    voronoi_regions = voronoi_manager.regions
    point_id = 0
    feature_list = []
    for point_index in voronoi_manager.point_region:
        vert_index_list = np.array(voronoi_regions[point_index])
        voronoi_poly_points = []
        if -1 not in vert_index_list and len(vert_index_list) > 3:
            voronoi_poly_points = voronoi_verticies[vert_index_list]
        elif vert_index_list.size>0:
            vert_index_list = vert_index_list[vert_index_list>=0]
            voronoi_poly_points = voronoi_verticies[vert_index_list]
            #CASE 1: 2 valid voronoi vertices
            if vert_index_list.size == 2:
                center_lon = voronoi_centroids[point_id][0]
                center_lat = voronoi_centroids[point_id][1]
                corner_lon1 = voronoi_poly_points[0][0]
                corner_lat1 = voronoi_poly_points[0][1]
                corner_lon2 = voronoi_poly_points[1][0]
                corner_lat2 = voronoi_poly_points[1][1]
                
                #check if need to add points in lon or lat
                if abs(corner_lon1-corner_lon2) > abs(corner_lat1-corner_lat2):
                    dLat = center_lat - corner_lat1
                    #append the corners in order
                    voronoi_poly_points = np.array([[corner_lon1, corner_lat1],
                                                    [corner_lon2, corner_lat2],
                                                    [corner_lon2, center_lat + dLat],
                                                    [corner_lon1, center_lat + dLat]])
                else:
                    dLon = center_lon - corner_lon1
                    #append the corners in order
                    voronoi_poly_points = np.array([[corner_lon1, corner_lat1],
                                                    [corner_lon2, corner_lat2],
                                                    [center_lon+dLon, corner_lat2],
                                                    [center_lon+dLon, corner_lat1]])
            #CASE 2: 1 valid voronoi vertex
            elif vert_index_list.size == 1:
                center_lon = voronoi_centroids[point_id][0]
                center_lat = voronoi_centroids[point_id][1]
                corner_lon = voronoi_poly_points[0][0]
                corner_lat = voronoi_poly_points[0][1]
                dLat = center_lat - corner_lat
                dLon = center_lon - corner_lon
                #append the corners in order
                voronoi_poly_points = np.array([[corner_lon, corner_lat],
                                                [center_lon + dLon, corner_lat],
                                                [center_lon + dLon, center_lat + dLat],
                                                [corner_lon, center_lat + dLat]])
                
        if len(voronoi_poly_points) == 4:
            feature_list.append({'polygon': Polygon(voronoi_poly_points),
                                 'lon' : voronoi_centroids[point_id][0],
                                 'lat' : voronoi_centroids[point_id][1]})
                
        point_id+=1
    return feature_list
