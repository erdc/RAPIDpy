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

def pointsToVoronoiGridShapefile(lat, lon, vor_shp_path, extent=None):
    """
    Converts points to grid via voronoi
    """
    if extent:
        Ybuffer = 2 * abs(lat[0]-lat[1])
        Xbuffer = 2 * abs(lon[0]-lon[1])
        # Extract the lat and lon within buffered extent (buffer with 2* interval degree)
        YMin = extent[2]
        YMax = extent[3]
        XMin = extent[0]
        XMax = extent[1]
        lat0 = lat[(lat >= (YMin - Ybuffer)) & (lat <= (YMax + Ybuffer))]
        lon0 = lon[(lon >= (XMin - Xbuffer)) & (lon <= (XMax + Xbuffer))]
    else:
        lat0 = lat
        lon0 = lon
    # collect point coordinates and ID
    ptList = []
    for lon_index, ptX in enumerate(lon0):
        for lat_index, ptY in enumerate(lat0):
            ptList.append([ptX, ptY]) 

    voronoi_centroids = np.array(ptList) # set-up for input to Delaunay

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
        elif vert_index_list.any():
            #ASSUME RECTANGLE
            vert_index_list = vert_index_list[vert_index_list>=0]
            num_verts = 0
            if vert_index_list.any():
                num_verts = 1
                if hasattr(vert_index_list, "__len__"):
                    num_verts = len(vert_index_list)
            voronoi_poly_points = voronoi_verticies[vert_index_list]
            #CASE 1: 2 valid voronoi vertices
            if num_verts == 2:
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
                #print "TWO", voronoi_poly_points
            #CADE 2: 1 valid voronoi vertex
            elif num_verts == 1:
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
                #print "ONE", voronoi_poly_points
                
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
    if extent:
        buffer = 2 * max(abs(lat[0]-lat[1]),abs(lon[0] - lon[1]))
        # Extract the lat and lon within buffered extent (buffer with 2* interval degree)
        YMin = extent[2]
        YMax = extent[3]
        XMin = extent[0]
        XMax = extent[1]
        lat0 = lat[(lat >= (YMin - buffer)) & (lat <= (YMax + buffer))]
        lon0 = lon[(lon >= (XMin - buffer)) & (lon <= (XMax + buffer))]
    else:
        lat0 = lat
        lon0 = lon
    # collect point coordinates and ID
    ptList = []
    ptGridInfo = []
    for lon_index, ptX in enumerate(lon0):
        for lat_index, ptY in enumerate(lat0):
            ptList.append([ptX, ptY]) 
            ptGridInfo.append([lon_index, lat_index])
    voronoi_centroids = np.array(ptList) # set-up for input to Delaunay
    
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
        elif vert_index_list.any():
            #ASSUME RECTANGLE
            vert_index_list = vert_index_list[vert_index_list>=0]
            num_verts = 0
            if vert_index_list.any():
                num_verts = 1
                if hasattr(vert_index_list, "__len__"):
                    num_verts = len(vert_index_list)
            voronoi_poly_points = voronoi_verticies[vert_index_list]
            #CASE 1: 2 valid voronoi vertices
            if num_verts == 2:
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
                #print "TWO", voronoi_poly_points
            #CADE 2: 1 valid voronoi vertex
            elif num_verts == 1:
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
                #print "ONE", voronoi_poly_points
                
        if len(voronoi_poly_points) == 4:
            feature_list.append({'polygon': Polygon(voronoi_poly_points),
                                 'lon' : voronoi_centroids[point_id][0],
                                 'lat' : voronoi_centroids[point_id][1]})
                
        point_id+=1
    return feature_list