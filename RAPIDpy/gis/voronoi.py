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
from itertools import product
try:
    from osgeo import ogr, osr
    from scipy.spatial import Voronoi
    from shapely.geometry import Polygon
    from shapely.wkt import loads as shapely_load_wkt
except Exception:
    raise Exception("You need scipy, gdal, and shapely python packages to run these tools ...")


def _get_lat_lon_indices(lsm_lat_array, lsm_lon_array, lat, lon):
    """
    Determines the index in the array (1D or 2D) where the
    lat/lon point is
    """
    if lsm_lat_array.ndim == 2 and lsm_lon_array.ndim == 2:
        lsm_lat_indices_from_lat, lsm_lon_indices_from_lat = np.where((lsm_lat_array == lat))
        lsm_lat_indices_from_lon, lsm_lon_indices_from_lon = np.where((lsm_lon_array == lon))

        index_lsm_grid_lat = np.intersect1d(lsm_lat_indices_from_lat, lsm_lat_indices_from_lon)[0]
        index_lsm_grid_lon = np.intersect1d(lsm_lon_indices_from_lat, lsm_lon_indices_from_lon)[0]

    elif lsm_lat_array.ndim == 1 and lsm_lon_array.ndim == 1:
        index_lsm_grid_lon = np.where(lsm_lon_array == lon)[0][0]
        index_lsm_grid_lat = np.where(lsm_lat_array == lat)[0][0]
    else:
        raise IndexError("Lat/Lon lists have invalid dimensions. Only 1D or 2D arrays allowed ...")

    return index_lsm_grid_lat, index_lsm_grid_lon


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

    ptList = []
    if (lsm_lat_array.ndim == 2) and (lsm_lon_array.ndim == 2):
        # generate point list with 2D lat lon lists
        if extent:
            #exctract subset within extent
            lsm_dx = np.max(np.absolute(np.diff(lsm_lon_array)))
            lsm_dy = np.max(np.absolute(np.diff(lsm_lat_array, axis=0)))
            
            #remove values with NaN
            lsm_lat_array = np.ma.filled(lsm_lat_array, fill_value=-9999)
            lsm_lon_array = np.ma.filled(lsm_lon_array, fill_value=-9999)
            
            lsm_lat_indices_from_lat, lsm_lon_indices_from_lat = np.where((lsm_lat_array >= (YMin - 2*lsm_dy)) & (lsm_lat_array <= (YMax + 2*lsm_dy)))
            lsm_lat_indices_from_lon, lsm_lon_indices_from_lon = np.where((lsm_lon_array >= (XMin - 2*lsm_dx)) & (lsm_lon_array <= (XMax + 2*lsm_dx)))

            lsm_lat_indices = np.intersect1d(lsm_lat_indices_from_lat, lsm_lat_indices_from_lon)
            lsm_lon_indices = np.intersect1d(lsm_lon_indices_from_lat, lsm_lon_indices_from_lon)

            lsm_lat_list = lsm_lat_array[lsm_lat_indices,:][:,lsm_lon_indices]
            lsm_lon_list = lsm_lon_array[lsm_lat_indices,:][:,lsm_lon_indices]

        else:
            lsm_lat_indices = lsm_lat_array
            lsm_lon_indices = lsm_lon_array
            lsm_lat_list = lsm_lat_array
            lsm_lon_list = lsm_lon_array
        # Create a list of geographic coordinate pairs
        for i in range(len(lsm_lat_indices)):
            for j in range(len(lsm_lon_indices)):
                ptList.append([lsm_lon_list[i][j], lsm_lat_list[i][j]])
                
    elif lsm_lat_array.ndim == 1 and lsm_lon_array.ndim == 1:
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

    if len(ptList) <=0:
        raise IndexError("The watershed is outside of the bounds of the land surface model grid ...")

    return np.array(ptList) # set-up for input to Delaunay    

def _get_voronoi_poly_points(vert_index_list, voronoi_vertices, voronoi_centroid):
    """
    This function returns the corner points for a polygon from scipy voronoi information
    """
    voronoi_poly_points = []
    if -1 not in vert_index_list and len(vert_index_list) > 3:
        voronoi_poly_points = voronoi_vertices[vert_index_list]
    elif vert_index_list.size>0:
        #ASSUME RECTANGLE
        vert_index_list = vert_index_list[vert_index_list>=0]
        voronoi_poly_points = voronoi_vertices[vert_index_list]
        #CASE 1: 2 valid voronoi vertices
        if vert_index_list.size == 2:
            center_lon = voronoi_centroid[0]
            center_lat = voronoi_centroid[1]
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
            center_lon = voronoi_centroid[0]
            center_lat = voronoi_centroid[1]
            corner_lon = voronoi_poly_points[0][0]
            corner_lat = voronoi_poly_points[0][1]
            dLat = center_lat - corner_lat
            dLon = center_lon - corner_lon
            #append the corners in order
            voronoi_poly_points = np.array([[corner_lon, corner_lat],
                                            [center_lon + dLon, corner_lat],
                                            [center_lon + dLon, center_lat + dLat],
                                            [corner_lon, center_lat + dLat]])


    return voronoi_poly_points

def pointsToVoronoiGridShapefile(lat, lon, vor_shp_path, extent=None):
    """
    Converts points to shapefile grid via voronoi
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
    voronoi_vertices = voronoi_manager.vertices
    voronoi_regions = voronoi_manager.regions
    for point_id, region_index in enumerate(voronoi_manager.point_region):
        vert_index_list = np.array(voronoi_regions[region_index])
        voronoi_centroid = voronoi_centroids[point_id]
        voronoi_poly_points = _get_voronoi_poly_points(vert_index_list, 
                                                       voronoi_vertices, 
                                                       voronoi_centroid)
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
            feat.SetField('GRID_LON', float(voronoi_centroid[0]))
            feat.SetField('GRID_LAT', float(voronoi_centroid[1]))
            feat.SetGeometry(poly)  
            layer.CreateFeature(feat)
            feat = poly = ring = None


def pointsToVoronoiGridArray(lat, lon, extent=None, coord_transform_in=None, coord_transforma_out=None, no_data_value=None):
    """
    Converts points to grid array via voronoi
    """
    voronoi_centroids = _get_voronoi_centroid_array(lat, lon, extent)

    lon_lat = np.vstack([(lon, lat)])
    number_of_rows = lon_lat.shape[1]
    number_of_columns = lon_lat.shape[2]

    # convert to shapely objects
    voronoi_centroids = []
    for i, j in product(range(number_of_rows), range(number_of_columns)):
        voronoi_centroids.append(list(lon_lat[:, i, j]))

    if coord_transform_in:
        proj_centroids = coord_transform_in.TransformPoints(voronoi_centroids)
        voronoi_centroids = np.asarray(proj_centroids)[:, 0:2]
    
    # find nodes surrounding polygon centroid
    # sort nodes in counterclockwise order
    # create polygon perimeter through nodes
    print("Building Voronoi polygons...")
    #compute voronoi
    voronoi_manager  = Voronoi(voronoi_centroids)
    voronoi_vertices = voronoi_manager.vertices
    voronoi_regions = voronoi_manager.regions
    feature_list = []
    for point_id, region_index in enumerate(voronoi_manager.point_region):
        vert_index_list = np.array(voronoi_regions[region_index])
        voronoi_centroid = voronoi_centroids[point_id]
        voronoi_poly_points = _get_voronoi_poly_points(vert_index_list, 
                                                       voronoi_vertices, 
                                                       voronoi_centroid)
                
        if len(voronoi_poly_points) == 4:
            poly = Polygon(voronoi_poly_points)

            if coord_transforma_out:
                ogr_poly = ogr.CreateGeometryFromWkt(poly.wkt)
                ogr_poly.Transform(coord_transforma_out)
                poly = shapely_load_wkt(ogr_poly.ExportToWkt())

                voronoi_centroid = coord_transforma_out.TransformPoint(*voronoi_centroid)

            lat_index, lon_index = _get_lat_lon_indices(lat, lon,
                                                        voronoi_centroid[1],
                                                        voronoi_centroid[0])

            feature_list.append({'polygon': poly,
                                 'lon': voronoi_centroid[0],
                                 'lat': voronoi_centroid[1],
                                 'lat_index': lat_index,
                                 'lon_index': lon_index,
                                 })

    write_shapefile(feature_list)

    return feature_list


def _fill_nan(a, axis):

    filled = np.copy(a)
    rolled = np.roll(a, shift=1, axis=axis)
    # rolled = np.nanmean(a, axis=axis)
    idx = np.where(np.isnan(filled))
    filled[idx] = np.take(rolled, idx[1])  # TODO why inx[1]? see (http://stackoverflow.com/questions/18689235/numpy-array-replace-nan-values-with-average-of-columns)

    rolled = np.roll(a, shift=-1, axis=axis)
    idx = np.where(np.isnan(filled))
    filled[idx] = np.take(rolled, idx[1])

    return filled


def pointsArrayToPolygonGridList(lat, lon, extent=None, coord_transform_in=None, coord_transform_out=None,
                                no_data_value=-9999):
    """
    Converts points to grid array via voronoi
    """
    if lat.ndim == 1 or lon.ndim == 1:
        number_of_rows = lat.shape[0]
        number_of_columns = lon.shape[0]

        # create a 2D matrix of lats and lons
        lon, lat = np.meshgrid(lon, lat)
        # lat = np.repeat(lat, number_of_columns, axis=0).reshape(number_of_rows, -1)
        # lon = np.repeat(lon, number_of_rows, axis=0).reshape(number_of_columns, -1).T

    # combine lat and lon arrays into 3D array
    lon_lat = np.vstack([(lon, lat)])

    lon_lat[lon_lat == no_data_value] = np.nan

    if coord_transform_in is not None:
        centroids = np.apply_along_axis(lambda p: coord_transform_in.TransformPoint(float(p[0]), float(p[1])), 0, lon_lat)
    else:
        centroids = lon_lat

    v_distances = np.linalg.norm(centroids[:, 0:-1, :] - centroids[:, 1:, :], axis=0) / 2.0
    h_distances = np.linalg.norm(centroids[:, :, 0:-1] - centroids[:, :, 1:], axis=0) / 2.0

    v_distances = _fill_nan(v_distances, axis=1)
    h_distances = _fill_nan(h_distances, axis=0)

    v_up = np.insert(v_distances, 0, v_distances[0, :], axis=0)
    v_down = np.insert(v_distances, -1, v_distances[-1, :], axis=0)

    h_left = np.insert(h_distances, 0, h_distances[:, 0], axis=1)
    h_right = np.insert(h_distances, -1, h_distances[:, -1], axis=1)

    poly_coords = np.array([(centroids[0] - h_left,  centroids[1] + v_up),
                            (centroids[0] + h_right, centroids[1] + v_up),
                            (centroids[0] + h_right, centroids[1] - v_down),
                            (centroids[0] - h_left,  centroids[1] - v_down)])

    # poly_coords_list = np.reshape(poly_coords, (4, 2, -1))
    number_of_rows = len(poly_coords[0, 0, :, 0])
    number_of_columns = len(poly_coords[0, 0, 0, :])

    # convert to shapely objects
    feature_list = []
    for i, j in product(range(number_of_rows), range(number_of_columns)):
        coords = poly_coords[:, :, i, j]
        if not np.any(np.isnan(coords)):  # this alternative may be faster: if not np.isnan(np.sum(coords)):
            poly = Polygon(coords)
            if coord_transform_out is not None:
                ogr_poly = ogr.CreateGeometryFromWkb(poly.wkb)
                ogr_poly.Transform(coord_transform_out)
                poly = shapely_load_wkt(ogr_poly.ExportToWkt())
            feature_list.append({'polygon': poly,
                                 'lon': poly.centroid.coords[0][0],
                                 'lat': poly.centroid.coords[0][1],
                                 'lon_index': i,
                                 'lat_index': j})

    write_shapefile(feature_list)

    return feature_list


def write_shapefile(feature_list):
    from shapely.geometry import mapping
    import fiona
    schema = {
        'geometry': 'Polygon',
        'properties': {'id': 'int'},
    }
    with fiona.open('voronoi_grid.shp', 'w', layer='grid7', driver='ESRI Shapefile', schema=schema) as c:
        for x, p in enumerate(feature_list):
            c.write({'geometry': mapping(p['polygon']),
                     'properties': {'id': x},
                     })