# -*- coding: utf-8 -*-
"""
    RAPIDpy.gis

    Created by Alan D Snow, 2016.
    Based on RAPID_Toolbox for ArcMap
    License: BSD 3-Clause
"""
from osgeo import ogr


def open_shapefile(shapefile_path, file_geodatabase=None):
    """Opens a shapefile using either a shapefile path
        or a file geodatabase
    """
    if file_geodatabase:
        gdb_driver = ogr.GetDriverByName("OpenFileGDB")
        ogr_shapefile = gdb_driver.Open(file_geodatabase)
        ogr_shapefile_lyr = ogr_shapefile.GetLayer(shapefile_path)
    else:
        ogr_shapefile = ogr.Open(shapefile_path)
        ogr_shapefile_lyr = ogr_shapefile.GetLayer()
    return ogr_shapefile_lyr, ogr_shapefile
