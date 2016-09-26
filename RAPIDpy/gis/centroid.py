# -*- coding: utf-8 -*-
##
##  centroid.py
##  RAPIDpy
##
##  Created by Alan D Snow.
##
##  Copyright © 2016 Alan D Snow. All rights reserved.
##  License: BSD 3-Clause

from csv import writer as csv_writer
try:
    from osgeo import gdal, ogr, osr
except Exception:
    raise Exception("You need the gdal python package to run this tool ...")
# Enable GDAL/OGR exceptions
gdal.UseExceptions()

#local
from ..helper_functions import open_csv

def FlowlineToPoint(in_drainage_line,
                    river_id,
                    out_csv_file,
                    file_geodatabase=None):
    """
    Converts flowline feature to a list of centroid points with their comid in EPSG:4326.

    Args:
        in_drainage_line(str): Path to the stream network (i.e. Drainage Line) shapefile.
        river_id(str): The name of the field with the river ID (Ex. 'HydroID', 'COMID', or 'LINKNO').
        out_csv_file(str): Path to the output csv file with the centroid points.
        file_geodatabase(Optional[str]): Path to the file geodatabase. If you use this option, in_drainage_line is the name of the stream network feature class. (WARNING: Not always stable with GDAL.)
    
    Example::
    
        from RAPIDpy.gis.centroid import FlowlineToPoint
        #------------------------------------------------------------------------------
        #main process
        #------------------------------------------------------------------------------
        if __name__ == "__main__":
            FlowlineToPoint(in_drainage_line='/path/to/drainageline.shp',
                            river_id='LINKNO',
                            out_csv_file='/path/to/comid_lat_lon_z.csv',
                            )
    
    """

    if file_geodatabase:
        gdb_driver = ogr.GetDriverByName("OpenFileGDB")
        ogr_file_geodatabase = gdb_driver.Open(file_geodatabase, 0)
        ogr_drainage_line_shapefile_lyr = ogr_file_geodatabase.GetLayer(in_drainage_line)
    else:
        ogr_drainage_line_shapefile = ogr.Open(in_drainage_line)
        ogr_drainage_line_shapefile_lyr = ogr_drainage_line_shapefile.GetLayer()
    
    ogr_drainage_line_shapefile_lyr_proj = ogr_drainage_line_shapefile_lyr.GetSpatialRef()
    osr_geographic_proj = osr.SpatialReference()
    osr_geographic_proj.ImportFromEPSG(4326)
    proj_transform = None
    if ogr_drainage_line_shapefile_lyr_proj != osr_geographic_proj:
        proj_transform = osr.CoordinateTransformation(ogr_drainage_line_shapefile_lyr_proj, osr_geographic_proj)

    #print valid field names to table
    with open_csv(out_csv_file, 'w') as outfile:
        writer = csv_writer(outfile)
        writer.writerow(['rivid','lat','lon','z'])
        for feature in ogr_drainage_line_shapefile_lyr:
            feat_geom = feature.GetGeometryRef()
            if proj_transform:
                feat_geom.Transform(proj_transform)
            centroid = feat_geom.Centroid()
            centroid_pt = centroid.GetPoint(0)
            writer.writerow([feature.GetField(river_id), centroid_pt[1], centroid_pt[0], centroid_pt[2]])