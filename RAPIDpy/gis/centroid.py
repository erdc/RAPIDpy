# -*- coding: utf-8 -*-
"""
    centroid.py
    RAPIDpy

    Created by Alan D Snow, 2016.
    License: BSD 3-Clause
"""
from csv import writer as csv_writer

from osgeo import gdal, osr

# local
from . import open_shapefile
from ..helper_functions import open_csv

# Enable GDAL/OGR exceptions
gdal.UseExceptions()


def FlowlineToPoint(in_drainage_line,
                    river_id,
                    out_csv_file,
                    file_geodatabase=None):
    """
    Converts flowline feature to a list of centroid points with their rivid
    in EPSG:4326.

    Parameters
    ----------
    in_drainage_line: str
        Path to the stream network (i.e. Drainage Line) shapefile.
    river_id: str
        The name of the field with the river ID
        (Ex. 'HydroID', 'COMID', or 'LINKNO').
    out_csv_file: str
        Path to the output csv file with the centroid points.
    file_geodatabase: str, optional
        Path to the file geodatabase. If you use this option, in_drainage_line
         is the name of the stream network feature class
         (WARNING: Not always stable with GDAL).


    Example::

        from RAPIDpy.gis.centroid import FlowlineToPoint

        FlowlineToPoint(
            in_drainage_line='/path/to/drainageline.shp',
            river_id='LINKNO',
            out_csv_file='/path/to/comid_lat_lon_z.csv')

    """
    ogr_drainage_line_shapefile_lyr, ogr_drainage_line_shapefile = \
        open_shapefile(in_drainage_line, file_geodatabase)

    ogr_drainage_line_shapefile_lyr_proj = \
        ogr_drainage_line_shapefile_lyr.GetSpatialRef()
    osr_geographic_proj = osr.SpatialReference()
    osr_geographic_proj.ImportFromEPSG(4326)
    proj_transform = None
    if ogr_drainage_line_shapefile_lyr_proj != osr_geographic_proj:
        proj_transform = osr.CoordinateTransformation(
            ogr_drainage_line_shapefile_lyr_proj, osr_geographic_proj)

    # print valid field names to table
    with open_csv(out_csv_file, 'w') as outfile:
        writer = csv_writer(outfile)
        writer.writerow(['rivid', 'lat', 'lon', 'z'])
        for feature in ogr_drainage_line_shapefile_lyr:
            feat_geom = feature.GetGeometryRef()
            if proj_transform:
                feat_geom.Transform(proj_transform)
            centroid = feat_geom.Centroid()
            centroid_pt = centroid.GetPoint(0)
            writer.writerow([
                feature.GetField(river_id),
                centroid_pt[1],
                centroid_pt[0],
                centroid_pt[2]
            ])

    del ogr_drainage_line_shapefile
