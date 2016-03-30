# -*- coding: utf-8 -*-
##
##  taudem.py
##  RAPIDpy
##
##  Created by Alan D Snow.
##  Based on ArcGIS python scripts by David Tarboton (https://github.com/dtarb/TauDEM)
##  Copyright © 2016 Alan D Snow. All rights reserved.
##  License: BSD 3-Clause

from datetime import datetime
from multiprocessing import cpu_count
import numpy as np
import os
from subprocess import PIPE, Popen

try:
    from osgeo import gdal, ogr
    from shapely.wkb import loads as shapely_loads
    from shapely.ops import cascaded_union
except ImportError:
    raise Exception("You need to install the gdal and shapely python packages to use this tool ...")
    
class TauDEM(object):
    """
    TauDEM process manager
    """
    def __init__(self,
                 taudem_exe_path="", 
                 num_processors=1,
                 use_all_processors=False):
        """
        Initializer
        """             
        if use_all_processors or num_processors > cpu_count:
            num_processors = cpu_count()

        self.taudem_exe_path = taudem_exe_path
        self.num_processors = num_processors

    def _run_mpi_cmd(self, cmd):
        """
        This runs the command you send in
        """
            # Get and describe the first argument
        #
    
        print("Number of Processes: {0}".format(self.num_processors))
        time_start = datetime.utcnow()

        # Construct the taudem command line.
        cmd = ['mpiexec', '-n', str(self.num_processors)] + cmd
        print("Command Line: {0}".format(" ".join(cmd)))
        process = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=False)
        out, err = process.communicate()
        if out:
            print("OUTPUT:")
            for line in out.split('\n'):
                print line
        if err:
            print("ERROR:")
            print(err)
            #raise Exception(err)
        print("Time to complete: {0}".format(datetime.utcnow()-time_start))
            
    def _add_prj_file(self, original_gis_file,
                      new_gis_file):
        """
        Adds prj file
        """
        out_prj_file = "{0}.prj".format(os.path.splitext(new_gis_file)[0])
        if original_gis_file.endswith(".shp"):
            dataset = ogr.Open(original_gis_file)
            layer = dataset.GetLayer()
            spatial_ref = layer.GetSpatialRef()
            spatial_ref.MorphToESRI()
            spatial_ref_str = spatial_ref.ExportToWkt()
        else:
            dataset = gdal.Open(original_gis_file)
            spatial_ref_str = dataset.GetProjection()
            
        with open(out_prj_file, 'w') as prj_file:
            prj_file.write(spatial_ref_str)
    
    def rasterToPolygon(self, raster_file, polygon_file):
        """
        Converts raster to polygon and then dissolves it
        """
        print("Process: Raster to Polygon ...")
        time_start = datetime.utcnow()
        temp_polygon_file = "{0}_temp.shp".format(os.path.splitext(os.path.basename(polygon_file))[0])
        cmd = ["gdal_polygonize.py", raster_file,
               "-f", "ESRI Shapefile", temp_polygon_file,  
               os.path.splitext(os.path.basename(temp_polygon_file))[0],
               "LINKNO"]

        process = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=False)
        out, err = process.communicate()
        if out:
            print("OUTPUT:")
            for line in out.split('\n'):
                print line
        if err:
            print("ERROR:")
            print(err)
            #raise Exception(err)
        print("Time to convert to polygon: {0}".format(datetime.utcnow()-time_start))
        
        print("Dissolving ...")
        time_start_dissolve = datetime.utcnow()
        ogr_polygin_shapefile = ogr.Open(temp_polygon_file)
        ogr_polygon_shapefile_lyr = ogr_polygin_shapefile.GetLayer()
        number_of_features = ogr_polygon_shapefile_lyr.GetFeatureCount()
        polygon_rivid_list = np.zeros(number_of_features, dtype=np.int32)
        for feature_idx, catchment_feature in enumerate(ogr_polygon_shapefile_lyr):
            polygon_rivid_list[feature_idx] = catchment_feature.GetField('LINKNO')
        
        
        shp_drv = ogr.GetDriverByName('ESRI Shapefile')
        dissolve_shapefile = shp_drv.CreateDataSource(polygon_file)
        dissolve_layer = dissolve_shapefile.CreateLayer('', ogr_polygon_shapefile_lyr.GetSpatialRef(), ogr.wkbPolygon)
        dissolve_layer.CreateField(ogr.FieldDefn('LINKNO', ogr.OFTInteger))
        dissolve_layer_defn = dissolve_layer.GetLayerDefn()

        for unique_rivid in np.unique(polygon_rivid_list):
            #get indices where it is in the polygon
            feature_indices = np.where(polygon_rivid_list==unique_rivid)[0]
            new_feat = ogr.Feature(dissolve_layer_defn)
            new_feat.SetField('LINKNO', int(unique_rivid))

            if len(feature_indices) == 1:
                ##write feature to file
                feature = ogr_polygon_shapefile_lyr.GetFeature(feature_indices[0])
                new_feat.SetGeometry(feature.GetGeometryRef())
            else:
                ##dissolve
                dissolve_poly_list = []
                for feature_index in feature_indices:
                    feature = ogr_polygon_shapefile_lyr.GetFeature(feature_index)
                    feat_geom = feature.GetGeometryRef()
                    dissolve_poly_list.append(shapely_loads(feat_geom.ExportToWkb()))                
                dissolve_polygon = cascaded_union(dissolve_poly_list)
                new_feat.SetGeometry(ogr.CreateGeometryFromWkb(dissolve_polygon.wkb))
            dissolve_layer.CreateFeature(new_feat)
        #clean up
        shp_drv.DeleteDataSource(temp_polygon_file)
        print("Time to dissolve: {0}".format(datetime.utcnow()-time_start_dissolve))
        print("Total time to convert: {0}".format(datetime.utcnow()-time_start))
        
    def pitRemove(self, 
                  elevation_grid,
                  pit_filled_elevation_grid,
                  input_depression_mask_grid=None,
                  consider4way=False,
                  ):
        """
        Remove low spots from DEM
        """
        print("PROCESS: PitRemove")
        self.pit_filled_elevation_grid = pit_filled_elevation_grid
             
        # Construct the taudem command line.
        cmd = [os.path.join(self.taudem_exe_path, 'pitremove'),
               '-z', elevation_grid, 
               '-fel', self.pit_filled_elevation_grid,
               ]
    
        if input_depression_mask_grid:
            cmd += ['-depmask', input_depression_mask_grid]
        if consider4way:
            cmd += ['-4way']
        
        self._run_mpi_cmd(cmd)

        #create projection file
        self._add_prj_file(elevation_grid,
                           self.pit_filled_elevation_grid)
        
    def dinfFlowDirection(self, 
                          flow_dir_grid,
                          slope_grid,
                          pit_filled_elevation_grid):
        """
        Calculates flow direction with Dinf method
        """                     
        print("PROCESS: DinfFlowDirection")
        # Construct the taudem command line.
        cmd = [os.path.join(self.taudem_exe_path, 'dinfflowdir'),
               '-fel', pit_filled_elevation_grid, 
               '-ang', flow_dir_grid, 
               '-slp', slope_grid,
               ]
    
        self._run_mpi_cmd(cmd)

        #create projection files
        self._add_prj_file(pit_filled_elevation_grid,
                           flow_dir_grid)
        self._add_prj_file(pit_filled_elevation_grid,
                           slope_grid)

    def d8FlowDirection(self, 
                        flow_dir_grid,
                        slope_grid,
                        pit_filled_elevation_grid=None):
        """
        Calculates flow direction with D8 method
        """                     
        print("PROCESS: D8FlowDirection")
        if pit_filled_elevation_grid:
            self.pit_filled_elevation_grid = pit_filled_elevation_grid
    
        self.flow_dir_grid = flow_dir_grid
        
        # Construct the taudem command line.
        cmd = [os.path.join(self.taudem_exe_path, 'd8flowdir'),
               '-fel', self.pit_filled_elevation_grid, 
               '-p', self.flow_dir_grid, 
               '-sd8', slope_grid,
               ]
    
        self._run_mpi_cmd(cmd)

        #create projection files
        self._add_prj_file(self.pit_filled_elevation_grid,
                           self.flow_dir_grid)
        self._add_prj_file(self.pit_filled_elevation_grid,
                           slope_grid)
        
    def dinfContributingArea(self, 
                             contributing_area_grid,
                             flow_dir_grid,
                             outlet_shapefile=None,
                             weight_grid=None,
                             edge_contamination=False,
                             ):
        """
        Calculates contributing area with Dinf method
        """                     
        print("PROCESS: DinfContributingArea")

        # Construct the taudem command line.
        cmd = [os.path.join(self.taudem_exe_path, 'areadinf'),
               '-ang', flow_dir_grid, 
               '-sca', contributing_area_grid,
               ]
               
        if outlet_shapefile:
            cmd += ['-o', outlet_shapefile]
        if weight_grid:
            cmd += ['-wg', weight_grid]
        if not edge_contamination:
            cmd = cmd + ['-nc']    
    
        self._run_mpi_cmd(cmd)

        #create projection file
        self._add_prj_file(flow_dir_grid,
                           contributing_area_grid)
                           
    def d8ContributingArea(self, 
                           contributing_area_grid,
                           outlet_shapefile=None,
                           weight_grid=None,
                           edge_contamination=False,
                           flow_dir_grid=None):
        """
        Calculates contributing area with D8 method
        """                     
        print("PROCESS: D8ContributingArea")
        if flow_dir_grid:
            self.flow_dir_grid = flow_dir_grid

        self.contributing_area_grid = contributing_area_grid
        # Construct the taudem command line.
        cmd = [os.path.join(self.taudem_exe_path, 'aread8'),
               '-p', self.flow_dir_grid, 
               '-ad8', self.contributing_area_grid,
               ]
               
        if outlet_shapefile:
            cmd += ['-o', outlet_shapefile]
        if weight_grid:
            cmd += ['-wg', weight_grid]
        if not edge_contamination:
            cmd = cmd + ['-nc']    
    
        self._run_mpi_cmd(cmd)

        #create projection file
        self._add_prj_file(self.flow_dir_grid,
                           self.contributing_area_grid)

    def streamDefByThreshold(self, 
                             stream_raster_grid,
                             threshold,
                             contributing_area_grid=None,
                             mask_grid=None,
                             ):
        """
        Calculates the stream definition by threshold
        """
        print("PROCESS: StreamDefByThreshold")
        if contributing_area_grid:
            self.contributing_area_grid = contributing_area_grid
            
        self.stream_raster_grid = stream_raster_grid
        
        # Construct the taudem command line.
        cmd = [os.path.join(self.taudem_exe_path, 'threshold'),
               '-ssa', self.contributing_area_grid, 
               '-src', self.stream_raster_grid,
               '-thresh', str(threshold),
               ]

        if mask_grid:
            cmd += ['-mask', mask_grid]
            
        self._run_mpi_cmd(cmd)
        
        #create projection file
        self._add_prj_file(self.contributing_area_grid,
                           self.stream_raster_grid)
                           
    def streamReachAndWatershed(self, 
                                delineate,
                                out_stream_order_grid,
                                out_network_connectivity_tree,
                                out_network_coordinates,
                                out_stream_reach_file,
                                out_watershed_grid,
                                pit_filled_elevation_grid=None,
                                flow_dir_grid=None,
                                contributing_area_grid=None,
                                stream_raster_grid=None,
                                outlet_shapefile=None
                                ):
        """
        Creates vector network and shapefile from stream raster grid
        """
        print("PROCESS: StreamReachAndWatershed")
        if pit_filled_elevation_grid:
            self.pit_filled_elevation_grid = pit_filled_elevation_grid
        if flow_dir_grid:
            self.flow_dir_grid = flow_dir_grid
        if contributing_area_grid:
            self.contributing_area_grid = contributing_area_grid
        if stream_raster_grid:
            self.stream_raster_grid = stream_raster_grid
            
        # Construct the taudem command line.
        cmd = [os.path.join(self.taudem_exe_path, 'streamnet'),
               '-fel', self.pit_filled_elevation_grid, 
               '-p', self.flow_dir_grid,
               '-ad8', self.contributing_area_grid,
               '-src', self.stream_raster_grid,
               '-ord', out_stream_order_grid,
               '-tree', out_network_connectivity_tree,
               '-coord', out_network_coordinates,
               '-net', out_stream_reach_file,
               '-w', out_watershed_grid,
               ]
               
        if outlet_shapefile:
            cmd += ['-o', outlet_shapefile]
        if delineate:
            cmd += ['-sw']
            
        self._run_mpi_cmd(cmd)

        #create projection file
        self._add_prj_file(self.pit_filled_elevation_grid,
                           out_stream_reach_file)
        self._add_prj_file(self.pit_filled_elevation_grid,
                           out_watershed_grid)
                           
    def demToStreamNetwork(self, elevation_dem, output_directory, 
                           threshold=1000, delineate=False):
        """
        This function will run all of the processes to genrate a stream network
        from an elevation dem
        """

        time_start = datetime.utcnow()
        pit_filled_elevation_grid = os.path.join(output_directory, 'pit_filled_elevation_grid.tif')
        self.pitRemove(elevation_dem,
                       pit_filled_elevation_grid)
        flow_dir_grid = os.path.join(output_directory, 'flow_dir_grid.tif')
        slope_grid = os.path.join(output_directory, 'slope_grid.tif')
        contributing_area_grid = os.path.join(output_directory, 'contributing_area_grid.tif')
        self.d8FlowDirection(flow_dir_grid,
                             slope_grid)
        self.d8ContributingArea(contributing_area_grid)
        stream_raster_grid = os.path.join(output_directory, 'stream_raster_grid.tif')
        self.streamDefByThreshold(stream_raster_grid,
                                  threshold)
        out_stream_order_grid = os.path.join(output_directory, 'out_stream_order_grid.tif')
        out_network_connectivity_tree = os.path.join(output_directory, 'out_network_connectivity_tree.txt')
        out_network_coordinates = os.path.join(output_directory, 'out_network_coordinates.txt')
        out_stream_reach_file = os.path.join(output_directory, 'out_stream_reach_file.shp')
        out_watershed_grid = os.path.join(output_directory, 'out_watershed_grid.tif')
        self.streamReachAndWatershed(delineate,
                                     out_stream_order_grid,
                                     out_network_connectivity_tree,
                                     out_network_coordinates,
                                     out_stream_reach_file,
                                     out_watershed_grid)
                                     
        #convert watersed grid to shapefile
        polygon_file = os.path.join(output_directory, 'catchments.shp')                          
        self.rasterToPolygon(out_watershed_grid, polygon_file)
        print("TOTAL time to complete: {0}".format(datetime.utcnow()-time_start))