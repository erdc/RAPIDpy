# -*- coding: utf-8 -*-
"""
    taudem.py
    RAPIDpy

    Created by Alan D Snow, 2016.
    Command line function based on ArcGIS python scripts by David Tarboton
    (https://github.com/dtarb/TauDEM)
    License: BSD 3-Clause
"""
from datetime import datetime
from multiprocessing import cpu_count
import os
from subprocess import PIPE, Popen
from sys import getrecursionlimit, setrecursionlimit

from gazar.grid import GDALGrid
import numpy as np
from past.builtins import xrange  # pylint: disable=redefined-builtin
from pyproj import Geod
from shapely.wkb import loads as shapely_loads
from shapely.ops import cascaded_union
from osgeo import gdal, ogr, osr

from ..helper_functions import log


# -----------------------------------------------------------------------------
# MAIN CLASS
# -----------------------------------------------------------------------------
class TauDEM(object):
    """
    TauDEM process manager.

    Attributes
    ----------
    taudem_exe_path: str, optional
        Path to TauDEM directory containing executables. This is requred to
        use TauDEM functionality.
    num_processors: int, optional
        Number of proessors to use with TauDEM. It only works if
        use_all_processors=False.
    use_all_processors: bool, optional
        If True, the TauDEM processes will use all avaialble processors.
    mpiexec_path: str, optional
        Path to mpiexec command. Default is 'mpiexec'.


    Initialization Example:

    .. code:: python

        from RAPIDpy.gis.taudem import TauDEM

        td = TauDEM("/path/to/scripts/TauDEM", use_all_processors=True)

    """
    def __init__(self,
                 taudem_exe_path="",
                 num_processors=1,
                 use_all_processors=False,
                 mpiexec_path="mpiexec"):
        """
        Initializer
        """
        if use_all_processors or num_processors > cpu_count():
            num_processors = cpu_count()

        self.taudem_exe_path = taudem_exe_path
        self.num_processors = num_processors
        self.mpiexec_path = mpiexec_path

        # other attributes
        self.pit_filled_elevation_grid = None
        self.flow_dir_grid = None
        self.contributing_area_grid = None
        self.stream_raster_grid = None

    def _run_mpi_cmd(self, cmd):
        """
        This runs the command you send in
        """
        log("Number of Processes: {0}".format(self.num_processors))
        time_start = datetime.utcnow()

        # Construct the taudem command line.
        cmd = [self.mpiexec_path, '-n', str(self.num_processors)] + cmd
        log("Command Line: {0}".format(" ".join(cmd)))
        process = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=False)
        out, err = process.communicate()
        if out:
            log("OUTPUT:")
            for line in out.split(b'\n'):
                log(line)
        if err:
            log(err, severity="WARNING")
        log("Time to complete: {0}".format(datetime.utcnow()-time_start))

    @staticmethod
    def _add_prj_file(original_gis_file, new_gis_file):
        """
        Adds projection file
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

    @staticmethod
    def extractSubNetwork(network_file,
                          out_subset_network_file,
                          outlet_ids,
                          river_id_field,
                          next_down_id_field,
                          river_magnitude_field,
                          safe_mode=True):
        """
        Extracts a subset river network from the main river network based on
        the outlet IDs.

        Parameters
        ----------
        network_file: str
            Path to the stream network shapefile.
        out_subset_network_file: str
            Path to the output subset stream network shapefile.
        outlet_ids: list
            List of integers reperesenting the outlet IDs to be included in
            the subset stream network.
        river_id_field: str
            Name of the river ID field in the stream network shapefile.
        next_down_id_field: str
            Name if the field with the river ID of the next downstream
            river segment in the stream network shapefile.
        river_magnitude_field: str
            Name of the river magnitude field in the stream network shapefile.
        safe_mode: bool, optional
            If True, it will kill the simulation early before over taxing
            your computer. If you are confident your computer can handle it,
            set it to False.


        Here is an example of how to use this:

        .. code:: python

            import os
            from RAPIDpy.gis.taudem import TauDEM


            output_directory = '/path/to/output/files'
            network_shp = os.path.join(output_directory,
                                       "stream_reach_file.shp")
            out_shp = os.path.join(output_directory,
                                   "stream_reach_file_subset.shp")

            TauDEM.extractSubNetwork(
                network_file=network_shp,
                out_subset_network_file=out_shp,
                outlet_ids=[60830],
                river_id_field="LINKNO",
                next_down_id_field="DSLINKNO",
                river_magnitude_field="Magnitude",
            )

        """
        network_shapefile = ogr.Open(network_file)
        network_layer = network_shapefile.GetLayer()
        number_of_features = network_layer.GetFeatureCount()
        network_layer_defn = network_layer.GetLayerDefn()
        rivid_list = np.zeros(number_of_features, dtype=np.int32)
        next_down_rivid_list = np.zeros(number_of_features, dtype=np.int32)
        for feature_idx, drainage_line_feature in enumerate(network_layer):
            rivid_list[feature_idx] = \
                drainage_line_feature.GetField(river_id_field)
            next_down_rivid_list[feature_idx] = \
                drainage_line_feature.GetField(next_down_id_field)

        def getSubNetworkIDList(outlet_river_id,
                                _rivid_list,
                                _next_down_rivid_list):
            """
            Adds ids upstream of the outlet to a list
            """
            sub_network_index_list = []
            try:
                for feature_ii in \
                        np.where(_next_down_rivid_list == outlet_river_id)[0]:
                    sub_network_index_list.append(feature_ii)
                    sub_network_index_list += \
                        getSubNetworkIDList(_rivid_list[feature_ii],
                                            _rivid_list,
                                            _next_down_rivid_list)
            except IndexError:
                pass
            return sub_network_index_list

        original_recursion_limit = getrecursionlimit()
        try:
            main_sub_network_index_list = []
            for outlet_id in outlet_ids:
                outlet_index = np.where(rivid_list == outlet_id)[0][0]
                outlet_feature = network_layer.GetFeature(outlet_index)
                outlet_magnitude = \
                    outlet_feature.GetField(river_magnitude_field)
                if outlet_magnitude > original_recursion_limit:
                    if not safe_mode:
                        setrecursionlimit(outlet_magnitude)
                    else:
                        raise Exception("Current recursion limit {0} will not "
                                        "allow extraction for stream magnitude"
                                        " {1}. To override, set safe_mode to "
                                        "False ..."
                                        .format(original_recursion_limit,
                                                outlet_magnitude))
                main_sub_network_index_list.append(outlet_index)
                main_sub_network_index_list += \
                    getSubNetworkIDList(outlet_id,
                                        rivid_list,
                                        next_down_rivid_list)
        except Exception:
            setrecursionlimit(original_recursion_limit)
            raise

        setrecursionlimit(original_recursion_limit)

        # Write out subset to new shapefile
        shp_drv = ogr.GetDriverByName('ESRI Shapefile')
        # Remove output shapefile if it already exists
        if os.path.exists(out_subset_network_file):
            shp_drv.DeleteDataSource(out_subset_network_file)

        network_subset_shp = shp_drv.CreateDataSource(out_subset_network_file)
        network_subset_layer = \
            network_subset_shp.CreateLayer('',
                                           network_layer.GetSpatialRef(),
                                           ogr.wkbLineString)
        # Add input Layer Fields to the output Layer if it is the one we want
        for iii in xrange(network_layer_defn.GetFieldCount()):
            network_subset_layer.CreateField(
                network_layer_defn.GetFieldDefn(iii))
        network_subset_layer_defn = network_subset_layer.GetLayerDefn()

        for feature_index in main_sub_network_index_list:
            subset_feature = network_layer.GetFeature(feature_index)
            # add to list
            new_feat = ogr.Feature(network_subset_layer_defn)

            # Add field values from input Layer
            for iii in xrange(network_layer_defn.GetFieldCount()):
                new_feat.SetField(
                    network_subset_layer_defn.GetFieldDefn(iii).GetNameRef(),
                    subset_feature.GetField(iii))

            # Set geometry as centroid
            geom = subset_feature.GetGeometryRef()
            new_feat.SetGeometry(geom.Clone())
            # Add new feature to output Layer
            network_subset_layer.CreateFeature(new_feat)

    @classmethod
    def extractLargestSubNetwork(cls,
                                 network_file,
                                 out_subset_network_file,
                                 river_id_field,
                                 next_down_id_field,
                                 river_magnitude_field,
                                 safe_mode=True):
        """
        Extracts the larges sub network from the watershed based on the
        magnitude parameter.

        Parameters
        ----------
        network_file: str
            Path to the stream network shapefile.
        out_subset_network_file: str
            Path to the output subset stream network shapefile.
        river_id_field: str
            Name of the river ID field in the stream network shapefile.
        next_down_id_field: str
            Name of the field with the river ID of the next downstream river
            segment in the stream network shapefile.
        river_magnitude_field: str
            Name of the river magnitude field in the stream network shapefile.
        safe_mode: bool, optional
            If True, it will kill the simulation early before over taxing
            your computer. If you are confident your computer can handle it,
            set it to False.


        Here is an example of how to use this:

        .. code:: python

            import os
            from RAPIDpy.gis.taudem import TauDEM

            output_directory = '/path/to/output/files'
            network_shp = os.path.join(output_directory,
                                       "stream_reach_file.shp")
            out_shp = os.path.join(output_directory,
                                   "stream_reach_file_subset.shp")

            TauDEM.extractLargestSubNetwork(
                network_file=network_shp,
                out_subset_network_file=out_shp,
                river_id_field="LINKNO",
                next_down_id_field="DSLINKNO",
                river_magnitude_field="Magnitude",
            )
        """
        network_shapefile = ogr.Open(network_file)
        network_layer = network_shapefile.GetLayer()
        number_of_features = network_layer.GetFeatureCount()
        riv_magnuitude_list = np.zeros(number_of_features, dtype=np.int32)
        for feature_idx, drainage_line_feature in enumerate(network_layer):
            riv_magnuitude_list[feature_idx] =\
                drainage_line_feature.GetField(river_magnitude_field)

        max_magnitude_feature = \
            network_layer.GetFeature(np.argmax(riv_magnuitude_list))
        cls.extractSubNetwork(network_file,
                              out_subset_network_file,
                              [max_magnitude_feature.GetField(river_id_field)],
                              river_id_field,
                              next_down_id_field,
                              river_magnitude_field,
                              safe_mode)

    @staticmethod
    def extractSubsetFromWatershed(subset_network_file,
                                   subset_network_river_id_field,
                                   watershed_file,
                                   watershed_network_river_id_field,
                                   out_watershed_subset_file):
        """
        Extract catchment by using subset network file.
        Use this after using either
        :func:`~RAPIDpy.gis.taudem.TauDEM.extractSubNetwork()`
        or :func:`~RAPIDpy.gis.taudem.TauDEM.extractLargestSubNetwork()`.

        Parameters
        ----------
        subset_network_file: str
            Path to the pre-subsetted stream network shapefile.
        subset_network_river_id_field: str
            The field name with the river ID in the stream network shapefile.
        watershed_file: str
            Path to the watershed shapefile.
        watershed_network_river_id_field: str
            Name of the field with the river ID in the watershed shapefile.
        out_watershed_subset_file: str
            The path to output the subset watershed shapefile.


        Here is an example of how to use this:

        .. code:: python

            import os
            from RAPIDpy.gis.taudem import TauDEM

            output_directory = '/path/to/output/files'
            network_shp = os.path.join(output_directory,
                                       "stream_reach_file.shp")
            water_shp = os.path.join(output_directory,
                                    "watershed_shapefile.shp")
            out_shp = os.path.join(output_directory,
                                   "watershed_shapefile_subset.shp")
            TauDEM.extractSubsetFromWatershed(
                subset_network_filenetwork_shp,
                subset_network_river_id_field="LINKNO",
                watershed_file=water_shp,
                watershed_network_river_id_field="LINKNO",
                out_watershed_subset_file=out_shp)

        """
        subset_network_shapefile = ogr.Open(subset_network_file)
        subset_network_layer = subset_network_shapefile.GetLayer()

        ogr_watershed_shapefile = ogr.Open(watershed_file)
        ogr_watershed_shapefile_lyr = ogr_watershed_shapefile.GetLayer()
        ogr_watershed_shapefile_lyr_defn = \
            ogr_watershed_shapefile_lyr.GetLayerDefn()

        number_of_features = ogr_watershed_shapefile_lyr.GetFeatureCount()
        watershed_rivid_list = np.zeros(number_of_features, dtype=np.int32)
        for feature_idx, watershed_feature in \
                enumerate(ogr_watershed_shapefile_lyr):
            watershed_rivid_list[feature_idx] = \
                watershed_feature.GetField(watershed_network_river_id_field)

        shp_drv = ogr.GetDriverByName('ESRI Shapefile')
        # Remove output shapefile if it already exists
        if os.path.exists(out_watershed_subset_file):
            shp_drv.DeleteDataSource(out_watershed_subset_file)

        subset_watershed_shapefile = \
            shp_drv.CreateDataSource(out_watershed_subset_file)
        subset_watershed_layer = \
            subset_watershed_shapefile.CreateLayer(
                '',
                ogr_watershed_shapefile_lyr.GetSpatialRef(),
                ogr.wkbPolygon)
        # Add input Layer Fields to the output Layer if it is the one we want
        for iii in xrange(ogr_watershed_shapefile_lyr_defn.GetFieldCount()):
            subset_watershed_layer.CreateField(
                ogr_watershed_shapefile_lyr_defn.GetFieldDefn(iii))
        subset_watershed_layer_defn = subset_watershed_layer.GetLayerDefn()

        for drainage_line_feature in subset_network_layer:
            try:
                watershed_feature_index = \
                    np.where(watershed_rivid_list ==
                             drainage_line_feature.GetField(
                                 subset_network_river_id_field))[0][0]
            except IndexError:
                log("{0} {1} not found ...".format(
                    subset_network_river_id_field,
                    drainage_line_feature.GetField(
                        subset_network_river_id_field)))
                continue

            subset_feature = \
                ogr_watershed_shapefile_lyr.GetFeature(watershed_feature_index)
            # add to list
            new_feat = ogr.Feature(subset_watershed_layer_defn)

            # Add field values from input Layer
            for iii in \
                    xrange(ogr_watershed_shapefile_lyr_defn.GetFieldCount()):
                new_feat.SetField(
                    subset_watershed_layer_defn.GetFieldDefn(iii).GetNameRef(),
                    subset_feature.GetField(iii))

            # Set geometry as centroid
            geom = subset_feature.GetGeometryRef()
            new_feat.SetGeometry(geom.Clone())
            # Add new feature to output Layer
            subset_watershed_layer.CreateFeature(new_feat)

    @staticmethod
    def rasterToPolygon(raster_file, polygon_file):
        """
        Converts watershed raster to polygon and then dissolves it.
        It dissolves features based on the LINKNO attribute.
        """
        log("Process: Raster to Polygon ...")
        time_start = datetime.utcnow()
        temp_polygon_file = \
            "{0}_temp.shp".format(
                os.path.splitext(os.path.basename(polygon_file))[0])

        GDALGrid(raster_file).to_polygon(out_shapefile=temp_polygon_file,
                                         fieldname="LINKNO",
                                         self_mask=True)

        log("Time to convert to polygon: {0}"
            .format(datetime.utcnow()-time_start))

        log("Dissolving ...")
        time_start_dissolve = datetime.utcnow()
        ogr_polygin_shapefile = ogr.Open(temp_polygon_file)
        ogr_polygon_shapefile_lyr = ogr_polygin_shapefile.GetLayer()
        number_of_features = ogr_polygon_shapefile_lyr.GetFeatureCount()
        polygon_rivid_list = np.zeros(number_of_features, dtype=np.int32)
        for feature_idx, catchment_feature in \
                enumerate(ogr_polygon_shapefile_lyr):
            polygon_rivid_list[feature_idx] = \
                catchment_feature.GetField('LINKNO')

        shp_drv = ogr.GetDriverByName('ESRI Shapefile')
        # Remove output shapefile if it already exists
        if os.path.exists(polygon_file):
            shp_drv.DeleteDataSource(polygon_file)

        dissolve_shapefile = shp_drv.CreateDataSource(polygon_file)
        dissolve_layer = \
            dissolve_shapefile.CreateLayer(
                '',
                ogr_polygon_shapefile_lyr.GetSpatialRef(),
                ogr.wkbPolygon)
        dissolve_layer.CreateField(ogr.FieldDefn('LINKNO', ogr.OFTInteger))
        dissolve_layer_defn = dissolve_layer.GetLayerDefn()

        for unique_rivid in np.unique(polygon_rivid_list):
            # get indices where it is in the polygon
            feature_indices = np.where(polygon_rivid_list == unique_rivid)[0]
            new_feat = ogr.Feature(dissolve_layer_defn)
            new_feat.SetField('LINKNO', int(unique_rivid))

            if len(feature_indices) == 1:
                # write feature to file
                feature = \
                    ogr_polygon_shapefile_lyr.GetFeature(feature_indices[0])
                new_feat.SetGeometry(feature.GetGeometryRef())
            else:
                # dissolve
                dissolve_poly_list = []
                for feature_index in feature_indices:
                    feature = \
                        ogr_polygon_shapefile_lyr.GetFeature(feature_index)
                    feat_geom = feature.GetGeometryRef()
                    dissolve_poly_list.append(
                        shapely_loads(feat_geom.ExportToWkb()))
                dissolve_polygon = cascaded_union(dissolve_poly_list)
                new_feat.SetGeometry(
                    ogr.CreateGeometryFromWkb(dissolve_polygon.wkb))
            dissolve_layer.CreateFeature(new_feat)
        # clean up
        shp_drv.DeleteDataSource(temp_polygon_file)
        log("Time to dissolve: {0}".format(datetime.utcnow() -
                                           time_start_dissolve))
        log("Total time to convert: {0}".format(datetime.utcnow() -
                                                time_start))

    @staticmethod
    def addLengthMeters(stream_network):
        """
        Adds length field in meters to network
        (The added field name will be 'LENGTH_M').

        .. note:: This may be needed for generating the kfac file
                  depending on the units of your raster. See: :doc:`gis_tools`.

        Parameters
        ----------
        stream_network: str
            Path to stream network file.


        Here is an example of how to use this:

        .. code:: python

            import os
            from RAPIDpy.gis.taudem import TauDEM

            output_directory = '/path/to/output/files'
            TauDEM.addLengthMeters(os.path.join(output_directory,
                                                "stream_reach_file.shp"))

        """
        network_shapefile = ogr.Open(stream_network, 1)
        network_layer = network_shapefile.GetLayer()
        network_layer_defn = network_layer.GetLayerDefn()

        # make sure projection EPSG:4326
        network_layer_proj = network_layer.GetSpatialRef()
        geographic_proj = osr.SpatialReference()
        geographic_proj.ImportFromEPSG(4326)
        proj_transform = None
        if network_layer_proj != geographic_proj:
            proj_transform = osr.CoordinateTransformation(network_layer_proj,
                                                          geographic_proj)

        # check for field
        create_field = True
        for i in xrange(network_layer_defn.GetFieldCount()):
            field_name = network_layer_defn.GetFieldDefn(i).GetName()
            if field_name == 'LENGTH_M':
                create_field = False
                break

        if create_field:
            network_layer.CreateField(ogr.FieldDefn('LENGTH_M', ogr.OFTReal))

        geo_manager = Geod(ellps="WGS84")
        for network_feature in network_layer:
            feat_geom = network_feature.GetGeometryRef()
            # make sure coordinates are geographic
            if proj_transform:
                feat_geom.Transform(proj_transform)

            line = shapely_loads(feat_geom.ExportToWkb())
            lon_list, lat_list = line.xy
            dist = geo_manager.inv(lon_list[:-1], lat_list[:-1],
                                   lon_list[1:], lat_list[1:])[2]
            network_feature.SetField('LENGTH_M', sum(dist))
            network_layer.SetFeature(network_feature)

    def pitRemove(self,
                  elevation_grid,
                  pit_filled_elevation_grid,
                  input_depression_mask_grid=None,
                  consider4way=False,
                  ):
        """
        Remove low spots from DEM.
        """
        log("PROCESS: PitRemove")
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

        # create projection file
        self._add_prj_file(elevation_grid,
                           self.pit_filled_elevation_grid)

    def dinfFlowDirection(self,
                          flow_dir_grid,
                          slope_grid,
                          pit_filled_elevation_grid=None):
        """
        Calculates flow direction with Dinf method
        """
        log("PROCESS: DinfFlowDirection")
        if pit_filled_elevation_grid:
            self.pit_filled_elevation_grid = pit_filled_elevation_grid

        # Construct the taudem command line.
        cmd = [os.path.join(self.taudem_exe_path, 'dinfflowdir'),
               '-fel', self.pit_filled_elevation_grid,
               '-ang', flow_dir_grid,
               '-slp', slope_grid,
               ]

        self._run_mpi_cmd(cmd)

        # create projection files
        self._add_prj_file(self.pit_filled_elevation_grid,
                           flow_dir_grid)
        self._add_prj_file(self.pit_filled_elevation_grid,
                           slope_grid)

    def d8FlowDirection(self,
                        flow_dir_grid,
                        slope_grid,
                        pit_filled_elevation_grid=None):
        """
        Calculates flow direction with D8 method.
        """
        log("PROCESS: D8FlowDirection")
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

        # create projection files
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
        Calculates contributing area with Dinf method.
        """
        log("PROCESS: DinfContributingArea")

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

        # create projection file
        self._add_prj_file(flow_dir_grid,
                           contributing_area_grid)

    def d8ContributingArea(self,
                           contributing_area_grid,
                           outlet_shapefile=None,
                           weight_grid=None,
                           edge_contamination=False,
                           flow_dir_grid=None):
        """
        Calculates contributing area with D8 method.
        """
        log("PROCESS: D8ContributingArea")
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

        # create projection file
        self._add_prj_file(self.flow_dir_grid,
                           self.contributing_area_grid)

    def streamDefByThreshold(self,
                             stream_raster_grid,
                             threshold,
                             contributing_area_grid,
                             mask_grid=None,
                             ):
        """
        Calculates the stream definition by threshold.
        """
        log("PROCESS: StreamDefByThreshold")
        self.stream_raster_grid = stream_raster_grid

        # Construct the taudem command line.
        cmd = [os.path.join(self.taudem_exe_path, 'threshold'),
               '-ssa', contributing_area_grid,
               '-src', self.stream_raster_grid,
               '-thresh', str(threshold),
               ]

        if mask_grid:
            cmd += ['-mask', mask_grid]

        self._run_mpi_cmd(cmd)

        # create projection file
        self._add_prj_file(contributing_area_grid,
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
        log("PROCESS: StreamReachAndWatershed")
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

        # create projection file
        self._add_prj_file(self.pit_filled_elevation_grid,
                           out_stream_order_grid)
        self._add_prj_file(self.pit_filled_elevation_grid,
                           out_stream_reach_file)
        self._add_prj_file(self.pit_filled_elevation_grid,
                           out_watershed_grid)

    def demToStreamNetwork(self,
                           output_directory,
                           raw_elevation_dem="",
                           pit_filled_elevation_grid="",
                           flow_dir_grid_d8="",
                           contributing_area_grid_d8="",
                           flow_dir_grid_dinf="",
                           contributing_area_grid_dinf="",
                           use_dinf=False,
                           threshold=1000,
                           delineate=False):
        """
        This function will run all of the TauDEM processes to generate
        a stream network from an elevation raster.

        .. note:: For information about the stream reach and watershed process:
                  http://hydrology.usu.edu/taudem/taudem5/help53/StreamReachAndWatershed.html

        .. note:: For information about the *threshold* parameter, see:
                  http://hydrology.usu.edu/taudem/taudem5/help53/StreamDefinitionByThreshold.html

        Parameters
        ----------
        output_directory: str
            Path to output generated files to.
        raw_elevation_dem: str, optional
            Path to original elevation DEM file. Required if
            *pit_filled_elevation_grid* is not used.
        pit_filled_elevation_grid: str, optional
            Path to pit filled elevation DEM file. Required if
            *raw_elevation_dem* is not used.
        flow_dir_grid_d8: str, optional
            Path to flow direction grid generated using TauDEM's D8 method.
        contributing_area_grid_d8: str, optional
            Path to contributing area grid generated using TauDEM's D8 method.
        flow_dir_grid_dinf: str, optional
            Path to flow direction grid generated using TauDEM's
            D-Infinity method (EXPERIMENTAL).
        contributing_area_grid_dinf: str, optional
            Path to contributing area grid generated using TauDEM's
            D-Infinity method (EXPERIMENTAL).
        use_dinf: bool, optional
            Use the D-Infinity method to get stream definition (EXPERIMENTAL).
        threshold: int, optional
            The stream threshold or maximum number of upstream grid cells.
            See above note.
        delineate: bool, optional
            If True, this will use the delineate option for theis method
            using TauDEM. Default is False.


        Here is an example of how to use this:

        .. code:: python

            from RAPIDpy.gis.taudem import TauDEM


            elevation_dem = '/path/to/dem.tif'
            output_directory = '/path/to/output/files'

            td = TauDEM("/path/to/scripts/TauDEM")
            td.demToStreamNetwork(output_directory,
                                  elevation_dem,
                                  threshold=1000)

        """
        time_start = datetime.utcnow()

        # FILL PITS IF NEEDED
        self.pit_filled_elevation_grid = pit_filled_elevation_grid
        if not pit_filled_elevation_grid:
            pit_filled_elevation_grid = \
                os.path.join(output_directory, 'pit_filled_elevation_grid.tif')
            self.pitRemove(raw_elevation_dem,
                           pit_filled_elevation_grid)

        # GENERATE D8 RASTERS
        self.flow_dir_grid = flow_dir_grid_d8
        if not flow_dir_grid_d8:
            flow_dir_grid_d8 = \
                os.path.join(output_directory, 'flow_dir_grid_d8.tif')
            slope_grid_d8 = os.path.join(output_directory, 'slope_grid_d8.tif')
            self.d8FlowDirection(flow_dir_grid_d8,
                                 slope_grid_d8)

        self.contributing_area_grid = contributing_area_grid_d8
        if not contributing_area_grid_d8:
            contributing_area_grid_d8 = \
                os.path.join(output_directory, 'contributing_area_grid_d8.tif')
            self.d8ContributingArea(contributing_area_grid_d8)

        stream_raster_grid = \
            os.path.join(output_directory, 'stream_raster_grid.tif')
        if use_dinf:
            log("USING DINF METHOD TO GET STREAM DEFINITION ...")
            if not flow_dir_grid_dinf:
                flow_dir_grid_dinf = \
                    os.path.join(output_directory, 'flow_dir_grid_dinf.tif')
                slope_grid_dinf = \
                    os.path.join(output_directory, 'slope_grid_dinf.tif')
                self.dinfFlowDirection(flow_dir_grid_dinf,
                                       slope_grid_dinf)
            if not contributing_area_grid_dinf:
                contributing_area_grid_dinf = \
                    os.path.join(output_directory,
                                 'contributing_area_grid_dinf.tif')
                self.dinfContributingArea(contributing_area_grid_dinf,
                                          flow_dir_grid_dinf)

            self.streamDefByThreshold(stream_raster_grid,
                                      threshold,
                                      contributing_area_grid_dinf)
        else:
            log("USING D8 METHOD TO GET STREAM DEFINITION ...")
            self.streamDefByThreshold(stream_raster_grid,
                                      threshold,
                                      contributing_area_grid_d8)

        # GENERATE STREAM NETWORK
        out_stream_order_grid = \
            os.path.join(output_directory, 'stream_order_grid.tif')
        out_network_connectivity_tree = \
            os.path.join(output_directory, 'network_connectivity_tree.txt')
        out_network_coordinates = \
            os.path.join(output_directory, 'network_coordinates.txt')
        out_stream_reach_file = \
            os.path.join(output_directory, 'stream_reach_file.shp')
        out_watershed_grid = \
            os.path.join(output_directory, 'watershed_grid.tif')
        self.streamReachAndWatershed(delineate,
                                     out_stream_order_grid,
                                     out_network_connectivity_tree,
                                     out_network_coordinates,
                                     out_stream_reach_file,
                                     out_watershed_grid)

        # convert watersed grid to shapefile
        out_watershed_shapefile = \
            os.path.join(output_directory, 'watershed_shapefile.shp')
        self.rasterToPolygon(out_watershed_grid, out_watershed_shapefile)
        log("Total time to complete: {0}".format(datetime.utcnow() -
                                                 time_start))
