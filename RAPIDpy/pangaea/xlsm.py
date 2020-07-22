# -*- coding: utf-8 -*-
#
#  xlsm.py
#  pangaea
#
#  Author : Alan D Snow, 2017.
#  License: BSD 3-Clause
"""pangea.xlsm
    This module is an extension for xarray for land surface models.
    (see: http://xarray.pydata.org/en/stable/internals.html#extending-xarray)
"""
from affine import Affine
import numpy as np
from osgeo import osr, gdalconst
import pandas as pd
from pyproj import Proj, transform
from gazar.grid import (geotransform_from_yx, resample_grid,
                        utm_proj_from_latlon, ArrayGrid)
import wrf
import xarray as xr


@xr.register_dataset_accessor('lsm')
class LSMGridReader(object):
    """
    This is an extension for xarray specifically
    designed for land surface models.

    Read with pangaea example::

        import pangaea as pa

        with pa.open_mfdataset('/path/to/ncfiles/*.nc',
                               lat_var='lat',
                               lon_var='lon',
                               time_var='time',
                               lat_dim='lat',
                               lon_dim='lon',
                               time_dim='time') as xds:
            print(xds.lsm.projection)

    Read with xarray example::

        import xarray as xr

        with pa.open_dataset('/path/to/file.nc') as xds:
            print(xds.lsm.projection)
    """
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._projection = None
        self._epsg = None
        self._geotransform = None
        self._affine = None
        self._center = None
        self._y_inverted = None

        # set variable information
        self.y_var = 'lat'
        self.x_var = 'lon'
        self.time_var = 'time'
        # set dimension information
        self.y_dim = 'y'
        self.x_dim = 'x'
        self.time_dim = 'time'
        # convert lon from [0 to 360] to [-180 to 180]
        self.lon_to_180 = False
        # coordinates are projected already
        self.coords_projected = False

    def to_datetime(self):
        """Converts time to datetime."""
        time_values = self._obj[self.time_var].values
        if 'datetime' not in str(time_values.dtype):
            try:
                time_values = [time_val.decode('utf-8') for
                               time_val in time_values]
            except AttributeError:
                pass

            try:
                datetime_values = pd.to_datetime(time_values)
            except ValueError:
                # WRF DATETIME FORMAT
                datetime_values = \
                    pd.to_datetime(time_values,
                                   format="%Y-%m-%d_%H:%M:%S")

            self._obj[self.time_var].values = datetime_values

    @property
    def y_inverted(self):
        """Is the y-coord inverted"""
        if self._y_inverted is None:
            y_coords = self._obj[self.y_var].values
            if y_coords.ndim == 3:
                y_coords = y_coords[0]
            if y_coords.ndim == 2:
                self._y_inverted = (y_coords[-1, 0] > y_coords[0, 0])
            else:
                self._y_inverted = (y_coords[-1] > y_coords[0])
        return self._y_inverted

    @property
    def datetime(self):
        """Get datetime object for time variable"""
        self.to_datetime()
        return pd.to_datetime(self._obj[self.time_var].values)

    def _load_wrf_projection(self):
        """Load the osgeo.osr projection for WRF Grid.

        - 'MAP_PROJ': The map projection type as an integer.
        - 'TRUELAT1': True latitude 1.
        - 'TRUELAT2': True latitude 2.
        - 'MOAD_CEN_LAT': Mother of all domains center latitude.
        - 'STAND_LON': Standard longitude.
        - 'POLE_LAT': Pole latitude.
        - 'POLE_LON': Pole longitude.
        """
        # load in params from WRF Global Attributes
        possible_proj_params = ('MAP_PROJ', 'TRUELAT1', 'TRUELAT2',
                                'MOAD_CEN_LAT', 'STAND_LON', 'POLE_LAT',
                                'POLE_LON', 'CEN_LAT', 'CEN_LON', 'DX', 'DY')
        proj_params = dict()
        for proj_param in possible_proj_params:
            if proj_param in self._obj.attrs:
                proj_params[proj_param] = self._obj.attrs[proj_param]

        # determine projection from WRF Grid
        proj = wrf.projection.getproj(**proj_params)

        # export to Proj4 and add as osr projection
        self._projection = osr.SpatialReference()
        self._projection.ImportFromProj4(str(proj.proj4()))

    def _load_grib_projection(self):
        """Get the osgeo.osr projection for Grib Grid.
            - grid_type:  Lambert Conformal
            - Latin1:     True latitude 1.
            - Latin2:     True latitude 2.
            - Lov:        Central meridian.
            - Lo1:        Pole longitude.
            - La1:        Pole latitude.
            - Dx:         [ 3.]
            - Dy:         [ 3.]
        """
        lat_var_attrs = self._obj[self.y_var].attrs
        if 'Lambert Conformal' in lat_var_attrs['grid_type']:
            mean_lat = self._obj[self.y_var].mean().values
            proj4_str = ("+proj=lcc "
                         "+lat_1={true_lat_1} "
                         "+lat_2={true_lat_2} "
                         "+lat_0={latitude_of_origin} "
                         "+lon_0={central_meridian} "
                         "+x_0=0 +y_0=0 "
                         "+ellps=WGS84 +datum=WGS84 "
                         "+units=m +no_defs") \
                .format(true_lat_1=lat_var_attrs['Latin1'][0],
                        true_lat_2=lat_var_attrs['Latin2'][0],
                        latitude_of_origin=mean_lat,
                        central_meridian=lat_var_attrs['Lov'][0])
        else:
            raise ValueError("Unsupported projection: {grid_type}"
                             .format(grid_type=lat_var_attrs['grid_type']))

        # export to Proj4 and add as osr projection
        self._projection = osr.SpatialReference()
        self._projection.ImportFromProj4(proj4_str)

    @property
    def projection(self):
        """:func:`osgeo.osr.SpatialReference`
            The projection for the dataset.
        """
        if self._projection is None:
            # read projection information from global attributes
            map_proj4 = self._obj.attrs.get('proj4')
            if map_proj4 is not None:
                self._projection = osr.SpatialReference()
                self._projection.ImportFromProj4(str(map_proj4))
            elif 'MAP_PROJ' in self._obj.attrs:
                self._load_wrf_projection()
            elif 'grid_type' in self._obj[self.y_var].attrs:
                self._load_grib_projection()
            elif 'ProjectionCoordinateSystem' in self._obj.keys():
                # national water model
                proj4_str = self._obj['ProjectionCoordinateSystem'] \
                                .attrs['proj4']
                self._projection = osr.SpatialReference()
                self._projection.ImportFromProj4(str(proj4_str))
            else:
                # default to EPSG 4326
                self._projection = osr.SpatialReference()
                self._projection.ImportFromEPSG(4326)
            # make sure EPSG loaded if possible
            self._projection.AutoIdentifyEPSG()
        return self._projection

    @property
    def epsg(self):
        """str: EPSG code"""
        if self._epsg is None:
            self._epsg = self.projection.GetAuthorityCode(None)
        return self._epsg

    @property
    def dx(self):
        """float: Pixel size in x direction."""
        return self.geotransform[1]

    @property
    def dy(self):
        """float: Pixel size in y direction."""
        return -self.geotransform[-1]

    @property
    def geotransform(self):
        """:obj:`tuple`: The geotransform for grid."""
        if self._geotransform is None:
            if self._obj.attrs.get('geotransform') is not None:
                self._geotransform = [float(g) for g in
                                      self._obj.attrs.get('geotransform')]

            elif str(self.epsg) != '4326':
                proj_y, proj_x = self.coords
                self._geotransform = geotransform_from_yx(proj_y,
                                                          proj_x)
            else:
                self._geotransform = geotransform_from_yx(*self.latlon)

        return self._geotransform

    @property
    def affine(self):
        """:func:`Affine`: The affine for the transformation."""
        if self._affine is None:
            self._affine = Affine.from_gdal(*self.geotransform)
        return self._affine

    @property
    def x_size(self):
        """int: Number of columns in the dataset."""
        return self._obj.dims[self.x_dim]

    @property
    def y_size(self):
        """int: Number of rows in the dataset."""
        return self._obj.dims[self.y_dim]

    @property
    def _raw_coords(self):
        """Gets the raw coordinated of dataset"""
        x_coords = self._obj[self.x_var].values
        y_coords = self._obj[self.y_var].values

        if x_coords.ndim == 3:
            x_coords = x_coords[0]
        if y_coords.ndim == 3:
            y_coords = y_coords[0]

        if x_coords.ndim < 2:
            x_coords, y_coords = np.meshgrid(x_coords, y_coords)

        # WRF & NWM Grids are upside down
        if self.y_inverted:
            x_coords = x_coords[::-1]
            y_coords = y_coords[::-1]

        return y_coords, x_coords

    @property
    def latlon(self):
        """Returns lat,lon arrays

            .. warning:: The grids always be returned with [0,0]
                as Northeast and [-1,-1] as Southwest.
        """
        if 'MAP_PROJ' in self._obj.attrs:
            lat, lon = wrf.latlon_coords(self._obj, as_np=True)
            if lat.ndim == 3:
                lat = lat[0]
            if lon.ndim == 3:
                lon = lon[0]
            # WRF Grid is upside down
            lat = lat[::-1]
            lon = lon[::-1]
        else:
            lat, lon = self._raw_coords

        if self.coords_projected:
            lon, lat = transform(Proj(self.projection.ExportToProj4()),
                                 Proj(init='epsg:4326'),
                                 lon,
                                 lat)

        if self.lon_to_180:
            lon = (lon + 180) % 360 - 180  # convert [0, 360] to [-180, 180]

        return lat, lon

    @property
    def coords(self):
        """Returns y, x coordinate arrays

            .. warning:: The grids always be returned with [0,0]
                as Northeast and [-1,-1] as Southwest.
        """
        if not self.coords_projected:
            lat, lon = self.latlon
            x_coords, y_coords = \
                transform(Proj(init='epsg:4326'),
                          Proj(self.projection.ExportToProj4()),
                          lon,
                          lat)
            return y_coords, x_coords
        return self._raw_coords

    @property
    def center(self):
        """Return the geographic center point of this dataset."""
        if self._center is None:
            # we can use a cache on our accessor objects, because accessors
            # themselves are cached on instances that access them.
            lat, lon = self.latlon
            self._center = (float(np.nanmean(lon)), float(np.nanmean(lat)))
        return self._center

    def _export_dataset(self, variable, new_data, grid):
        """Export subset of dataset."""
        lats, lons = grid.latlon

        return xr.Dataset({variable: (['time', 'y', 'x'],
                                      new_data,
                                      self._obj[variable].attrs),
                           },
                          coords={'lat': (['y', 'x'],
                                          lats,
                                          self._obj[variable]
                                          .coords[self.y_var].attrs),
                                  'lon': (['y', 'x'],
                                          lons,
                                          self._obj[variable]
                                          .coords[self.x_var].attrs),
                                  'time': (['time'],
                                           self._obj[self.time_var].values,
                                           self._obj[self.time_var].attrs),
                                  },
                          attrs={'proj4': grid.proj4,
                                 'geotransform': grid.geotransform,
                                 }
                          )

    def resample(self, variable, match_grid):
        """Resample data to grid.

            Parameters
            ----------
            variable: :obj:`str`
                Name of variable in dataset.
            match_grid: :func:`gdal.Dataset` or :func:`sloot.grid.GDALGrid`
                Grid you want the data resampled to match resolution.
                You can also pass the path to the grid.
        """
        new_data = []
        for band in range(self._obj.dims[self.time_dim]):
            data = self._obj[variable][band].values
            arr_grid = ArrayGrid(in_array=data,
                                 wkt_projection=self.projection.ExportToWkt(),
                                 geotransform=self.geotransform)
            resampled_data_grid = resample_grid(original_grid=arr_grid,
                                                match_grid=match_grid,
                                                as_gdal_grid=True)
            new_data.append(resampled_data_grid.np_array())

        self.to_datetime()
        return self._export_dataset(variable, np.array(new_data),
                                    resampled_data_grid)

    def _getvar(self, variable, yslice, xslice):
        """Get the variable either directly or calculated"""
        # FAILED ATTEMPT TO USE wrf.getvar
        # if 'MAP_PROJ' in self._obj.attrs:
        #    try:
        #        nc_file = self._obj._file_obj.ds
        #    except AttributeError:
        #        nc_file = self._obj._file_obj.file_objs
        #    var = wrf.getvar(nc_file, variable)
        def extract_slice(var, slice_arr):
            """extract by slice"""
            if var.ndim == 4:
                var = var[slice_arr[0], slice_arr[1],
                          slice_arr[2], slice_arr[3]]
            if var.ndim == 3:
                var = var[slice_arr[0], slice_arr[1], slice_arr[2]]
            else:
                var = var[slice_arr[0], slice_arr[1]]
            return var

        var = self._obj[variable]
        slc = [slice(None)] * var.ndim
        # flip in y-direction
        if self.y_inverted:
            slc[var.get_axis_num(self.y_dim)] = slice(None, None, -1)
            var = extract_slice(var, slc)
        # get data out
        slc[var.get_axis_num(self.x_dim)] = xslice
        slc[var.get_axis_num(self.y_dim)] = yslice
        return extract_slice(var, slc)

    def getvar(self, variable,
               yslice=slice(None),
               xslice=slice(None),
               calc_4d_method=None,
               calc_4d_dim=None):
        """Get variable from model with subset options.

            .. warning:: The grids will always be returned with [0,0]
                as Northeast and [-1,-1] as Southwest.

            Parameters
            ----------
            variable: :obj:`str`
                Name of variable in dataset.
            yslice: :obj:`slice`, optional
                Slice in y-direction of grid to extract data from.
            xslice: :obj:`slice`, optional
                Slice in x-direction of grid to extract data from.
            calc_4d_method: :obj:`str`
                Method to convert 4D variables to 3D variables
                (Ex. 'mean', 'min', or 'max').
            calc_4d_dim: :obj:`str`
                Dimension to reduce grid from 4D to 3D (Ex. 'top_bottom').

            Returns
            -------
            :func:`xarray.DataArray`
        """
        data = self._getvar(variable, yslice, xslice)

        if data.ndim == 4:
            if calc_4d_method is None or calc_4d_dim is None:
                raise ValueError("The variable {var} has 4 dimension. "
                                 "Need 'calc_4d_method' and 'calc_4d_dim' "
                                 "to proceed ...".format(var=variable))
            data = getattr(data, calc_4d_method)(dim=calc_4d_dim)

        data[self.time_var] = self._obj[self.time_var]

        return data

    def to_projection(self, variable, projection):
        """Convert Grid to New Projection.

            Parameters
            ----------
            variable: :obj:`str`
                Name of variable in dataset.
            projection: :func:`osr.SpatialReference`
                Projection to convert data to.

            Returns
            -------
            :func:`xarray.Dataset`
        """
        new_data = []
        for band in range(self._obj.dims[self.time_dim]):
            arr_grid = ArrayGrid(in_array=self._obj[variable][band].values,
                                 wkt_projection=self.projection.ExportToWkt(),
                                 geotransform=self.geotransform)
            ggrid = arr_grid.to_projection(projection, gdalconst.GRA_Average)
            new_data.append(ggrid.np_array())

        self.to_datetime()
        return self._export_dataset(variable, np.array(new_data),
                                    ggrid)

    def to_utm(self, variable):
        """Convert Grid to UTM projection at center of grid.

            Parameters
            ----------
            variable: :obj:`str`
                Name of variable in dataset.

            Returns
            -------
            :func:`xarray.Dataset`
        """
        # get utm projection
        center_lon, center_lat = self.center
        dst_proj = utm_proj_from_latlon(center_lat, center_lon,
                                        as_osr=True)
        return self.to_projection(variable, dst_proj)

    def to_tif(self, variable, time_index, out_path):
        """Dump a variable at a time index to a geotiff.

            Parameters
            ----------
            variable: :obj:`str`
                Name of variable in dataset.
            time_index: int
                0-based time index,
            out_path: :obj:`str`
                Path to output geotiff file,
        """
        arr_grid = ArrayGrid(in_array=self._obj[variable][time_index].values,
                             wkt_projection=self.projection.ExportToWkt(),
                             geotransform=self.geotransform)
        arr_grid.to_tif(out_path)
