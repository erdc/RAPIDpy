# -*- coding: utf-8 -*-
#
#  read.py
#  pangaea
#
#  Author : Alan D Snow, 2017.
#  License: BSD 3-Clause
"""pangea.read

    This module provides helper functions to read in
    land surface model datasets.
"""
import numpy as np
import pandas as pd
import xarray as xr


def open_mfdataset(path_to_lsm_files,
                   lat_var,
                   lon_var,
                   time_var,
                   lat_dim,
                   lon_dim,
                   time_dim,
                   lon_to_180=False,
                   coords_projected=False,
                   loader=None,
                   engine=None,
                   autoclose=True):
    """
    Wrapper to open land surface model netcdf files
    using :func:`xarray.open_mfdataset`.

    .. warning:: The time dimension and variable will both be
        renamed to 'time' to enable slicing.

    Parameters
    ----------
    path_to_lsm_files: :obj:`str`
        Path to land surface model files with wildcard.
        (Ex. '/path/to/files/*.nc')
    lat_var: :obj:`str`
        Latitude variable (Ex. lat).
    lon_var: :obj:`str`
        Longitude variable (Ex. lon).
    time_var: :obj:`str`
        Time variable (Ex. time).
    lat_dim: :obj:`str`
        Latitude dimension (Ex. lat).
    lon_dim: :obj:`str`
        Longitude dimension (Ex. lon).
    time_dim: :obj:`str`
        Time dimension (ex. time).
    lon_to_180: bool, optional, default=False
        It True, will convert longitude from [0 to 360]
        to [-180 to 180].
    coords_projected: bool, optional, default=False
        It True, it will assume the coordinates are already
        in the projected coordinate system.
    loader: str, optional, default=None
        If 'hrrr', it will load in the HRRR dataset.
    engine: str, optional
        See: :func:`xarray.open_mfdataset` documentation.
    autoclose: :obj:`str`, optional, default=True
        If True, will use autoclose option with
        :func:`xarray.open_mfdataset`.

    Returns
    -------
    :func:`xarray.Dataset`


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
    """
    def define_coords(xds):
        """xarray loader to ensure coordinates are loaded correctly"""
        # remove time dimension from lat, lon coordinates
        if xds[lat_var].ndim == 3:
            xds[lat_var] = xds[lat_var].squeeze(time_dim)
        # make sure coords are defined as coords
        # prior version was xds.set_coords, changed with deprecation of inplace routine
        if lat_var not in xds.coords \
                or lon_var not in xds.coords \
                or time_var not in xds.coords:
            xds = xds.set_coords([lat_var, lon_var, time_var])
        #                   inplace=True)
        return xds

    def extract_hrrr_date(xds):
        """xarray loader for HRRR"""
        for var in xds.variables:
            if 'initial_time' in xds[var].attrs.keys():
                grid_time = pd.to_datetime(xds[var].attrs['initial_time'],
                                           format="%m/%d/%Y (%H:%M)")
                if 'forecast_time' in xds[var].attrs.keys():
                    time_units = 'h'
                    if 'forecast_time_units' in xds[var].attrs.keys():
                        time_units = \
                            str(xds[var].attrs['forecast_time_units'][0])
                    time_dt = int(xds[var].attrs['forecast_time'][0])
                    grid_time += np.timedelta64(time_dt, time_units)

                return xds.assign(time=grid_time)
        return xds

    if loader == 'hrrr':
        preprocess = extract_hrrr_date
        engine = 'pynio' if engine is None else engine
    else:
        preprocess = define_coords

    xds = xr.open_mfdataset(path_to_lsm_files,
                            autoclose=autoclose,
                            preprocess=preprocess,
                            concat_dim=time_dim,
                            engine=engine,
                            combine = 'nested') 

    xds.lsm.y_var = lat_var
    xds.lsm.x_var = lon_var
    xds.lsm.y_dim = lat_dim
    xds.lsm.x_dim = lon_dim
    xds.lsm.lon_to_180 = lon_to_180
    xds.lsm.coords_projected = coords_projected

    # make sure time dimensions are same for slicing
    # prior version only xds.rename(...
    xds = xds.rename(
        {
            time_dim: 'time',
            time_var: 'time',
        }
    #    inplace=True
    )

    xds.lsm.to_datetime()
    return xds
