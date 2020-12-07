# -*- coding: utf-8 -*-
"""
   CreateInflowFileFromERA5Runoff.py
   RAPIDpy

   Adapted from CreateInflowFileFromERAInterimRunoff.py.
   License: BSD-3-Clause
"""
from netCDF4 import Dataset

from .CreateInflowFileFromGriddedRunoff import \
    CreateInflowFileFromGriddedRunoff


class CreateInflowFileFromERA5Runoff(CreateInflowFileFromGriddedRunoff):
    """Create Inflow File From ERA Interim Runoff

    Creates RAPID NetCDF input of water inflow based on
    ERA5 runoff and previously created weight table.
    """
    land_surface_model_name = "ERA5"
    header_wt = ['rivid', 'area_sqm', 'lon_index', 'lat_index', 'npoints']
    dims_oi = ['time', 'longitude', 'latitude']
    vars_oi = ['time', 'longitude', 'latitude', 'ro']
    runoff_vars = ['ro']
    length_time = {"Daily": 1, "3-Hourly": 8}

    def __init__(self):
        """Define the attributes to look for"""
        self.runoff_vars = ['ro']
        super(CreateInflowFileFromERA5Runoff, self).__init__()

    def data_validation(self, in_nc):
        """Check the necessary dimensions and variables in the input
        netcdf data"""
        data_nc = Dataset(in_nc)

        dims = list(data_nc.dimensions)
        ndim = len(dims)

	ndim_intersect = len(set(dims) & set(self.dims_oi))
        if not ndim_intersect == ndim:
            data_nc.close()
            raise Exception("{0} {1}".format(self.error_messages[1], dims))

        nc_vars = list(data_nc.variables)
        nvar = len(nc_vars)
 
        nvar_intersect = len(set(nc_vars) & set(self.vars_oi))
        if nvar_intersect == nvar:
            runoff_intersect = set(nc_vars) & set(self.runoff_vars)
            self.runoff_vars = [runoff_intersect.pop()]
        else:
            data_nc.close()
            raise Exception("{0} {1}".format(self.error_messages[2], nc_vars))
        data_nc.close()
