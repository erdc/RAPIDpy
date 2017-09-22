# -*- coding: utf-8 -*-
"""
   CreateInflowFileFromERAInterimRunoff.py
   RAPIDpy

   Created by Alan D. Snow, 2015
   Adapted from CreateInflowFileFromECMWFRunoff.py.
   License: BSD-3-Clause
"""
from netCDF4 import Dataset

from .CreateInflowFileFromGriddedRunoff import \
    CreateInflowFileFromGriddedRunoff


class CreateInflowFileFromERAInterimRunoff(CreateInflowFileFromGriddedRunoff):
    """Create Inflow File From ERA Interim Runoff

    Creates RAPID NetCDF input of water inflow based on
    ERA Interim runoff and previously created weight table.
    """
    land_surface_model_name = "ERA Interim"
    header_wt = ['rivid', 'area_sqm', 'lon_index', 'lat_index', 'npoints']
    dims_oi = [['lon', 'lat', 'time'], ['longitude', 'latitude', 'time']]
    vars_oi = [["lon", "lat", "time", "RO"],
               ['longitude', 'latitude', 'time', 'ro']]
    length_time = {"Daily": 1, "3-Hourly": 8}

    def __init__(self):
        """Define the attributes to look for"""
        self.runoff_vars = ['ro']
        super(CreateInflowFileFromERAInterimRunoff, self).__init__()

    def data_validation(self, in_nc):
        """Check the necessary dimensions and variables in the input
        netcdf data"""
        data_nc = Dataset(in_nc)

        dims = list(data_nc.dimensions)

        if dims not in self.dims_oi:
            data_nc.close()
            raise Exception("{0} {1}".format(self.error_messages[1], dims))

        nc_vars = list(data_nc.variables)

        if nc_vars == self.vars_oi[0]:
            self.runoff_vars = [self.vars_oi[0][-1]]
        elif nc_vars == self.vars_oi[1]:
            self.runoff_vars = [self.vars_oi[1][-1]]
        else:
            data_nc.close()
            raise Exception("{0} {1}".format(self.error_messages[2], nc_vars))
        data_nc.close()
