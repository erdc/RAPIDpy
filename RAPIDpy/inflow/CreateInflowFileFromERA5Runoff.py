# -*- coding: utf-8 -*-
"""
   CreateInflowFileFromERA5Runoff.py
   RAPIDpy

   Created originally by Alan D. Snow, 2015
   Adapted for ERA-5 from CreateInflowFileFromERAInterimRunoff.py
        by Chris Edwards, 2018
   License: BSD-3-Clause
"""
from netCDF4 import Dataset

from .CreateInflowFileFromGriddedRunoff import \
    CreateInflowFileFromGriddedRunoff


class CreateInflowFileFromERA5Runoff(CreateInflowFileFromGriddedRunoff):
    """Create Inflow File From ERA 5 Runoff

    Creates RAPID NetCDF input of water inflow based on
    ERA 5 runoff and previously created weight table.
    """
    land_surface_model_name = "ERA 5"
    header_wt = ['rivid', 'area_sqm', 'lon_index', 'lat_index', 'npoints']
    dims_oi = [['lon', 'lat', 'time'], ['longitude', 'latitude', 'time']]
    vars_oi = [["lon", "lat", "time", "RO"],
               ['longitude', 'latitude', 'time', 'ro']]
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
