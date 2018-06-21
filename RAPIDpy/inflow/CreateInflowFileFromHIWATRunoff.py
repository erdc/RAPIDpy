# -*- coding: utf-8 -*-
"""
   CreateInflowFileFromHIWAT.py
   RAPIDpy

   Created by Xiaohui Qiao, 2018
   Adapted from CreateInflowFileFromECMWFRunoff.py.
   License: BSD-3-Clause
"""
from netCDF4 import Dataset

from .CreateInflowFileFromGriddedRunoff import \
    CreateInflowFileFromGriddedRunoff


class CreateInflowFileFromHIWATRunoff(CreateInflowFileFromGriddedRunoff):
    """Create Inflow File From HIWAT precipitation

    Creates RAPID NetCDF input of water inflow based on
    HIWAT precipitation and previously created weight table.
    """
    land_surface_model_name = "HIWAT"
    header_wt = ['rivid', 'area_sqm', 'lon_index', 'lat_index', 'npoints']
    dims_oi = ['longitude', 'latitude', 'time']
    vars_oi = ['longitude', 'latitude', 'time', 'PCP']
    length_time = {"Hourly": 1}

    def __init__(self):
        """Define the attributes to look for"""
        self.runoff_vars = ['PCP']
        super(CreateInflowFileFromHIWATRunoff, self).__init__()

    def data_validation(self, in_nc):
        """Check the necessary dimensions and variables in the input
        netcdf data"""
        data_nc = Dataset(in_nc)

        for dim in self.dims_oi:
            if dim not in data_nc.dimensions.keys():
                data_nc.close()
                raise Exception(self.error_messages[1])

        for var in self.vars_oi:
            if var not in data_nc.variables.keys():
                data_nc.close()
                raise Exception(self.error_messages[2])
        data_nc.close()
