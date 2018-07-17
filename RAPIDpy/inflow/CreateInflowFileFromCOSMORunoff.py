# -*- coding: utf-8 -*-
"""
   CreateInflowFileFromCOSMORunoff.py
   RAPIDpy

   Created by Xiaohui Qiao, 2018
   Adapted from CreateInflowFileFromLDASRunoff.py.
   License: BSD-3-Clause
"""
from netCDF4 import Dataset

from .CreateInflowFileFromGriddedRunoff import \
    CreateInflowFileFromGriddedRunoff


class CreateInflowFileFromCOSMORunoff(CreateInflowFileFromGriddedRunoff):
    """Create Inflow File From COSMO Runoff

    Base class for creating RAPID NetCDF input
    of water inflow based on COSMO model
    runoff and previously created weight table.
    """
    land_surface_model_name = "COSMO"

    def __init__(self,
                 lat_dim,  # "rlat",
                 lon_dim,  # "rlon",
                 time_dim,  # "time"
                 lat_var,  # "lat",
                 lon_var,  # "lon",
                 runoff_vars):  # ["RUNOFF_G", "RUNOFF_S"],
        """Define the attributes to look for"""
        self.dims_oi = [lon_dim, lat_dim, time_dim]
        self.vars_oi = [lon_var, lat_var] + runoff_vars
        self.runoff_vars = runoff_vars
        self.length_time = {"Hourly": 1}

        super(CreateInflowFileFromCOSMORunoff, self).__init__()

    def data_validation(self, in_nc):
        """Check the necessary dimensions and variables in the
        input netcdf data"""
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
        return
