# -*- coding: utf-8 -*-
"""
   CreateInflowFileFromLDASRunoff.py
   RAPIDpy

   Created by Alan D. Snow, 2015
   Adapted from CreateInflowFileFromECMWFRunoff.py.
   License: BSD-3-Clause
"""
from netCDF4 import Dataset

from .CreateInflowFileFromGriddedRunoff import \
    CreateInflowFileFromGriddedRunoff


class CreateInflowFileFromLDASRunoff(CreateInflowFileFromGriddedRunoff):
    """Create Inflow File From LDAS Runoff

    Base class for creating RAPID NetCDF input
    of water inflow based on LDAS land surface model
    runoff and previously created weight table.
    """
    land_surface_model_name = "LDAS"

    def __init__(self,
                 lat_dim,  # "g0_lat_0",
                 lon_dim,  # "g0_lon_1",
                 lat_var,  # "g0_lat_0",
                 lon_var,  # "g0_lon_1",
                 runoff_vars):  # ["Qsb_GDS0_SFC_ave1h", "Qs_GDS0_SFC_ave1h"],
        """Define the attributes to look for"""
        self.dims_oi = [lon_dim, lat_dim]
        self.vars_oi = [lon_var, lat_var] + runoff_vars
        self.runoff_vars = runoff_vars
        self.length_time = {"Hourly": 1}

        super(CreateInflowFileFromLDASRunoff, self).__init__()

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
