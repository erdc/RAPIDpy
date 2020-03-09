# -*- coding: utf-8 -*-
"""
   CreateInflowFileFromGALWEMRunoff.py
   RAPIDpy

   Created by Matthew P. Geheran, 2018
   Adapted from CreateInflowFileFromECMWFRunoff.py.
"""
from netCDF4 import Dataset

from .CreateInflowFileFromGriddedRunoff import \
    CreateInflowFileFromGriddedRunoff


class CreateInflowFileFromGALWEMRunoff(CreateInflowFileFromGriddedRunoff):
    """Create Inflow File From LDAS Runoff

    Base class for creating RAPID NetCDF input
    of water inflow based on LDAS land surface model
    runoff and previously created weight table.
    """
    land_surface_model_name = "LDAS"

    def __init__(self,
                 lat_dim,  # "lat",
                 lon_dim,  # "lon",
                 lat_var,  # "lat",
                 lon_var,  # "lon",
                 runoff_vars):  # ["ssrun", "bgrun"],
        """Define the attributes to look for"""
        self.dims_oi = [lon_dim, lat_dim]
        self.vars_oi = [lon_var, lat_var] + runoff_vars
        self.runoff_vars = runoff_vars
        self.length_time = {"3-Hourly": 1}

        super(CreateInflowFileFromGALWEMRunoff, self).__init__()

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
