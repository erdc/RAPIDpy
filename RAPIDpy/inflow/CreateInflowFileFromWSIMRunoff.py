# -*- coding: utf-8 -*-
"""
   CreateInflowFileFromLDASRunoff.py
   RAPIDpy
   Created by Matthew P. Geheran, 2018.
   Adapted from CreateInflowFileFromLDASRunoff.py.
"""
from netCDF4 import Dataset

from .CreateInflowFileFromGriddedRunoff import \
    CreateInflowFileFromGriddedRunoff


class CreateInflowFileFromWSIMRunoff(CreateInflowFileFromGriddedRunoff):
    """Create Inflow File From WSIM Runoff

    Base class for creating RAPID NetCDF input
    of water inflow based on the WSIM land surface model
    runoff and previously created weight table.
    """
    land_surface_model_name = "WSIM"

    def __init__(self,
                 lat_dim,  # "lat",
                 lon_dim,  # "lon",
                 lat_var,  # "latitude",
                 lon_var,  # "longitude",
                 runoff_vars):  # ["surface_runoff", "sub_surface_runoff"],
        """Define the attributes to look for"""
        self.dims_oi = [lon_dim, lat_dim]
        self.vars_oi = [lon_var, lat_var] + runoff_vars
        self.runoff_vars = runoff_vars
        self.length_time = {"Hourly": 6}

        super(CreateInflowFileFromWSIMRunoff, self).__init__()

    def data_validation(self, in_nc):
        """Check the necessary dimensions and variables in the
        input netcdf data"""
        data_nc = Dataset(in_nc)
        for dim in self.dims_oi:
            # MPG DEBUG:
            print dim
            print data_nc.dimensions.keys()
            if dim not in data_nc.dimensions.keys():
                data_nc.close()
                raise Exception(self.error_messages[1])

        for var in self.vars_oi:
            # MPG DEBUG:
            print var
            print data_nc.variables.keys()
            if var not in data_nc.variables.keys():
                data_nc.close()
                raise Exception(self.error_messages[2])

        data_nc.close()
        return
