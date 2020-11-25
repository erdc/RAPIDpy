# -*- coding: utf-8 -*-
"""
   CreateInflowFileFromERARunoff.py
   RAPIDpy

   Created by Alan D. Snow, 2015
   Adapted from CreateInflowFileFromECMWFRunoff.py.
   License: BSD-3-Clause
"""
from netCDF4 import Dataset

from .CreateInflowFileFromGriddedRunoff import \
    CreateInflowFileFromGriddedRunoff


class CreateInflowFileFromERARunoff(CreateInflowFileFromGriddedRunoff):
    """Create Inflow File From ERA Runoff

    Creates RAPID NetCDF input of water inflow based on
    ERA runoff and previously created weight table.
    """
    land_surface_model_name = "ERA"
    header_wt = ['rivid', 'area_sqm', 'lon_index', 'lat_index', 'npoints']
    dims_oi = [['lon', 'lat', 'time'],
               ['longitude', 'latitude', 'time'],
               ['time','lon','lat'],
               ['time','longitude','latitude'],
               ['latitude','longitude','time'],
               ['time','latitude','longitude']]
    vars_oi = [['lon', 'lat', 'time', 'RO'],['lon','lat','time','ro'],
               ['time','lon','lat','RO'],['time','lon','lat','ro'],
               ['time','longitude','latitude','RO'],['time','longitude','latitude','ro'],
               ['longitude', 'latitude', 'time', 'RO'],['longitude', 'latitude', 'time', 'ro'],
               ['latitude','longitude','time','RO'],['latitude','longitude','time','ro'],
               ['latitude','longitude','RO','time'],['latitude','longitude','ro','time'],
               ['time','latitude','longitude','RO'],['time','latitude','longitude','ro']]
    length_time = {"Daily": 1, "3-Hourly": 8, "1-Hourly":24}

    def __init__(self):
        """Define the attributes to look for"""
        self.runoff_vars = ['ro']
        super(CreateInflowFileFromERARunoff, self).__init__()

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
        elif nc_vars == self.vars_oi[2]:
            self.runoff_vars = [self.vars_oi[2][-1]]
        elif nc_vars == self.vars_oi[3]:
            self.runoff_vars = [self.vars_oi[3][-1]]
        elif nc_vars == self.vars_oi[4]:
            self.runoff_vars = [self.vars_oi[4][-1]]
        elif nc_vars == self.vars_oi[5]:
            self.runoff_vars = [self.vars_oi[5][-1]]
        elif nc_vars == self.vars_oi[6]:
            self.runoff_vars = [self.vars_oi[6][-1]]
        elif nc_vars == self.vars_oi[7]:
            self.runoff_vars = [self.vars_oi[7][-1]]
        elif nc_vars == self.vars_oi[8]:
            self.runoff_vars = [self.vars_oi[8][-1]]
        elif nc_vars == self.vars_oi[9]:
            self.runoff_vars = [self.vars_oi[9][-1]]
        elif nc_vars == self.vars_oi[10]:
            self.runoff_vars = [self.vars_oi[10][-2]]
        elif nc_vars == self.vars_oi[11]:
            self.runoff_vars = [self.vars_oi[11][-2]]
        elif nc_vars == self.vars_oi[12]:
            self.runoff_vars = [self.vars_oi[12][-1]]
        elif nc_vars == self.vars_oi[13]:
            self.runoff_vars = [self.vars_oi[13][-1]]
        else:
            data_nc.close()
            raise Exception("{0} {1}".format(self.error_messages[2], nc_vars))
        data_nc.close()
