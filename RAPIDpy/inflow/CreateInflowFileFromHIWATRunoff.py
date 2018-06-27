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
import os


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

def make_hiwat_fake_data(in_nc_path):

    print("creating fake zero data for 2 days")
    lsm_file_list = []
    for walkdir_info in os.walk(in_nc_path,
                                followlinks=True):
        for lsm_file in walkdir_info[2]:
            if lsm_file.endswith(".nc") or \
                    lsm_file.endswith(".nc4"):
                lsm_file_list.append(
                    os.path.join(walkdir_info[0], lsm_file))
    lsm_file_list = sorted(lsm_file_list)

    if len(lsm_file_list) > 1:
        raise Exception("More than one hiwat input file found.")
    in_nc_file_path = lsm_file_list[0]
    file_folder_path, filename_full = os.path.split(lsm_file_list[0])
    filename, ext = os.path.splitext(filename_full)
    out_nc_file_path = os.path.join(file_folder_path, filename + "_fake" + ext)
    print(out_nc_file_path)

    from shutil import copyfile
    copyfile(in_nc_file_path, out_nc_file_path)

    with Dataset(out_nc_file_path, "a") as out_nc:

         out_nc.variables['time'][: ] += 24*3600*2
         out_nc.variables['PCP'][:,:,:] = 0.0

