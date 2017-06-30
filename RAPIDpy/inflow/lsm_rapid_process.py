# -*- coding: utf-8 -*-
#
#  lsm_rapid_process.py
#  RAPIDpy
#
#  Created by Alan D. Snow.
#  Copyright © 2015-2016 Alan D Snow. All rights reserved.
#  License: BSD 3-Clause

from datetime import datetime, timedelta
import multiprocessing
import os
import re
import traceback

# external packages
import pandas as pd
import pangaea
from netCDF4 import Dataset
import numpy as np

# local imports
from ..rapid import RAPID
from .CreateInflowFileFromERAInterimRunoff import CreateInflowFileFromERAInterimRunoff
from .CreateInflowFileFromLDASRunoff import CreateInflowFileFromLDASRunoff
from .CreateInflowFileFromWRFHydroRunoff import CreateInflowFileFromWRFHydroRunoff
from ..postprocess.generate_return_periods import generate_return_periods
from ..postprocess.generate_seasonal_averages import generate_seasonal_averages
from ..utilities import (case_insensitive_file_search,
                         get_valid_directory_list,
                         partition)


# ------------------------------------------------------------------------------
# MULTIPROCESSING FUNCTION
# ------------------------------------------------------------------------------
def generate_inflows_from_runoff(args):
    """
    prepare runoff inflow file for rapid
    """
    runoff_file_list = args[0]
    file_index_list = args[1]
    weight_table_file = args[2]
    grid_type = args[3]
    rapid_inflow_file = args[4]
    RAPID_Inflow_Tool = args[5]
    mp_lock = args[6]

    time_start_all = datetime.utcnow()

    if not isinstance(runoff_file_list, list):
        runoff_file_list = [runoff_file_list]
    else:
        runoff_file_list = runoff_file_list

    if not isinstance(file_index_list, list):
        file_index_list = [file_index_list]
    else:
        file_index_list = file_index_list
    if runoff_file_list and file_index_list:
        # prepare ECMWF file for RAPID
        index_string = "Index: {0}".format(file_index_list[0])
        if len(file_index_list) > 1:
            index_string += " to {0}".format(file_index_list[-1])
        print(index_string)
        runoff_string = "File(s): {0}".format(runoff_file_list[0])
        if len(runoff_file_list) > 1:
            runoff_string += " to {0}".format(runoff_file_list[-1])
        print(runoff_string)
        print("Converting inflow ...")
        try:
            RAPID_Inflow_Tool.execute(nc_file_list=runoff_file_list,
                                      index_list=file_index_list,
                                      in_weight_table=weight_table_file,
                                      out_nc=rapid_inflow_file,
                                      grid_type=grid_type,
                                      mp_lock=mp_lock,
                                      )
        except Exception:
            # This prints the type, value, and stack trace of the
            # current exception being handled.
            traceback.print_exc()
            raise

        time_finish_ecmwf = datetime.utcnow()
        print("Time to convert inflows: {0}".format(time_finish_ecmwf-time_start_all))

# ------------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------------------------

DEFAULT_LSM_INPUTS = {
    't255': {
        'file_datetime_re_pattern': r'\d{8}',
        'file_datetime_pattern': "%Y%m%d",
    },
    't511': {
        'file_datetime_re_pattern': r'\d{8}',
        'file_datetime_pattern': "%Y%m%d",
    },
    't159': {
        'file_datetime_re_pattern': r'\d{8}',
        'file_datetime_pattern': "%Y%m%d",
    },
    'gldas2': {
        'file_datetime_re_pattern': r'\d{8}\.\d{2}',
        'file_datetime_pattern': "%Y%m%d.%H",
    },
    'gldas': {
        'file_datetime_re_pattern': r'\d{8}\.\d{2}',
        'file_datetime_pattern': "%Y%m%d.%H",
    },
    'nldas': {
        'file_datetime_re_pattern': r'\d{8}\.\d{2}',
        'file_datetime_pattern': "%Y%m%d.%H",
    },
    'cmip5': {
        'file_datetime_re_pattern': r'\d{4}',
        'file_datetime_pattern': "%Y",
    },
    'lis': {
        'file_datetime_re_pattern': r'\d{10}',
        'file_datetime_pattern': "%Y%m%d%H",
    },
    'joules': {
        'file_datetime_re_pattern': r'\d{10}',
        'file_datetime_pattern': "%Y%m%d%H",
    },
    'wrf': {
        'file_datetime_re_pattern': r'\d{10}',
        'file_datetime_pattern': "%Y%m%d%H",
    },
}


def identify_lsm_grid(lsm_grid_path):
    """
    This is used to idenfity the input LSM grid
    """
    # check to see what kind of file we are dealing with
    lsm_example_file = Dataset(lsm_grid_path)

    # INDENTIFY LAT/LON DIMENSIONS
    dim_list = lsm_example_file.dimensions.keys()

    latitude_dim = "lat"
    if 'latitude' in dim_list:
        latitude_dim = 'latitude'
    elif 'g0_lat_0' in dim_list:
        # GLDAS/NLDAS MOSAIC
        latitude_dim = 'g0_lat_0'
    elif 'lat_110' in dim_list:
        # NLDAS NOAH/VIC
        latitude_dim = 'lat_110'
    elif 'north_south' in dim_list:
        # LIS/Joules
        latitude_dim = 'north_south'
    elif 'south_north' in dim_list:
        # WRF Hydro
        latitude_dim = 'south_north'
    elif 'Y' in dim_list:
        # FLDAS
        latitude_dim = 'Y'

    longitude_dim = "lon"
    if 'longitude' in dim_list:
        longitude_dim = 'longitude'
    elif 'g0_lon_1' in dim_list:
        # GLDAS/NLDAS MOSAIC
        longitude_dim = 'g0_lon_1'
    elif 'lon_110' in dim_list:
        # NLDAS NOAH/VIC
        longitude_dim = 'lon_110'
    elif 'east_west' in dim_list:
        # LIS/Joules
        longitude_dim = 'east_west'
    elif 'west_east' in dim_list:
        # WRF Hydro
        longitude_dim = 'west_east'
    elif 'X' in dim_list:
        # FLDAS
        longitude_dim = 'X'

    time_dim = None
    if 'time' in dim_list:
        time_dim = 'time'
    elif 'Time' in dim_list:
        time_dim = 'Time'
    elif 'Times' in dim_list:
        time_dim = 'Times'
    elif 'times' in dim_list:
        time_dim = 'times'


    lat_dim_size = len(lsm_example_file.dimensions[latitude_dim])
    lon_dim_size = len(lsm_example_file.dimensions[longitude_dim])

    # IDENTIFY VARIABLES
    var_list = lsm_example_file.variables.keys()

    latitude_var = "lat"
    if 'latitude' in var_list:
        latitude_var = 'latitude'
    elif 'g0_lat_0' in var_list:
        latitude_var = 'g0_lat_0'
    elif 'lat_110' in var_list:
        latitude_var = 'lat_110'
    elif 'north_south' in var_list:
        latitude_var = 'north_south'
    elif 'XLAT' in var_list:
        # WRF
        latitude_var = 'XLAT'
    elif 'Y' in var_list:
        # FLDAS
        latitude_var = 'Y'

    longitude_var = "lon"
    if 'longitude' in var_list:
        longitude_var = 'longitude'
    elif 'g0_lon_1' in var_list:
        longitude_var = 'g0_lon_1'
    elif 'lon_110' in var_list:
        longitude_var = 'lon_110'
    elif 'east_west' in var_list:
        longitude_var = 'east_west'
    elif 'XLONG' in var_list:
        # WRF
        longitude_var = 'XLONG'
    elif 'X' in var_list:
        # FLDAS
        longitude_var = 'X'

    time_var = None
    if 'time' in var_list:
        time_var = 'time'
    elif 'Time' in var_list:
        time_var = 'Time'
    elif 'Times' in var_list:
        time_var = 'Times'
    elif 'times' in var_list:
        time_var = 'times'

    surface_runoff_var = ""
    subsurface_runoff_var = ""
    total_runoff_var = ""
    for var in var_list:
        if var.startswith("SSRUN"):
            # NLDAS/GLDAS
            surface_runoff_var = var
        elif var.startswith("BGRUN"):
            # NLDAS/GLDAS
            subsurface_runoff_var = var
        elif var == "Qs_acc":
            # GLDAS v2
            surface_runoff_var = var
        elif var == "Qsb_acc":
            # GLDAS v2
            subsurface_runoff_var = var
        elif var == "Qs_tavg":
            # FLDAS
            surface_runoff_var = var
        elif var == "Qsb_tavg":
            # FLDAS
            subsurface_runoff_var = var
        elif var == "Qs_inst":
            # LIS
            surface_runoff_var = var
        elif var == "Qsb_inst":
            # LIS
            subsurface_runoff_var = var
        elif var == "SFROFF":
            # WRF Hydro
            surface_runoff_var = var
        elif var == "UDROFF":
            # WRF Hydro
            subsurface_runoff_var = var
        elif var.lower() == "ro":
            # ERA Interim
            total_runoff_var = var
        elif var == "total runoff":
            # CMIP5 data
            total_runoff_var = var

    # IDENTIFY GRID TYPE
    lsm_file_data = {
        "weight_file_name": "",
        "grid_type": "",
        "model_name": "",
        "description": "",
        "rapid_inflow_tool": None,
        "latitude_var": latitude_var,
        "longitude_var": longitude_var,
        "time_var": time_var,
        "latitude_dim": latitude_dim,
        "longitude_dim": longitude_dim,
        "time_dim": time_dim,
    }

    institution = ""
    title = ""
    try:
        institution = lsm_example_file.getncattr("institution")
    except AttributeError:
        pass
    try:
        title = lsm_example_file.getncattr("title")
    except AttributeError:
        pass

    runoff_vars = [surface_runoff_var, subsurface_runoff_var]

    if institution == "European Centre for Medium-Range Weather Forecasts" \
            or total_runoff_var.lower() == "ro":
        # these are the ECMWF models
        if lat_dim_size == 361 and lon_dim_size == 720:
            print("Runoff file identified as ERA Interim Low Res (T255) GRID")
            # A) ERA Interim Low Res (T255)
            # Downloaded as 0.5 degree grid
            #  dimensions:
            # 	 longitude = 720 ;
            # 	 latitude = 361 ;
            lsm_file_data["description"] = "ERA Interim (T255 Grid)"
            lsm_file_data["model_name"] = "erai"
            lsm_file_data["weight_file_name"] = r'weight_era_t255\.csv'
            lsm_file_data["grid_type"] = 't255'

        elif lat_dim_size == 512 and lon_dim_size == 1024:
            print("Runoff file identified as ERA Interim High Res (T511) GRID")
            # B) ERA Interim High Res (T511)
            #  dimensions:
            #   lon = 1024 ;
            #   lat = 512 ;
            lsm_file_data["description"] = "ERA Interim (T511 Grid)"
            lsm_file_data["weight_file_name"]= r'weight_era_t511\.csv'
            lsm_file_data["model_name"] = "erai"
            lsm_file_data["grid_type"] = 't511'
        elif lat_dim_size == 161 and lon_dim_size == 320:
            print("Runoff file identified as ERA 20CM (T159) GRID")
            # C) ERA 20CM (T159) - 3hr - 10 ensembles
            # Downloaded as 1.125 degree grid
            #  dimensions:
            #   longitude = 320 ;
            #   latitude = 161 ;
            lsm_file_data["description"] = "ERA 20CM (T159 Grid)"
            lsm_file_data["weight_file_name"] = r'weight_era_t159\.csv'
            lsm_file_data["model_name"] = "era_20cm"
            lsm_file_data["grid_type"] = 't159'
        else:
            lsm_example_file.close()
            raise Exception("Unsupported ECMWF grid.")

        lsm_file_data["rapid_inflow_tool"] = \
            CreateInflowFileFromERAInterimRunoff()

    elif institution == "NASA GSFC":
        if title == "GLDAS2.0 LIS land surface model output":
            print("Runoff file identified as GLDAS v2 LIS GRID")
            # this is the LIS model
            lsm_file_data["weight_file_name"] = r'weight_gldas2\.csv'
            lsm_file_data["grid_type"] = 'gldas2'
            lsm_file_data["description"] = "GLDAS2.0 LIS"
            lsm_file_data["model_name"] = "nasa"

        else:
            print("Runoff file identified as LIS GRID")
            # this is the LIS model (can be FLDAS)
            # THIS CASE CAN ALSO BE FOR FLDAS, however you will need to add
            # the file_datetime_pattern && file_datetime_re_pattern for it to
            # work if it is not 3-hourly time step.
            lsm_file_data["weight_file_name"] = r'weight_lis\.csv'
            lsm_file_data["grid_type"] = 'lis'
            lsm_file_data["description"] = "NASA GSFC LIS"
            lsm_file_data["model_name"] = "nasa"

    elif institution == "Met Office, UK":
        print("Runoff file identified as Joules GRID")
        lsm_file_data["weight_file_name"] = r'weight_joules\.csv'
        lsm_file_data["grid_type"] = 'joules'
        lsm_file_data["description"] = "Met Office Joules"
        lsm_file_data["model_name"] = "met_office"

    elif institution == "NCAR, USACE, USBR":
        print("Runoff file identified as CMIP5")
        lsm_file_data["weight_file_name"] = r'weight_cmip5\.csv'
        lsm_file_data["grid_type"] = 'cmip5'
        lsm_file_data["description"] = "CMIP5 Runoff"
        lsm_file_data["model_name"] = "cmip5"

        runoff_vars = [total_runoff_var]

    elif surface_runoff_var.startswith("SSRUN") \
            and subsurface_runoff_var.startswith("BGRUN"):

        lsm_file_data["model_name"] = "nasa"
        if lat_dim_size == 600 and lon_dim_size == 1440:
            print("Runoff file identified as GLDAS GRID")
            # GLDAS NC FILE
            # dimensions:
            #     g0_lat_0 = 600 ;
            #     g0_lon_1 = 1440 ;
            # variables
            # SSRUN_GDS0_SFC_ave1h (surface), BGRUN_GDS0_SFC_ave1h (subsurface)
            #  or
            # SSRUNsfc_GDS0_SFC_ave1h (surface), BGRUNsfc_GDS0_SFC_ave1h (subsurface)
            lsm_file_data["description"] = "GLDAS"
            lsm_file_data["weight_file_name"] = r'weight_gldas\.csv'
            lsm_file_data["grid_type"] = 'gldas'

        elif lat_dim_size <= 224 and lon_dim_size <= 464:
            print("Runoff file identified as NLDAS GRID")
            # NLDAS MOSAIC FILE
            # dimensions:
            #     g0_lat_0 = 224 ;
            #     g0_lon_1 = 464 ;
            # NLDAS NOAH/VIC FILE
            # dimensions:
            #     lat_110 = 224 ;
            #     lon_110 = 464 ;

            lsm_file_data["description"] = "NLDAS"
            lsm_file_data["weight_file_name"] = r'weight_nldas\.csv'
            lsm_file_data["grid_type"] = 'nldas'
        else:
            lsm_example_file.close()
            raise Exception("Unsupported runoff grid.")

    else:
        title = ""
        try:
            title = lsm_example_file.getncattr("TITLE")
        except AttributeError:
            pass

        if "WRF" in title:
            lsm_file_data["description"] = "WRF/WRF-Hydro Runoff"
            lsm_file_data["weight_file_name"] = r'weight_wrf\.csv'
            lsm_file_data["model_name"] = 'wrf'
            lsm_file_data["grid_type"] = 'wrf'

            lsm_file_data['rapid_inflow_tool'] = \
                CreateInflowFileFromWRFHydroRunoff(
                    latitude_dim,
                    longitude_dim,
                    latitude_var,
                    longitude_var,
                    surface_runoff_var,
                    subsurface_runoff_var,
                )
        else:
            lsm_example_file.close()
            raise Exception("Unsupported LSM grid.")

    lsm_example_file.close()

    # set the inflow tool to use the LDAS tool by default
    if lsm_file_data["rapid_inflow_tool"] is None:
        lsm_file_data["rapid_inflow_tool"] = \
            CreateInflowFileFromLDASRunoff(
                latitude_dim,
                longitude_dim,
                latitude_var,
                longitude_var,
                runoff_vars,
            )

    return lsm_file_data


def determine_start_end_timestep(lsm_file_list, file_re_match=None, file_datetime_pattern=None,
                                 expected_time_step=None, lsm_grid_info=None):
    """
    Determine the start and end date from LSM input files
    """

    if lsm_grid_info is None:
        lsm_grid_info = identify_lsm_grid(lsm_file_list[0])

    if None in (lsm_grid_info['time_var'], lsm_grid_info['time_dim'])\
            or lsm_grid_info['model_name'] in ('era_20cm', 'erai'):
        # NOTE: the ERA20CM and ERA 24hr time variables in the tests are erroneous
        if None in (file_re_match, file_datetime_pattern):
            raise ValueError("LSM files missing time dimension and/or variable."
                             "To mitigate this, add the 'file_re_match' and "
                             "'file_datetime_pattern' arguments.")

        if lsm_grid_info['time_dim'] is None:
            print("Assuming time dimension is 1")
            file_size_time = 1
        else:
            lsm_example_file = Dataset(lsm_file_list[0])
            file_size_time = len(lsm_example_file.dimensions[lsm_grid_info['time_dim']])
            lsm_example_file.close()

        total_num_time_steps = int(file_size_time * len(lsm_file_list))

        # determine the start time from the existing files
        actual_simulation_start_datetime = datetime.strptime(file_re_match.search(lsm_file_list[0]).group(0),
                                                             file_datetime_pattern)

        # check to see if the time step matches expected
        if len(lsm_file_list) > 1:
            time_step = int((datetime.strptime(file_re_match.search(lsm_file_list[1]).group(0), file_datetime_pattern)
                             - actual_simulation_start_datetime).total_seconds()
                            / float(file_size_time))

        elif expected_time_step is not None:
            time_step = int(expected_time_step)
        else:
            raise ValueError("Only one LSM file with one timestep present. "
                             "'expected_time_step' parameter required to continue.")

        # determine the end datetime
        actual_simulation_end_datetime = \
            datetime.strptime(file_re_match.search(lsm_file_list[-1]).group(0),
                              file_datetime_pattern) \
            + timedelta(seconds=(file_size_time-1) * time_step)
    else:
        with pangaea.open_mfdataset(lsm_file_list,
                                    lat_var=lsm_grid_info['latitude_var'],
                                    lon_var=lsm_grid_info['longitude_var'],
                                    time_var=lsm_grid_info['time_var'],
                                    lat_dim=lsm_grid_info['latitude_dim'],
                                    lon_dim=lsm_grid_info['longitude_dim'],
                                    time_dim=lsm_grid_info['time_dim']) as xds:

            datetime_arr = [pd.to_datetime(dval) for dval in xds.lsm.datetime.values]
            actual_simulation_start_datetime = datetime_arr[0]
            actual_simulation_end_datetime = datetime_arr[-1]
            total_num_time_steps = len(datetime_arr)

            if total_num_time_steps <= 1:
                if expected_time_step is not None:
                    time_step = int(expected_time_step)
                else:
                    raise ValueError("Only one LSM file with one timestep present. "
                                     "'expected_time_step' parameter required to continue.")

            else:
                time_step = int(np.diff(xds.lsm.datetime.values)[0] / np.timedelta64(1, 's'))

    if expected_time_step is not None:
        if time_step != int(expected_time_step):
            print("WARNING: The time step used {0} is different than expected {1}".format(time_step,
                                                                                          expected_time_step))

    return actual_simulation_start_datetime, actual_simulation_end_datetime, time_step, total_num_time_steps

# ------------------------------------------------------------------------------
# MAIN PROCESS
# ------------------------------------------------------------------------------
def run_lsm_rapid_process(rapid_executable_location,
                          lsm_data_location,
                          rapid_io_files_location=None,
                          rapid_input_location=None,
                          rapid_output_location=None,
                          simulation_start_datetime=None,
                          simulation_end_datetime=datetime.utcnow(),
                          file_datetime_pattern=None,
                          file_datetime_re_pattern=None,
                          initial_flows_file=None,
                          ensemble_list=[None],
                          generate_rapid_namelist_file=True,
                          run_rapid_simulation=True,
                          generate_return_periods_file=False,
                          return_period_method='weibul',
                          generate_seasonal_averages_file=False,
                          generate_seasonal_initialization_file=False,
                          generate_initialization_file=False,
                          use_all_processors=True,
                          num_processors=1,
                          mpiexec_command="mpiexec",
                          cygwin_bin_location="",
                          modeling_institution="US Army Engineer Research and Development Center",
                          convert_one_hour_to_three=False,
                          expected_time_step=None,
                          ):
    """
    This is the main process to generate inflow for RAPID and to run RAPID.

    Args:
        rapid_executable_location(str): Path to the RAPID executable.
        lsm_data_location(str): Path to the directory containing the Land Surface Model output files.
        rapid_io_files_location(Optional[str]): Path to the directory containing the input and output folders for RAPID. This is for running multiple watersheds.
        rapid_input_location(Optional[str]): Path to directory with RAPID simulation input data. Required if `rapid_io_files_location` is not set.
        rapid_output_location(Optional[str]): Path to directory to put output. Required if `rapid_io_files_location` is not set.
        simulation_start_datetime(Optional[datetime]): Datetime object with date bound of earliest simulation start.
        simulation_end_datetime(Optional[datetime]): Datetime object with date bound of latest simulation end. Defaults to datetime.utcnow().
        file_datetime_pattern(Optional[str]): Datetime pattern for files (Ex. '%Y%m%d%H'). If set, file_datetime_re_pattern is required. Various defaults used by each model.
        file_datetime_re_pattern(Optional[raw str]): Regex pattern to extract datetime (Ex. r'\d{10}'). If set, file_datetime_pattern is required. Various defaults used by each model.
        initial_flows_file(Optional[str]): If given, this is the path to a file with initial flows for the simulaion.
        ensemble_list(Optional[list]): This is the expexted ensemble name appended to the end of the file name.
        generate_rapid_namelist_file(Optional[bool]): If True, this will create a RAPID namelist file for the run in your RAPID input directory. Default is True.
        run_rapid_simulation(Optional[bool]): If True, the RAPID simulation will run after generating the inflow file. Default is True.
        generate_return_periods_file(Optional[bool]): If True, the return period file will be generated in the output. Default is False.
        return_period_method(Optional[str]): If True, the return period file will be generated in the output. Default is False.
        generate_seasonal_averages_file(Optional[bool]): If True, the season average file will be generated. Default is False.
        generate_seasonal_initialization_file(Optional[bool]): If True, an intialization based on the seasonal average for the current day of the year will be created. Default is False.
        generate_initialization_file(Optional[bool]): If True, an initialization file from the last time step of the simulation willl be created. Default is False.
        use_all_processors(Optional[bool]): If True, it will use all available processors to perform this operation. Default is True.
        num_processors(Optional[int]): If use_all_processors is False, this argument will determine the number of processors to use. Default is 1.
        mpiexec_command(Optional[str]): This is the command to execute RAPID. Default is "mpiexec".
        cygwin_bin_location(Optional[str]): If using Windows, this is the path to the Cygwin bin location. Default is "".
        modeling_institution(Optional[str]): This is the institution performing the modeling and is in the output files. Default is "US Army Engineer Research and Development Center".
        convert_one_hour_to_three(Optional[bool]): If the time step is expected to be 1-hr it will convert to 3. Set to False if the LIS, NLDAS, or Joules grid time step is greater than 1-hr.
        expected_time_step(Optional[int]): The time step in seconds of your LSM input data if only one file is given. Required if only one file is present.

    Returns:
        list: A list of output file information.

    Example of regular run:

    .. code:: python

        from datetime import datetime
        from RAPIDpy.inflow import run_lsm_rapid_process
        #------------------------------------------------------------------------------
        #main process
        #------------------------------------------------------------------------------
        if __name__ == "__main__":
            run_lsm_rapid_process(
                rapid_executable_location='/home/alan/rapid/src/rapid',
                rapid_io_files_location='/home/alan/rapid-io',
                lsm_data_location='/home/alan/era_data',
            )

    Example of single input/output run:

    .. code:: python

        from datetime import datetime
        from RAPIDpy.inflow import run_lsm_rapid_process
        #------------------------------------------------------------------------------
        #main process
        #------------------------------------------------------------------------------
        if __name__ == "__main__":
            run_lsm_rapid_process(
                rapid_executable_location='/home/alan/rapid/src/rapid',
                rapid_input_location='/home/alan/rapid-io/input/provo_watershed',
                rapid_output_location='/home/alan/rapid-io/output/provo_watershed',
                lsm_data_location='/home/alan/era_data',
            )
    Example of run with FLDAS and datetime filter:

    .. note:: http://disc.sci.gsfc.nasa.gov/uui/datasets?keywords=FLDAS

    .. code:: python

        from datetime import datetime
        from RAPIDpy.inflow import run_lsm_rapid_process
        #------------------------------------------------------------------------------
        #main process
        #------------------------------------------------------------------------------
        if __name__ == "__main__":
            run_lsm_rapid_process(
                rapid_executable_location='/home/alan/rapid/src/rapid',
                rapid_io_files_location='/home/alan/rapid-io',
                lsm_data_location='/home/alan/lsm_data',
                simulation_start_datetime=datetime(1980, 1, 1),
                file_datetime_re_pattern = r'\d{8}',
                file_datetime_pattern = "%Y%m%d",
            )

    Example of run with CMIP5:

    .. note:: http://gdo-dcp.ucllnl.org/downscaled_cmip_projections/techmemo/BCSD5HydrologyMemo.pdf

    .. code:: python

        from datetime import datetime
        from RAPIDpy.inflow import run_lsm_rapid_process
        #------------------------------------------------------------------------------
        #main process
        #------------------------------------------------------------------------------
        if __name__ == "__main__":
            run_lsm_rapid_process(
                rapid_executable_location='/home/jimwlewis/rapid/src/rapid',
                rapid_io_files_location='/data/rapid-io4',
                lsm_data_location='/data/rapid-io4/input/cmip5-jun01',
                simulation_start_datetime=datetime(2001, 1, 1),
                simulation_end_datetime=datetime(2002, 12, 31),
                file_datetime_pattern="%Y",
                file_datetime_re_pattern=r'\d{4}',
            )
    """
    time_begin_all = datetime.utcnow()

    # use all processors makes precedent over num_processors arg
    if use_all_processors == True:
        NUM_CPUS = multiprocessing.cpu_count()
    elif num_processors > multiprocessing.cpu_count():
        print("WARNING: Num processors requested exceeded max. Set to max ...")
        NUM_CPUS = multiprocessing.cpu_count()
    else:
        NUM_CPUS = num_processors

    # get list of correclty formatted rapid input directories in rapid directory

    rapid_directories = []
    if rapid_io_files_location is not None:
        main_rapid_input_directory = os.path.join(rapid_io_files_location, 'input')
        for watershed_directory in get_valid_directory_list(main_rapid_input_directory):
            watershed_input_path = os.path.join(main_rapid_input_directory,
                                                watershed_directory)
            watershed_output_path = os.path.join(rapid_io_files_location,
                                                 'output',
                                                 watershed_directory)
            rapid_directories.append((watershed_input_path, watershed_output_path))
    elif None not in (rapid_input_location, rapid_output_location):
        rapid_directories = [(rapid_input_location, rapid_output_location)]
    else:
        raise ValueError("Need 'rapid_io_files_location' or 'rapid_input_location' "
                         "and 'rapid_output_location' set to continue.")

    all_output_file_information = []

    for ensemble in ensemble_list:
        output_file_information = {
            'ensemble' : ensemble,
        }
        ensemble_file_ending = ".nc"
        ensemble_file_ending4 = ".nc4"
        if ensemble != None:
            ensemble_file_ending = "_{0}.nc".format(ensemble)
            ensemble_file_ending4 = "_{0}.nc4".format(ensemble)

        # get list of files
        lsm_file_list = []
        for subdir, dirs, files in os.walk(lsm_data_location, followlinks=True):
            for lsm_file in files:
                if lsm_file.endswith(ensemble_file_ending) or lsm_file.endswith(ensemble_file_ending4):
                    lsm_file_list.append(os.path.join(subdir, lsm_file))
        lsm_file_list = sorted(lsm_file_list)

        # IDENTIFY THE GRID
        lsm_file_data = identify_lsm_grid(lsm_file_list[0])

        # load in the datetime pattern
        if file_datetime_pattern is None or file_datetime_re_pattern is None:
            file_datetime_re_pattern = DEFAULT_LSM_INPUTS[lsm_file_data['grid_type']]['file_datetime_re_pattern']
            file_datetime_pattern = DEFAULT_LSM_INPUTS[lsm_file_data['grid_type']]['file_datetime_pattern']
        file_re_match = re.compile(file_datetime_re_pattern)

        # get subset based on time bounds
        if simulation_start_datetime is not None:
            print("Filtering files by datetime ...")
            lsm_file_list_subset = []
            for lsm_file in lsm_file_list:
                match = file_re_match.search(lsm_file)
                file_date = datetime.strptime(match.group(0), file_datetime_pattern)
                if file_date > simulation_end_datetime:
                    break
                if file_date >= simulation_start_datetime:
                    lsm_file_list_subset.append(os.path.join(subdir, lsm_file))

            lsm_file_list = sorted(lsm_file_list_subset)

        print("Running from {0} to {1}".format(lsm_file_list[0],
                                               lsm_file_list[-1]))

        # get number of time steps in file
        actual_simulation_start_datetime, actual_simulation_end_datetime, time_step, total_num_time_steps = \
            determine_start_end_timestep(lsm_file_list,
                                         file_re_match=file_re_match,
                                         file_datetime_pattern=file_datetime_pattern,
                                         expected_time_step=expected_time_step,
                                         lsm_grid_info=lsm_file_data)

        # VALIDATING INPUT IF DIVIDING BY 3
        time_step_multiply_factor = 1
        if (lsm_file_data['grid_type'] in ('nldas', 'lis', 'joules')) and convert_one_hour_to_three:
            num_extra_files = total_num_time_steps % 3
            if num_extra_files != 0:
                print("WARNING: Number of files needs to be divisible by 3. Remainder is {0}".format(num_extra_files))
                print("This means your simulation will be truncated")
            total_num_time_steps /= 3
            time_step *= 3

        # compile the file ending
        out_file_ending = "{0}_{1}_{2}hr_{3:%Y%m%d}to{4:%Y%m%d}{5}".format(lsm_file_data['model_name'],
                                                             lsm_file_data['grid_type'],
                                                             int(time_step/3600),
                                                             actual_simulation_start_datetime,
                                                             actual_simulation_end_datetime,
                                                             ensemble_file_ending)

        # run LSM processes
        for master_watershed_input_directory, master_watershed_output_directory in rapid_directories:
            print("Running from: {0}".format(master_watershed_input_directory))
            try:
                os.makedirs(master_watershed_output_directory)
            except OSError:
                pass

            # create inflow to dump data into
            master_rapid_runoff_file = os.path.join(master_watershed_output_directory,
                                                    'm3_riv_bas_{0}'.format(out_file_ending))

            weight_table_file = case_insensitive_file_search(master_watershed_input_directory,
                                                             lsm_file_data['weight_file_name'])

            try:
                in_rivid_lat_lon_z_file = case_insensitive_file_search(master_watershed_input_directory,
                                                                       r'comid_lat_lon_z\.csv')
            except Exception:
                in_rivid_lat_lon_z_file = ""
                print("WARNING: comid_lat_lon_z file not found. The lat/lon will not be added ...")
                pass

            print("Writing inflow file to: {0}".format(master_rapid_runoff_file))
            lsm_file_data['rapid_inflow_tool'].generateOutputInflowFile(
                out_nc=master_rapid_runoff_file,
                start_datetime_utc=actual_simulation_start_datetime,
                number_of_timesteps=total_num_time_steps,
                simulation_time_step_seconds=time_step,
                in_rapid_connect_file=case_insensitive_file_search(master_watershed_input_directory,
                                                                   r'rapid_connect\.csv'),
                in_rivid_lat_lon_z_file=in_rivid_lat_lon_z_file,
                land_surface_model_description=lsm_file_data['description'],
                modeling_institution=modeling_institution
            )

            job_combinations = []
            if (lsm_file_data['grid_type'] in ('nldas', 'lis', 'joules')) and convert_one_hour_to_three:
                print("Grouping {0} in threes".format(lsm_file_data['grid_type']))
                lsm_file_list = [lsm_file_list[nldas_index:nldas_index+3]
                                 for nldas_index in range(0, len(lsm_file_list), 3)
                                 if len(lsm_file_list[nldas_index:nldas_index+3]) == 3]

            if len(lsm_file_list) < NUM_CPUS:
                NUM_CPUS = len(lsm_file_list)
            mp_lock = multiprocessing.Manager().Lock()
            partition_list, partition_index_list = partition(lsm_file_list, NUM_CPUS)

            for loop_index, cpu_grouped_file_list in enumerate(partition_list):
                if cpu_grouped_file_list and partition_index_list[loop_index]:
                    job_combinations.append((cpu_grouped_file_list,
                                             partition_index_list[loop_index],
                                             weight_table_file,
                                             lsm_file_data['grid_type'],
                                             master_rapid_runoff_file,
                                             lsm_file_data['rapid_inflow_tool'],
                                             mp_lock))
                    # COMMENTED CODE IS FOR DEBUGGING
#                    generate_inflows_from_runoff((cpu_grouped_file_list,
#                                                  partition_index_list[loop_index],
#                                                  lsm_file_data['weight_table_file'],
#                                                  lsm_file_data['grid_type'],
#                                                  master_rapid_runoff_file,
#                                                  lsm_file_data['rapid_inflow_tool'],
#                                                  mp_lock))
            pool = multiprocessing.Pool(NUM_CPUS)
            pool.map(generate_inflows_from_runoff,
                     job_combinations)
            pool.close()
            pool.join()

            # set up RAPID manager
            rapid_manager = RAPID(rapid_executable_location=rapid_executable_location,
                                  cygwin_bin_location=cygwin_bin_location,
                                  num_processors=NUM_CPUS,
                                  mpiexec_command=mpiexec_command,
                                  ZS_TauR=time_step,  # duration of routing procedure (time step of runoff data)
                                  ZS_dtR=15 * 60,  # internal routing time step
                                  ZS_TauM=total_num_time_steps * time_step,  # total simulation time
                                  ZS_dtM=time_step  # RAPID recommended internal time step (1 day)
                                  )

            if initial_flows_file and os.path.exists(initial_flows_file):
                rapid_manager.update_parameters(
                    Qinit_file=initial_flows_file,
                    BS_opt_Qinit=True
                )

            # run RAPID for the watershed
            lsm_rapid_output_file = os.path.join(master_watershed_output_directory,
                                                 'Qout_{0}'.format(out_file_ending))
            rapid_manager.update_parameters(
                rapid_connect_file=case_insensitive_file_search(master_watershed_input_directory,
                                                                r'rapid_connect\.csv'),
                Vlat_file=master_rapid_runoff_file,
                riv_bas_id_file=case_insensitive_file_search(master_watershed_input_directory,
                                                             r'riv_bas_id\.csv'),
                k_file=case_insensitive_file_search(master_watershed_input_directory,
                                                    r'k\.csv'),
                x_file=case_insensitive_file_search(master_watershed_input_directory,
                                                    r'x\.csv'),
                Qout_file=lsm_rapid_output_file
            )

            rapid_manager.update_reach_number_data()


            output_file_information[os.path.basename(master_watershed_input_directory)] = {
                'm3_riv': master_rapid_runoff_file,
                'qout': lsm_rapid_output_file
            }


            if generate_rapid_namelist_file:
                rapid_manager.generate_namelist_file(os.path.join(master_watershed_input_directory,
                                                                  "rapid_namelist_{}".format(out_file_ending[:-3])))
            if run_rapid_simulation:
                rapid_manager.run()

                rapid_manager.make_output_CF_compliant(
                    simulation_start_datetime=actual_simulation_start_datetime,
                    comid_lat_lon_z_file=in_rivid_lat_lon_z_file,
                    project_name="{0} Based Historical flows by {1}".format(lsm_file_data['description'],
                                                                            modeling_institution)
                )

                # generate return periods
                if generate_return_periods_file and os.path.exists(lsm_rapid_output_file) and lsm_rapid_output_file:
                    return_periods_file = os.path.join(master_watershed_output_directory,
                                                       'return_periods_{0}'.format(out_file_ending))
                    # assume storm has 3 day length
                    storm_length_days = 3
                    generate_return_periods(qout_file=lsm_rapid_output_file,
                                            return_period_file=return_periods_file,
                                            num_cpus=NUM_CPUS,
                                            storm_duration_days=storm_length_days,
                                            method=return_period_method)
                                            
                # generate seasonal averages file
                if generate_seasonal_averages_file and os.path.exists(lsm_rapid_output_file) and lsm_rapid_output_file:
                    seasonal_averages_file = os.path.join(master_watershed_output_directory,
                                                          'seasonal_averages_{0}'.format(out_file_ending))
                    generate_seasonal_averages(lsm_rapid_output_file,
                                               seasonal_averages_file,
                                               NUM_CPUS)

                # generate seasonal initialization file
                if generate_seasonal_initialization_file and os.path.exists(lsm_rapid_output_file) and lsm_rapid_output_file:
                    seasonal_qinit_file = os.path.join(master_watershed_input_directory,
                                                       'seasonal_qinit_{0}.csv'.format(out_file_ending[:-3]))
                    rapid_manager.generate_seasonal_intitialization(seasonal_qinit_file)

                # generate initialization file
                if generate_initialization_file and os.path.exists(lsm_rapid_output_file) and lsm_rapid_output_file:
                    qinit_file = os.path.join(master_watershed_input_directory,
                                              'qinit_{0}.csv'.format(out_file_ending[:-3]))
                    rapid_manager.generate_qinit_from_past_qout(qinit_file)

        all_output_file_information.append(output_file_information)

    # print info to user
    time_end = datetime.utcnow()
    print("Time Begin All: {0}".format(time_begin_all))
    print("Time Finish All: {0}".format(time_end))
    print("TOTAL TIME: {0}".format(time_end-time_begin_all))

    return all_output_file_information