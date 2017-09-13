# -*- coding: utf-8 -*-
##
##  test_rapidpy.py
##  RAPIDpy
##
##  Created by Alan D. Snow 2016.
##  Copyright © 2015 Alan D Snow. All rights reserved.
##

from datetime import datetime
from filecmp import cmp as fcmp
from netCDF4 import Dataset
import numpy as np
from numpy.testing import assert_array_equal, assert_almost_equal
import os
import pytest
from pytz import timezone
from shutil import copy

#local import
from RAPIDpy import RAPID
from RAPIDpy import RAPIDDataset
from RAPIDpy.dataset import compare_qout_files
from RAPIDpy.helper_functions import (compare_csv_decimal_files,
                                      compare_csv_timeseries_files,
                                      remove_files)

from RAPIDpy.postprocess import find_goodness_of_fit
from RAPIDpy.postprocess import ConvertRAPIDOutputToCF

#GLOBAL VARIABLES
MAIN_TESTS_FOLDER = os.path.dirname(os.path.abspath(__file__))
COMPARE_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'compare')
INPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'data')
OUTPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'output')
RAPID_EXE_PATH = os.path.join(MAIN_TESTS_FOLDER,
                              "..", "..",
                              "rapid", "src", "rapid")
CYGWIN_BIN_PATH = 'C:\\cygwin64\\bin'

#------------------------------------------------------------------------------
# MAIN TEST SCRIPTS
#------------------------------------------------------------------------------
def test_generate_rapid_input_file():
    """
    Checks RAPID input file generation with valid input
    """
    print("TEST 1: GENERATE NAMELIST FILE")
    rapid_manager = RAPID(rapid_executable_location=RAPID_EXE_PATH,
                          cygwin_bin_location=CYGWIN_BIN_PATH,
                          use_all_processors=True,
                          ZS_TauR = 24*3600, #duration of routing procedure (time step of runoff data)
                          ZS_dtR = 15*60, #internal routing time step
                          ZS_TauM = 12*24*3600, #total simulation time
                          ZS_dtM = 24*3600 #input time step
                         )
    rapid_manager.update_parameters(rapid_connect_file='rapid_connect.csv',
                                    Vlat_file='m3_riv.nc',
                                    riv_bas_id_file='riv_bas_id.csv',
                                    k_file='k.csv',
                                    x_file='x.csv',
                                    Qout_file='Qout.nc'
                                    )
    generated_input_file = os.path.join(OUTPUT_DATA_PATH,
                                        "rapid_namelist-GENERATE")
    rapid_manager.generate_namelist_file(generated_input_file)
    generated_input_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                 "rapid_namelist-GENERATE")


    assert (fcmp(generated_input_file, generated_input_file_solution))

    remove_files(generated_input_file)

def test_update_rapid_input_file():
    """
    Checks RAPID input file update with valid input
    """
    print("TEST 2: UPDATE NAMELIST FILE")
    rapid_manager = RAPID(rapid_executable_location=RAPID_EXE_PATH,
                          cygwin_bin_location=CYGWIN_BIN_PATH,
                          use_all_processors=True,
                         )
    rapid_manager.update_parameters(rapid_connect_file='rapid_connect.csv',
                                    Vlat_file='m3_riv.nc',
                                    riv_bas_id_file='riv_bas_id.csv',
                                    k_file='k.csv',
                                    x_file='x.csv',
                                    Qout_file='Qout.nc'
                                    )

    original_input_file = os.path.join(INPUT_DATA_PATH,
                                      "rapid_namelist_valid")

    updated_input_file = os.path.join(OUTPUT_DATA_PATH,
                                      "rapid_namelist-UPDATE")

    copy(original_input_file, updated_input_file)
    rapid_manager.update_namelist_file(updated_input_file)
    updated_input_file_solution = os.path.join(COMPARE_DATA_PATH,
                                               "rapid_namelist-UPDATE")


    assert (fcmp(updated_input_file, updated_input_file_solution))

    remove_files(updated_input_file)

def test_update_rapid_invalid_input_file():
    """
    Checks RAPID input file update with valid input
    """
    print("TEST 3: UPDATE WITH INVALID NAMELIST FILE")
    rapid_manager = RAPID(rapid_executable_location=RAPID_EXE_PATH,
                          cygwin_bin_location=CYGWIN_BIN_PATH,
                          use_all_processors=True,
                         )
    rapid_manager.update_parameters(rapid_connect_file='rapid_connect.csv',
                                    Vlat_file='m3_riv.nc',
                                    riv_bas_id_file='riv_bas_id.csv',
                                    k_file='k.csv',
                                    x_file='x.csv',
                                    Qout_file='Qout.nc'
                                    )

    original_input_file = os.path.join(INPUT_DATA_PATH,
                                      "rapid_namelist_invalid")

    updated_input_file = os.path.join(OUTPUT_DATA_PATH,
                                      "rapid_namelist-UPDATE-INVALID")

    copy(original_input_file, updated_input_file)
    rapid_manager.update_namelist_file(updated_input_file)
    updated_input_file_solution = os.path.join(COMPARE_DATA_PATH,
                                               "rapid_namelist-UPDATE-INVALID")


    assert (fcmp(updated_input_file, updated_input_file_solution))

    remove_files(updated_input_file)

def test_update_rapid_numbers_input_file():
    """
    Checks RAPID input file update with number validation
    """
    print("TEST 4: GENERATE NUMBERS FOR NAMELIST FILE")
    rapid_manager = RAPID(rapid_executable_location=RAPID_EXE_PATH,
                          cygwin_bin_location=CYGWIN_BIN_PATH,
                          use_all_processors=True,
                          rapid_connect_file=os.path.join(INPUT_DATA_PATH, 'rapid_connect.csv'),
                          riv_bas_id_file=os.path.join(INPUT_DATA_PATH, 'riv_bas_id.csv'),
                         )
    rapid_manager.update_reach_number_data()

    rapid_manager.update_parameters(rapid_connect_file='rapid_connect.csv',
                                    Vlat_file='m3_nasa_lis_3hr_20020830.nc',
                                    riv_bas_id_file='riv_bas_id.csv',
                                    k_file='k.csv',
                                    x_file='x.csv',
                                    Qout_file='Qout.nc'
                                    )

    generated_input_file = os.path.join(OUTPUT_DATA_PATH,
                                      "rapid_namelist-GENERATE-NUMBERS")

    rapid_manager.generate_namelist_file(generated_input_file)

    generated_input_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                 "rapid_namelist-GENERATE-NUMBERS")


    assert (fcmp(generated_input_file, generated_input_file_solution))

    remove_files(generated_input_file)

def test_update_rapid_numbers_forcing_input_file():
    """
    Checks RAPID input file update with forcing data and number validation
    """
    rapid_manager = RAPID(rapid_executable_location=RAPID_EXE_PATH,
                          cygwin_bin_location=CYGWIN_BIN_PATH,
                          use_all_processors=True,
                          rapid_connect_file=os.path.join(INPUT_DATA_PATH,
                                                          'rapid_connect.csv'),
                          riv_bas_id_file=os.path.join(INPUT_DATA_PATH,
                                                       'riv_bas_id.csv'),
                          for_tot_id_file=os.path.join(INPUT_DATA_PATH,
                                                       'for_tot_id.csv'),
                          for_use_id_file=os.path.join(INPUT_DATA_PATH,
                                                       'for_use_id.csv'),
                          ZS_dtF=3 * 60 * 60,
                          BS_opt_for=True
                          )

    rapid_manager.update_reach_number_data()

    rapid_manager.update_parameters(rapid_connect_file='rapid_connect.csv',
                                    Vlat_file='m3_nasa_lis_3hr_20020830.nc',
                                    riv_bas_id_file='riv_bas_id.csv',
                                    k_file='k.csv',
                                    x_file='x.csv',
                                    Qout_file='Qout.nc',
                                    Qfor_file='qfor.csv',
                                    for_tot_id_file='for_tot_id.csv',
                                    for_use_id_file='for_use_id.csv',
                                    )

    generated_input_file = os.path.join(OUTPUT_DATA_PATH,
                                      "rapid_namelist-GENERATE-NUMBERS-FORCING")

    rapid_manager.generate_namelist_file(generated_input_file)

    generated_input_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                 "rapid_namelist-GENERATE-NUMBERS-FORCING")


    assert (fcmp(generated_input_file, generated_input_file_solution))

    remove_files(generated_input_file)

def test_update_rapid_runtime():
    """
    Checks RAPID input file update with number validation
    """
    print("TEST 5: GENERATE RUNTIME FOR NAMELIST FILE")
    rapid_manager = RAPID(rapid_executable_location=RAPID_EXE_PATH,
                          cygwin_bin_location=CYGWIN_BIN_PATH,
                          use_all_processors=True,
                          rapid_connect_file=os.path.join(INPUT_DATA_PATH, 'rapid_connect.csv'),
                          riv_bas_id_file=os.path.join(INPUT_DATA_PATH, 'riv_bas_id.csv'),
                          Vlat_file=os.path.join(INPUT_DATA_PATH, 'm3_nasa_lis_3hr_20020830.nc'),
                          ZS_TauR=3*3600,
                         )

    rapid_manager.update_reach_number_data()

    rapid_manager.update_simulation_runtime()

    rapid_manager.update_parameters(rapid_connect_file='rapid_connect.csv',
                                    Vlat_file='m3_nasa_lis_3hr_20020830.nc',
                                    riv_bas_id_file='riv_bas_id.csv',
                                    k_file='k.csv',
                                    x_file='x.csv',
                                    Qout_file='Qout.nc'
                                    )

    generated_input_file = os.path.join(OUTPUT_DATA_PATH,
                                        "rapid_namelist-GENERATE-TIME")

    rapid_manager.generate_namelist_file(generated_input_file)

    generated_input_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                 "rapid_namelist-GENERATE-TIME")


    assert (fcmp(generated_input_file, generated_input_file_solution))

    remove_files(generated_input_file)

def test_qout_same():
    """
    Test function to compare RAPID simulation Qout
    """
    print("TEST 6: TEST COMPARE RAPID QOUT")
    input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830.nc')
    input_qout_file_cf = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_CF.nc')

    assert (compare_qout_files(input_qout_file, input_qout_file_cf))

@pytest.mark.skipif(not os.path.exists(RAPID_EXE_PATH), reason='Only run if RAPID installed')
def test_run_rapid_simulation():
    """
    Test Running RAPID Simulation
    """

    print("TEST 7: TEST RUNNING RAPID SIMULATION")
    generated_qout_file = os.path.join(OUTPUT_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_generated.nc')



    rapid_manager = RAPID(rapid_executable_location=RAPID_EXE_PATH,
                          cygwin_bin_location=CYGWIN_BIN_PATH,
                          num_processors=1,
                          rapid_connect_file=os.path.join(INPUT_DATA_PATH, 'rapid_connect.csv'),
                          riv_bas_id_file=os.path.join(INPUT_DATA_PATH, 'riv_bas_id.csv'),
                          Vlat_file=os.path.join(INPUT_DATA_PATH, 'm3_nasa_lis_3hr_20020830.nc'),
                          k_file=os.path.join(INPUT_DATA_PATH, 'k.csv'),
                          x_file=os.path.join(INPUT_DATA_PATH, 'x.csv'),
                          ZS_dtM=10800,
                          ZS_dtR=900,
                          ZS_TauM=2*86400,
                          ZS_TauR=10800,
                          Qout_file=generated_qout_file
                         )
    rapid_manager.update_reach_number_data()
    rapid_manager.run()

    generated_qout_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                'Qout_nasa_lis_3hr_20020830.nc')

    #check Qout
    assert (compare_qout_files(generated_qout_file, generated_qout_file_solution))

    #check other info in netcdf file
    d1 = Dataset(generated_qout_file)
    d2 = Dataset(generated_qout_file_solution)
    assert (d1.dimensions.keys() == d2.dimensions.keys())
    assert (d1.variables.keys() == d2.variables.keys())
    assert ((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    d1.close()
    d2.close()

    remove_files(generated_qout_file)

def test_convert_file_to_be_cf_compliant_new_format_comid_lat_lon_z():
    """
    Test Convert RAPID Output to be CF Compliant for new format with COMID_LAT_LON_Z
    """
    print("TEST 8: TEST CONVERT RAPID OUTPUT TO CF COMPLIANT (COMID_LAT_LON_Z)")

    input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830.nc')
    temp_qout_file = os.path.join(OUTPUT_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_test_cf_lat_lon_z.nc')
    copy(input_qout_file, temp_qout_file)

    rapid_manager = RAPID(rapid_executable_location=RAPID_EXE_PATH,
                          cygwin_bin_location=CYGWIN_BIN_PATH,
                          Qout_file=temp_qout_file,
                          rapid_connect_file=os.path.join(INPUT_DATA_PATH, 'rapid_connect.csv'),
                          ZS_TauR=3*3600)

    rapid_manager.make_output_CF_compliant(simulation_start_datetime=datetime(2002, 8, 30),
                                           comid_lat_lon_z_file=os.path.join(INPUT_DATA_PATH, 'comid_lat_lon_z.csv'),
                                           project_name="ERA Interim (T511 Grid) 3 Hourly Runoff Based Historical flows by US Army ERDC")

    cf_qout_file_solution = os.path.join(COMPARE_DATA_PATH,
                                         'Qout_nasa_lis_3hr_20020830_CF.nc')

    #check Qout
    assert (compare_qout_files(temp_qout_file, cf_qout_file_solution))

    #check other info in netcdf file
    d1 = Dataset(temp_qout_file)
    d2 = Dataset(cf_qout_file_solution)
    assert (d1.dimensions.keys() == d2.dimensions.keys())
    assert (d1.variables.keys() == d2.variables.keys())
    assert ((d1.variables['time'][:] == d2.variables['time'][:]).all())
    assert ((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    assert ((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
    assert ((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
    d1.close()
    d2.close()

    remove_files(temp_qout_file)

def test_convert_file_to_be_cf_compliant_new_format():
    """
    Test Convert RAPID Output to be CF Compliant for new format without COMID_LAT_LON_Z
    """
    print("TEST 9: TEST CONVERT RAPID OUTPUT TO CF COMPLIANT (NO COMID_LAT_LON_Z)")

    input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830.nc')
    temp_qout_file = os.path.join(OUTPUT_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_test_cf.nc')
    copy(input_qout_file, temp_qout_file)

    rapid_manager = RAPID(rapid_executable_location=RAPID_EXE_PATH,
                          cygwin_bin_location=CYGWIN_BIN_PATH,
                          Qout_file=temp_qout_file,
                          rapid_connect_file=os.path.join(INPUT_DATA_PATH, 'rapid_connect.csv'),
                          ZS_TauR=3*3600)

    rapid_manager.make_output_CF_compliant(simulation_start_datetime=datetime(2002, 8, 30),
                                           comid_lat_lon_z_file="",
                                           project_name="ERA Interim (T511 Grid) 3 Hourly Runoff Based Historical flows by US Army ERDC")

    cf_qout_file_solution = os.path.join(COMPARE_DATA_PATH,
                                         'Qout_nasa_lis_3hr_20020830_CF_no_lat_lon_z.nc')

    #check Qout
    assert (compare_qout_files(temp_qout_file, cf_qout_file_solution))

    #check other info in netcdf file
    d1 = Dataset(temp_qout_file)
    d2 = Dataset(cf_qout_file_solution)
    assert (d1.dimensions.keys() == d2.dimensions.keys())
    assert (d1.variables.keys() == d2.variables.keys())
    assert ((d1.variables['time'][:] == d2.variables['time'][:]).all())
    assert ((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    d1.close()
    d2.close()

    remove_files(temp_qout_file)

def test_convert_file_to_be_cf_compliant_original_format():
    """
    Test Convert RAPID Output to be CF Compliant for original format
    """
    print("TEST 10: TEST CONVERT RAPID OUTPUT TO CF COMPLIANT - ORIGINAL (COMID_LAT_LON_Z)")

    input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_original.nc')
    temp_qout_file = os.path.join(OUTPUT_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_original_test_cf.nc')
    copy(input_qout_file, temp_qout_file)

    rapid_manager = RAPID(rapid_executable_location=RAPID_EXE_PATH,
                          cygwin_bin_location=CYGWIN_BIN_PATH,
                          Qout_file=temp_qout_file,
                          rapid_connect_file=os.path.join(INPUT_DATA_PATH, 'rapid_connect.csv'),
                          ZS_TauR=3*3600)

    rapid_manager.make_output_CF_compliant(simulation_start_datetime=datetime(2002, 8, 30),
                                           comid_lat_lon_z_file=os.path.join(INPUT_DATA_PATH, 'comid_lat_lon_z.csv'),
                                           project_name="ERA Interim (T511 Grid) 3 Hourly Runoff Based Historical flows by US Army ERDC")

    cf_qout_file_solution = os.path.join(COMPARE_DATA_PATH,
                                         'Qout_nasa_lis_3hr_20020830_CF.nc')

    #check Qout
    assert (compare_qout_files(temp_qout_file, cf_qout_file_solution))

    #check other info in netcdf file
    d1 = Dataset(temp_qout_file)
    d2 = Dataset(cf_qout_file_solution)
    assert (d1.dimensions.keys() == d2.dimensions.keys())
    assert (d1.variables.keys() == d2.variables.keys())
    assert ((d1.variables['time'][:] == d2.variables['time'][:]).all())
    assert ((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    assert ((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
    assert ((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
    d1.close()
    d2.close()

    remove_files(temp_qout_file)

def test_generate_qinit_file():
    """
    This tests the qinit file function to create an input qinit file for RAPID
    """
    print("TEST 11: TEST GENERATE QINIT FILE")
    rapid_manager = RAPID(rapid_executable_location=RAPID_EXE_PATH,
                          cygwin_bin_location=CYGWIN_BIN_PATH,
                          rapid_connect_file=os.path.join(INPUT_DATA_PATH, 'rapid_connect.csv')
                         )

    #test with original rapid outpur
    input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830.nc')
    original_qout_file = os.path.join(OUTPUT_DATA_PATH, 'Qout_nasa_lis_3hr_20020830.nc')
    copy(input_qout_file, original_qout_file)

    qinit_original_rapid_qout = os.path.join(OUTPUT_DATA_PATH, 'qinit_original_rapid_qout.nc')
    rapid_manager.update_parameters(Qout_file=original_qout_file)
    rapid_manager.generate_qinit_from_past_qout(qinit_file=qinit_original_rapid_qout)

    qinit_original_rapid_qout_solution = os.path.join(COMPARE_DATA_PATH, 'qinit_original_rapid_qout.nc')
    d1 = Dataset(qinit_original_rapid_qout)
    d2 = Dataset(qinit_original_rapid_qout_solution)
    assert (d1.dimensions.keys() == d2.dimensions.keys())
    assert (d1.variables.keys() == d2.variables.keys())
    assert (d1.variables['time'][:] is np.ma.masked)
    assert (d2.variables['time'][:] is np.ma.masked)
    assert_array_equal(d1.variables['rivid'][:], d2.variables['rivid'][:])
    assert_almost_equal(d1.variables['Qout'][:], d2.variables['Qout'][:])
    d1.close()
    d2.close()

    #test with CF rapid output and alternate time index
    cf_input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_CF.nc')
    cf_qout_file = os.path.join(OUTPUT_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_CF.nc')
    copy(cf_input_qout_file, cf_qout_file)

    qinit_cf_rapid_qout = os.path.join(OUTPUT_DATA_PATH, 'qinit_cf_rapid_qout.nc')
    rapid_manager.update_parameters(Qout_file=cf_qout_file)
    rapid_manager.generate_qinit_from_past_qout(qinit_file=qinit_cf_rapid_qout,
                                                time_index=5)
    qinit_cf_rapid_qout_solution = os.path.join(COMPARE_DATA_PATH, 'qinit_cf_rapid_qout.nc')

    d1 = Dataset(qinit_cf_rapid_qout)
    d2 = Dataset(qinit_cf_rapid_qout_solution)
    assert (d1.dimensions.keys() == d2.dimensions.keys())
    assert (d1.variables.keys() == d2.variables.keys())
    assert_array_equal(d1.variables['time'][:], d2.variables['time'][:])
    assert_array_equal(d1.variables['rivid'][:], d2.variables['rivid'][:])
    assert_almost_equal(d1.variables['Qout'][:], d2.variables['Qout'][:])
    d1.close()
    d2.close()

    remove_files(original_qout_file,
                 qinit_original_rapid_qout,
                 cf_qout_file,
                 qinit_cf_rapid_qout
                 )


def test_download_usgs_daily_avg():
    """
    This tests downloading USGS daily avg data
    """
    print("TEST 12: TEST DOWNLOAD USGS DAILY AVERAGE DATA")

    out_streamflow_file=os.path.join(OUTPUT_DATA_PATH,"gage_streamflow.csv")
    out_stream_id_file=os.path.join(OUTPUT_DATA_PATH,"gage_rivid.csv")

    rapid_manager = RAPID(rapid_executable_location=RAPID_EXE_PATH,
                          cygwin_bin_location=CYGWIN_BIN_PATH)

    rapid_manager.generate_usgs_avg_daily_flows_opt(reach_id_gage_id_file=os.path.join(INPUT_DATA_PATH,"usgs_gage_id_rivid.csv"),
                                                    start_datetime=datetime(2000,1,1),
                                                    end_datetime=datetime(2000,1,3),
                                                    out_streamflow_file=out_streamflow_file,
                                                    out_stream_id_file=out_stream_id_file)

    compare_streamflow_file=os.path.join(COMPARE_DATA_PATH,"gage_streamflow.csv")
    assert (compare_csv_decimal_files(out_streamflow_file, compare_streamflow_file, header=False))

    compare_stream_id_file=os.path.join(COMPARE_DATA_PATH,"gage_rivid.csv")
    assert (compare_csv_decimal_files(out_stream_id_file, compare_stream_id_file, header=False))

    remove_files(out_streamflow_file,
                 out_stream_id_file)


def test_extract_timeseries():
    """
    This tests extracting a timeseries from RAPID Qout file
    """
    print("TEST 13: TEST EXTRACT TIMESERIES FROM QINIT FILE")

    #for writing entire time series to file from new rapid output
    input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830.nc')
    new_qout_file = os.path.join(OUTPUT_DATA_PATH, 'Qout_nasa_lis_3hr_20020830.nc')
    copy(input_qout_file, new_qout_file)
    new_timeseries_file = os.path.join(OUTPUT_DATA_PATH, 'new_timeseries_file.csv')

    with RAPIDDataset(new_qout_file) as qout_nc:
        qout_nc.write_flows_to_csv(new_timeseries_file,
                                   river_id=75224)

        if qout_nc.is_time_variable_valid():
            original_timeseries_file_solution = os.path.join(COMPARE_DATA_PATH, 'original_timeseries.csv')
        else:
            original_timeseries_file_solution = os.path.join(COMPARE_DATA_PATH, 'original_timeseries-notime.csv')

    assert (compare_csv_timeseries_files(new_timeseries_file, original_timeseries_file_solution, header=False))

    #for writing entire time series to file from original rapid output
    input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_original.nc')
    original_qout_file = os.path.join(OUTPUT_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_original.nc')
    copy(input_qout_file, original_qout_file)
    original_timeseries_file = os.path.join(OUTPUT_DATA_PATH, 'original_timeseries.csv')

    with RAPIDDataset(original_qout_file) as qout_nc:
        qout_nc.write_flows_to_csv(original_timeseries_file,
                                   river_id=75224)
    original_timeseries_file_solution = os.path.join(COMPARE_DATA_PATH, 'original_timeseries-notime.csv')

    assert (compare_csv_timeseries_files(original_timeseries_file, original_timeseries_file_solution, header=False))

    #if file is CF compliant, you can write out daily average
    cf_input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_CF.nc')
    cf_qout_file = os.path.join(OUTPUT_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_CF.nc')
    copy(cf_input_qout_file, cf_qout_file)
    cf_timeseries_daily_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries_daily.csv')

    with RAPIDDataset(cf_qout_file) as qout_nc:
        qout_nc.write_flows_to_csv(cf_timeseries_daily_file,
                                   river_index=20,
                                   daily=True)

    cf_timeseries_daily_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries_daily.csv')
    assert (compare_csv_timeseries_files(cf_timeseries_daily_file, cf_timeseries_daily_file_solution, header=False))

    #if file is CF compliant, check write out timeseries
    cf_timeseries_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries.csv')
    with RAPIDDataset(cf_qout_file) as qout_nc:
        qout_nc.write_flows_to_csv(cf_timeseries_file,
                                   river_index=20)

    cf_timeseries_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries.csv')
    assert (compare_csv_timeseries_files(cf_timeseries_file, cf_timeseries_file_solution, header=False))

    #if file is CF compliant, you can write out daily average, filter by date, and use max mode
    cf_timeseries_daily_date_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries_daily_date.csv')

    with RAPIDDataset(cf_qout_file) as qout_nc:
        qout_nc.write_flows_to_csv(cf_timeseries_daily_date_file,
                                   river_id=75224,
                                   date_search_start=datetime(2002, 8, 31),
                                   date_search_end=datetime(2002, 8, 31, 23, 59, 59),
                                   daily=True,
                                   mode='max')

    cf_timeseries_daily_date_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries_daily_date.csv')
    assert (compare_csv_timeseries_files(cf_timeseries_daily_date_file, cf_timeseries_daily_date_file_solution, header=False))

    #if file is CF compliant, check write out timeseries and filter by date
    cf_timeseries_date_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries_date.csv')
    with RAPIDDataset(cf_qout_file) as qout_nc:
        qout_nc.write_flows_to_csv(cf_timeseries_date_file,
                                   date_search_start=datetime(2002, 8, 31),
                                   #date_search_end=None,
                                   river_id=75224)

    cf_timeseries_date_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries_date.csv')
    assert (compare_csv_timeseries_files(cf_timeseries_date_file, cf_timeseries_date_file_solution, header=False))

    remove_files(new_timeseries_file,
                 new_qout_file,
                 original_timeseries_file,
                 original_qout_file,
                 cf_timeseries_file,
                 cf_timeseries_date_file,
                 cf_timeseries_daily_file,
                 cf_timeseries_daily_date_file,
                 cf_qout_file,
                 )


def test_goodness_of_fit():
    """
    This tests the goodness of fit functions
    """
    print("TEST 14: TEST GOODNESS OF FIT FUNCTIONS")

    reach_id_file = os.path.join(INPUT_DATA_PATH, 'obs_reach_id.csv')
    observed_file = os.path.join(INPUT_DATA_PATH, 'obs_flow.csv')
    #using CF-compliant file
    cf_input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_CF.nc')
    cf_out_analysis_file = os.path.join(OUTPUT_DATA_PATH, 'cf_goodness_of_fit_results-daily.csv')
    find_goodness_of_fit(cf_input_qout_file, reach_id_file, observed_file,
                         cf_out_analysis_file, daily=True)

    cf_goodness_of_fit_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_goodness_of_fit_analysis.csv')
    assert (compare_csv_decimal_files(cf_out_analysis_file, cf_goodness_of_fit_file_solution))
    #using original RAPID file
    raw_goodness_of_fit_file_solution = os.path.join(COMPARE_DATA_PATH, 'raw_goodness_of_fit_analysis.csv')
    original_input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_original.nc')
    original_out_analysis_file = os.path.join(OUTPUT_DATA_PATH, 'original_goodness_of_fit_results-daily.csv')
    find_goodness_of_fit(original_input_qout_file, reach_id_file, observed_file,
                         original_out_analysis_file, steps_per_group=8)

    assert (compare_csv_decimal_files(original_out_analysis_file, raw_goodness_of_fit_file_solution))

    #using new RAPID file
    new_input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830.nc')
    new_out_analysis_file = os.path.join(OUTPUT_DATA_PATH, 'goodness_of_fit_results-daily.csv')
    find_goodness_of_fit(new_input_qout_file, reach_id_file, observed_file,
                         new_out_analysis_file, steps_per_group=8)

    assert (compare_csv_decimal_files(new_out_analysis_file, raw_goodness_of_fit_file_solution))

    reach_id_file = os.path.join(INPUT_DATA_PATH, 'obs_reach_id_1.csv')
    observed_file = os.path.join(INPUT_DATA_PATH, 'obs_flow_1.csv')
    #using CF-compliant file single input
    cf_out_analysis_file_1 = os.path.join(OUTPUT_DATA_PATH, 'cf_goodness_of_fit_results_1-daily.csv')
    find_goodness_of_fit(cf_input_qout_file, reach_id_file, observed_file,
                         cf_out_analysis_file_1, daily=True)

    cf_goodness_of_fit_file_solution_1 = os.path.join(COMPARE_DATA_PATH, 'cf_goodness_of_fit_analysis_1.csv')
    assert (compare_csv_decimal_files(cf_out_analysis_file_1, cf_goodness_of_fit_file_solution_1))

    remove_files(cf_out_analysis_file,
                 original_out_analysis_file,
                 new_out_analysis_file,
                 cf_out_analysis_file_1)

def test_cf_merge():
    """
    This tests merging two qout files
    """
    print("TEST 15: TEST MERGE QOUT")

    orig_qout_1 = os.path.join(INPUT_DATA_PATH, 'Qout_merge_3hr.nc')
    orig_qout_2 = os.path.join(INPUT_DATA_PATH, 'Qout_merge_6hr.nc')

    qout_1 = os.path.join(OUTPUT_DATA_PATH, 'Qout_merge_3hr.nc')
    qout_2 = os.path.join(OUTPUT_DATA_PATH, 'Qout_merge_6hr.nc')

    copy(orig_qout_1, qout_1)
    copy(orig_qout_2, qout_2)
    #Merge all files together at the end
    cv = ConvertRAPIDOutputToCF(rapid_output_file=[qout_1, qout_2],
                                start_datetime=datetime(2016, 2, 12),
                                time_step=[3*3600, 6*3600],
                                qinit_file="",
                                comid_lat_lon_z_file="",
                                rapid_connect_file="",
                                project_name="ECMWF-RAPID Predicted flows by US Army ERDC",
                                output_id_dim_name='rivid',
                                output_flow_var_name='Qout',
                                print_debug=False)
    cv.convert()

    cf_merge_qout_file_solution = os.path.join(COMPARE_DATA_PATH,
                                               'Qout_merge.nc')

    #check Qout
    assert (compare_qout_files(qout_1, cf_merge_qout_file_solution))

    #check other info in netcdf file
    d1 = Dataset(qout_1)
    d2 = Dataset(cf_merge_qout_file_solution)
    assert (d1.dimensions.keys() == d2.dimensions.keys())
    assert (d1.variables.keys() == d2.variables.keys())
    assert ((d1.variables['time'][:] == d2.variables['time'][:]).all())
    assert ((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    d1.close()
    d2.close()

    remove_files(qout_1,
                 qout_2)

def test_extract_timeseries_to_gssha_xys():
    """
    This tests extracting a timeseries from RAPID Qout file to GSHHA xys file
    """
    print("TEST 16: TEST EXTRACT TIMESERIES FROM Qout file to GSSHA xys file")

    cf_input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_CF.nc')
    cf_qout_file = os.path.join(OUTPUT_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_CF.nc')
    copy(cf_input_qout_file, cf_qout_file)

    #if file is CF compliant, you can write out daily average
    cf_timeseries_daily_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries_daily.xys')

    with RAPIDDataset(cf_qout_file) as qout_nc:
        qout_nc.write_flows_to_gssha_time_series_xys(cf_timeseries_daily_file,
                                                     series_name="RAPID_TO_GSSHA",
                                                     series_id=25,
                                                     river_index=20,
                                                     daily=True)

    cf_timeseries_daily_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries_daily.xys')
    assert (compare_csv_timeseries_files(cf_timeseries_daily_file, cf_timeseries_daily_file_solution))

    #if file is CF compliant, check write out timeseries
    cf_timeseries_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries.xys')
    with RAPIDDataset(cf_qout_file) as qout_nc:
        qout_nc.write_flows_to_gssha_time_series_xys(cf_timeseries_file,
                                                     series_name="RAPID_TO_GSSHA",
                                                     series_id=25,
                                                     river_index=20)

    cf_timeseries_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries.xys')
    assert (compare_csv_timeseries_files(cf_timeseries_file, cf_timeseries_file_solution, header=True))

    #if file is CF compliant, you can write out daily average, filter by date, and use max mode
    cf_timeseries_daily_date_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries_daily_date.xys')

    with RAPIDDataset(cf_qout_file) as qout_nc:
        qout_nc.write_flows_to_gssha_time_series_xys(cf_timeseries_daily_date_file,
                                                     series_name="RAPID_TO_GSSHA",
                                                     series_id=25,
                                                     river_id=75224,
                                                     date_search_start=datetime(2002, 8, 31),
                                                     date_search_end=datetime(2002, 8, 31, 23, 59, 59),
                                                     daily=True,
                                                     mode='max')

    cf_timeseries_daily_date_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries_daily_date.xys')
    assert (compare_csv_timeseries_files(cf_timeseries_daily_date_file, cf_timeseries_daily_date_file_solution))

    #if file is CF compliant, check write out timeseries and filter by date
    cf_timeseries_date_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries_date.xys')
    with RAPIDDataset(cf_qout_file) as qout_nc:
        qout_nc.write_flows_to_gssha_time_series_xys(cf_timeseries_date_file,
                                                     series_name="RAPID_TO_GSSHA",
                                                     series_id=25,
                                                     date_search_start=datetime(2002, 8, 31),
                                                     #date_search_end=None,
                                                     river_id=75224)

    cf_timeseries_date_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries_date.xys')
    assert (compare_csv_timeseries_files(cf_timeseries_date_file, cf_timeseries_date_file_solution))

    remove_files(cf_timeseries_file,
                 cf_qout_file,
                 cf_timeseries_daily_file,
                 cf_timeseries_daily_date_file,
                 cf_timeseries_date_file,
                 )

def test_extract_timeseries_to_gssha_ihg():
    """
    This tests extracting a timeseries from RAPID Qout file to GSHHA ihg file
    """
    print("TEST 16: TEST EXTRACT TIMESERIES FROM Qout file to GSSHA ihg file")

    cf_input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_CF.nc')
    cf_qout_file = os.path.join(OUTPUT_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_CF.nc')
    copy(cf_input_qout_file, cf_qout_file)


    #if file is CF compliant, you can write out daily average
    connection_list_file = os.path.join(INPUT_DATA_PATH, 'rapid_gssha_connect_file3.csv')
    cf_timeseries_daily_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries_daily.ihg')

    with RAPIDDataset(cf_qout_file) as qout_nc:
        qout_nc.write_flows_to_gssha_time_series_ihg(cf_timeseries_daily_file,
                                                     connection_list_file=connection_list_file,
                                                     daily=True)

    cf_timeseries_daily_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries_daily.ihg')
    assert (compare_csv_timeseries_files(cf_timeseries_daily_file, cf_timeseries_daily_file_solution, header=False))

    #if file is CF compliant, check write out timeseries
    connection_list_file = os.path.join(INPUT_DATA_PATH, 'rapid_gssha_connect_file1.csv')
    cf_timeseries_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries.ihg')
    with RAPIDDataset(cf_qout_file) as qout_nc:
        qout_nc.write_flows_to_gssha_time_series_ihg(cf_timeseries_file,
                                                     connection_list_file=connection_list_file,
                                                     )

    cf_timeseries_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries.ihg')
    assert (compare_csv_timeseries_files(cf_timeseries_file, cf_timeseries_file_solution, header=False))

    #if file is CF compliant, you can write out daily average, filter by date, and use max mode
    connection_list_file = os.path.join(INPUT_DATA_PATH, 'rapid_gssha_connect_file1.csv')
    cf_timeseries_daily_date_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries_daily_date.ihg')

    with RAPIDDataset(cf_qout_file) as qout_nc:
        qout_nc.write_flows_to_gssha_time_series_ihg(cf_timeseries_daily_date_file,
                                                     connection_list_file=connection_list_file,
                                                     date_search_start=datetime(2002, 8, 31),
                                                     date_search_end=datetime(2002, 8, 31, 23, 59, 59),
                                                     daily=True,
                                                     mode='max')

    cf_timeseries_daily_date_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries_daily_date.ihg')
    assert (compare_csv_timeseries_files(cf_timeseries_daily_date_file, cf_timeseries_daily_date_file_solution, header=False))

    #if file is CF compliant, check write out timeseries and filter by date
    connection_list_file = os.path.join(INPUT_DATA_PATH, 'rapid_gssha_connect_file3.csv')
    cf_timeseries_date_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries_date.ihg')
    with RAPIDDataset(cf_qout_file) as qout_nc:
        qout_nc.write_flows_to_gssha_time_series_ihg(cf_timeseries_date_file,
                                                     connection_list_file=connection_list_file,
                                                     date_search_start=datetime(2002, 8, 31),
                                                     #date_search_end=None,
                                                     )

    cf_timeseries_date_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries_date.ihg')
    assert (compare_csv_timeseries_files(cf_timeseries_date_file, cf_timeseries_date_file_solution, header=False))

    remove_files(cf_timeseries_file,
                 cf_qout_file,
                 cf_timeseries_daily_file,
                 cf_timeseries_daily_date_file,
                 cf_timeseries_date_file,
                 )

def test_extract_timeseries_to_gssha_ihg_tzinfo():
    """
    This tests extracting a timeseries from RAPID Qout file to GSHHA ihg file
    with different time zone output
    """
    print("TEST 17: TEST EXTRACT TIMESERIES FROM Qout file to GSSHA ihg file tzinfo")

    CENTRAL_TZ = timezone('US/Central')

    cf_input_qout_file = os.path.join(COMPARE_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_CF.nc')
    cf_qout_file = os.path.join(OUTPUT_DATA_PATH, 'Qout_nasa_lis_3hr_20020830_CF.nc')
    copy(cf_input_qout_file, cf_qout_file)


    #if file is CF compliant, you can write out daily average
    connection_list_file = os.path.join(INPUT_DATA_PATH, 'rapid_gssha_connect_file3.csv')
    cf_timeseries_daily_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries_daily_tz.ihg')

    with RAPIDDataset(cf_qout_file, out_tzinfo=CENTRAL_TZ) as qout_nc:
        qout_nc.write_flows_to_gssha_time_series_ihg(cf_timeseries_daily_file,
                                                     connection_list_file=connection_list_file,
                                                     daily=True)

    cf_timeseries_daily_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries_daily_tz.ihg')
    assert (compare_csv_timeseries_files(cf_timeseries_daily_file, cf_timeseries_daily_file_solution, header=False))

    #if file is CF compliant, check write out timeseries
    connection_list_file = os.path.join(INPUT_DATA_PATH, 'rapid_gssha_connect_file1.csv')
    cf_timeseries_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries_tz.ihg')
    with RAPIDDataset(cf_qout_file, out_tzinfo=CENTRAL_TZ) as qout_nc:
        qout_nc.write_flows_to_gssha_time_series_ihg(cf_timeseries_file,
                                                     connection_list_file=connection_list_file,
                                                     )

    cf_timeseries_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries_tz.ihg')
    assert (compare_csv_timeseries_files(cf_timeseries_file, cf_timeseries_file_solution, header=False))

    #if file is CF compliant, you can write out daily average, filter by date, and use max mode
    connection_list_file = os.path.join(INPUT_DATA_PATH, 'rapid_gssha_connect_file1.csv')
    cf_timeseries_daily_date_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries_daily_date_tz.ihg')

    with RAPIDDataset(cf_qout_file, out_tzinfo=CENTRAL_TZ) as qout_nc:
        qout_nc.write_flows_to_gssha_time_series_ihg(cf_timeseries_daily_date_file,
                                                     connection_list_file=connection_list_file,
                                                     date_search_start=datetime(2002, 8, 31),
                                                     date_search_end=datetime(2002, 8, 31, 23, 59, 59),
                                                     daily=True,
                                                     mode='max')

    cf_timeseries_daily_date_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries_daily_date_tz.ihg')
    assert (compare_csv_timeseries_files(cf_timeseries_daily_date_file, cf_timeseries_daily_date_file_solution, header=False))

    #if file is CF compliant, check write out timeseries and filter by date
    connection_list_file = os.path.join(INPUT_DATA_PATH, 'rapid_gssha_connect_file3.csv')
    cf_timeseries_date_file = os.path.join(OUTPUT_DATA_PATH, 'cf_timeseries_date_tz.ihg')
    with RAPIDDataset(cf_qout_file, out_tzinfo=CENTRAL_TZ) as qout_nc:
        qout_nc.write_flows_to_gssha_time_series_ihg(cf_timeseries_date_file,
                                                     connection_list_file=connection_list_file,
                                                     date_search_start=datetime(2002, 8, 31),
                                                     )

    cf_timeseries_date_file_solution = os.path.join(COMPARE_DATA_PATH, 'cf_timeseries_date_tz.ihg')
    assert (compare_csv_timeseries_files(cf_timeseries_date_file, cf_timeseries_date_file_solution, header=False))

    remove_files(cf_timeseries_file,
                 cf_qout_file,
                 cf_timeseries_daily_file,
                 cf_timeseries_daily_date_file,
                 cf_timeseries_date_file,
                 )
