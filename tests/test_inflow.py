# -*- coding: utf-8 -*-
##
##  test_inflow.py
##  RAPIDpy
##
##  Created by Alan D. Snow.
##  Copyright © 2016 Alan D Snow. All rights reserved.
##

from datetime import datetime
from netCDF4 import Dataset
from nose.tools import ok_
from numpy.testing import assert_almost_equal
import os
from shutil import rmtree


#local import
from RAPIDpy.inflow import run_lsm_rapid_process
from RAPIDpy.helper_functions import (compare_csv_decimal_files,
                                      remove_files)
#GLOBAL VARIABLES
MAIN_TESTS_FOLDER = os.path.dirname(os.path.abspath(__file__))
COMPARE_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'compare')
INPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'data')
RAPID_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'input', 'x-x')
OUTPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'output')
RAPID_EXE_PATH = os.path.join(MAIN_TESTS_FOLDER,
                              "..", "..",
                              "rapid", "src", "rapid")

#------------------------------------------------------------------------------
# MAIN TEST SCRIPTS
#------------------------------------------------------------------------------
def test_run_era_interim_inflow():
    """
    Checks generating inflow file from ERA Interim LSM
    """
    print "TEST 1: TEST GENERATE INFLOW FILE FROM ERA INTERIM DATA"
    run_lsm_rapid_process(
        rapid_executable_location=RAPID_EXE_PATH,
        rapid_io_files_location=MAIN_TESTS_FOLDER,
        lsm_data_location=os.path.join(INPUT_DATA_PATH, 'erai3'), #path to folder with LSM data
        simulation_start_datetime=datetime(1980, 1, 1),
        simulation_end_datetime=datetime(2014, 1, 31),
        generate_rapid_namelist_file=False,
        run_rapid_simulation=True,
        generate_return_periods_file=False,
        generate_seasonal_initialization_file=False,
        generate_initialization_file=True,
        use_all_processors=True,
    )
    
    #CHECK OUTPUT    
    #m3_riv
    generated_m3_file = os.path.join(OUTPUT_DATA_PATH, "x-x",
                                     "m3_riv_bas_erai_t511_3hr_20030121to20030122.nc")
    generated_m3_file_solution = os.path.join(COMPARE_DATA_PATH,
                                              "m3_riv_bas_erai_t511_3hr_20030121to20030122.nc")
    #check other info in netcdf file
    d1 = Dataset(generated_m3_file)
    d2 = Dataset(generated_m3_file_solution)
    assert_almost_equal(d1.variables['m3_riv'][:], d2.variables['m3_riv'][:])
    if 'rivid' in d2.variables.keys():
        ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    if 'lat' in d2.variables.keys():
        ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
    if 'lon' in d2.variables.keys():
        ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
    d1.close()
    d2.close()
    
    #qout file
    generated_qout_file = os.path.join(OUTPUT_DATA_PATH, "x-x",
                                       "Qout_erai_t511_3hr_20030121to20030122.nc")
    generated_qout_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                "Qout_erai_t511_3hr_20030121to20030122.nc")
    d1 = Dataset(generated_qout_file)
    d2 = Dataset(generated_qout_file_solution)
    assert_almost_equal(d1.variables['Qout'][:], d2.variables['Qout'][:])
    ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    ok_((d1.variables['time'][:] == d2.variables['time'][:]).all())
    if 'lat' in d2.variables.keys():
        ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
    if 'lon' in d2.variables.keys():
        ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
    d1.close()
    d2.close()
                                                 
    #initialization file
    generated_qinit_file = os.path.join(RAPID_DATA_PATH, 
                                        "qinit_erai_t511_3hr_20030121to20030122.csv")
    generated_qinit_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                "qinit_erai_t511_3hr_20030121to20030122.csv")
    ok_(compare_csv_decimal_files(generated_qinit_file, generated_qinit_file_solution))

    
    remove_files(generated_qinit_file)
    rmtree(os.path.join(OUTPUT_DATA_PATH, "x-x"))

def test_run_nldas2_inflow():
    """
    Checks generating inflow file from NLDAS V2 LSM
    """
    print "TEST 1: TEST GENERATE INFLOW FILE FROM NLDAS V2 DATA"
    run_lsm_rapid_process(
        rapid_executable_location=RAPID_EXE_PATH,
        rapid_io_files_location=MAIN_TESTS_FOLDER,
        lsm_data_location=os.path.join(INPUT_DATA_PATH, 'nldas2'), #path to folder with LSM data
        simulation_start_datetime=datetime(1980, 1, 1),
        simulation_end_datetime=datetime(2014, 1, 31),
        generate_rapid_namelist_file=False,
        run_rapid_simulation=False,
        use_all_processors=True,
    )
    
    #CHECK OUTPUT    
    #m3_riv
    generated_m3_file = os.path.join(OUTPUT_DATA_PATH, "x-x",
                                     "m3_riv_bas_nasa_nldas_3hr_20030121to20030121.nc")
    generated_m3_file_solution = os.path.join(COMPARE_DATA_PATH,
                                              "m3_riv_bas_nasa_nldas_3hr_20030121to20030121.nc")
    #check other info in netcdf file
    d1 = Dataset(generated_m3_file)
    d2 = Dataset(generated_m3_file_solution)
    assert_almost_equal(d1.variables['m3_riv'][:], d2.variables['m3_riv'][:])
    if 'rivid' in d2.variables.keys():
        ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    if 'lat' in d2.variables.keys():
        ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
    if 'lon' in d2.variables.keys():
        ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
    d1.close()
    d2.close()

    rmtree(os.path.join(OUTPUT_DATA_PATH, "x-x"))

if __name__ == '__main__':
    import nose
    nose.main()