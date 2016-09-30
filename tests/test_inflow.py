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
from shutil import copytree, rmtree

#local import
from RAPIDpy.inflow import run_lsm_rapid_process
from RAPIDpy.helper_functions import (compare_csv_decimal_files,
                                      remove_files)
#GLOBAL VARIABLES
MAIN_TESTS_FOLDER = os.path.dirname(os.path.abspath(__file__))
COMPARE_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'compare')
INFLOW_COMPARE_DATA_PATH = os.path.join(COMPARE_DATA_PATH, 'inflow')
LSM_INPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'data','lsm_grids')
OUTPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'output')
RAPID_DATA_PATH = os.path.join(OUTPUT_DATA_PATH, 'input')
RAPID_EXE_PATH = os.path.join(MAIN_TESTS_FOLDER,
                              "..", "..",
                              "rapid", "src", "rapid")
CYGWIN_BIN_PATH = 'C:\\cygwin64\\bin'

#------------------------------------------------------------------------------
# MAIN TEST SCRIPTS
#------------------------------------------------------------------------------
def test_run_era_interim_inflow():
    """
    Checks generating inflow file from ERA Interim LSM
    """
    print("TEST 1: TEST GENERATE INFLOW FILE FROM ERA INTERIM DATA")
    
    rapid_input_path = os.path.join(RAPID_DATA_PATH, "x-x")
    rapid_output_path = os.path.join(OUTPUT_DATA_PATH, "output", "x-x")
    #Create testing environment
    try:
        os.mkdir(RAPID_DATA_PATH)
    except OSError:
        pass
    try:
        os.mkdir(os.path.join(OUTPUT_DATA_PATH, "output"))
    except OSError:
        pass
    
    try:
        copytree(os.path.join(COMPARE_DATA_PATH, "gis","x-x"),
                 rapid_input_path)    
    except OSError:
        pass
    
    #run main process    
    run_lsm_rapid_process(
        rapid_executable_location=RAPID_EXE_PATH,
        cygwin_bin_location=CYGWIN_BIN_PATH,
        rapid_io_files_location=OUTPUT_DATA_PATH,
        lsm_data_location=os.path.join(LSM_INPUT_DATA_PATH, 'erai3'), 
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
    m3_file_name = "m3_riv_bas_erai_t511_3hr_20030121to20030122.nc"
    generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
    generated_m3_file_solution = os.path.join(INFLOW_COMPARE_DATA_PATH, m3_file_name)
    #check other info in netcdf file
    d1 = Dataset(generated_m3_file)
    d2 = Dataset(generated_m3_file_solution)
    assert_almost_equal(d1.variables['m3_riv'][:], d2.variables['m3_riv'][:], decimal=5)
    if 'rivid' in d2.variables.keys():
        ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    if 'lat' in d2.variables.keys():
        ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
    if 'lon' in d2.variables.keys():
        ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
    d1.close()
    d2.close()
    
    #qout file
    qout_file_name = "Qout_erai_t511_3hr_20030121to20030122.nc"
    generated_qout_file = os.path.join(rapid_output_path, qout_file_name)
    generated_qout_file_solution = os.path.join(INFLOW_COMPARE_DATA_PATH, qout_file_name)
    d1 = Dataset(generated_qout_file)
    d2 = Dataset(generated_qout_file_solution)
    assert_almost_equal(d1.variables['Qout'][:], d2.variables['Qout'][:], decimal=5)
    ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    ok_((d1.variables['time'][:] == d2.variables['time'][:]).all())
    if 'lat' in d2.variables.keys():
        ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
    if 'lon' in d2.variables.keys():
        ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
    d1.close()
    d2.close()
                                                 
    #initialization file
    qinit_file_name = "qinit_erai_t511_3hr_20030121to20030122.csv"
    generated_qinit_file = os.path.join(rapid_input_path, qinit_file_name)
    generated_qinit_file_solution = os.path.join(INFLOW_COMPARE_DATA_PATH, qinit_file_name)

    ok_(compare_csv_decimal_files(generated_qinit_file, generated_qinit_file_solution))
    
    #cleanup
    remove_files(generated_qinit_file)
    rmtree(os.path.join(OUTPUT_DATA_PATH, "input"))
    rmtree(os.path.join(OUTPUT_DATA_PATH, "output"))

def test_generate_nldas2_inflow():
    """
    Checks generating inflow file from NLDAS V2 LSM
    """
    print("TEST 2: TEST GENERATE INFLOW FILE FROM NLDAS V2 DATA")
    rapid_output_path = os.path.join(OUTPUT_DATA_PATH, "output", "x-x")
    #Create testing environment
    try:
        os.mkdir(RAPID_DATA_PATH)
    except OSError:
        pass
    try:
        os.mkdir(os.path.join(OUTPUT_DATA_PATH, "output"))
    except OSError:
        pass
    
    try:
        copytree(os.path.join(COMPARE_DATA_PATH, "gis","x-x"),
                 os.path.join(RAPID_DATA_PATH, "x-x"))    
    except OSError:
        pass
    
    #run main process    
    run_lsm_rapid_process(
        rapid_executable_location=RAPID_EXE_PATH,
        cygwin_bin_location=CYGWIN_BIN_PATH,
        rapid_io_files_location=OUTPUT_DATA_PATH,
        lsm_data_location=os.path.join(LSM_INPUT_DATA_PATH, 'nldas2'), 
        simulation_start_datetime=datetime(1980, 1, 1),
        simulation_end_datetime=datetime(2014, 1, 31),
        generate_rapid_namelist_file=False,
        run_rapid_simulation=False,
        use_all_processors=True,
    )
    
    #CHECK OUTPUT    
    #m3_riv
    m3_file_name = "m3_riv_bas_nasa_nldas_3hr_20030121to20030121.nc"
    generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
    generated_m3_file_solution = os.path.join(INFLOW_COMPARE_DATA_PATH, m3_file_name)
    #check other info in netcdf file
    d1 = Dataset(generated_m3_file)
    d2 = Dataset(generated_m3_file_solution)
    assert_almost_equal(d1.variables['m3_riv'][:], d2.variables['m3_riv'][:], decimal=5)
    if 'rivid' in d2.variables.keys():
        ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    if 'lat' in d2.variables.keys():
        ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
    if 'lon' in d2.variables.keys():
        ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
    d1.close()
    d2.close()

    #cleanup
    rmtree(os.path.join(OUTPUT_DATA_PATH, "input"))
    rmtree(os.path.join(OUTPUT_DATA_PATH, "output"))

def test_generate_era20cm_inflow():
    """
    Checks generating inflow file from ERA 20CM LSM
    """
    print("TEST 3: TEST GENERATE INFLOW FILE FROM ERA 20CM DATA")

    rapid_output_path = os.path.join(OUTPUT_DATA_PATH, "output", "x-x")
    #Create testing environment
    try:
        os.mkdir(RAPID_DATA_PATH)
    except OSError:
        pass
    try:
        os.mkdir(os.path.join(OUTPUT_DATA_PATH, "output"))
    except OSError:
        pass
    
    try:
        copytree(os.path.join(COMPARE_DATA_PATH, "gis","x-x"),
                 os.path.join(RAPID_DATA_PATH, "x-x"))    
    except OSError:
        pass

    run_lsm_rapid_process(
        rapid_executable_location=RAPID_EXE_PATH,
        cygwin_bin_location=CYGWIN_BIN_PATH,
        rapid_io_files_location=OUTPUT_DATA_PATH,
        lsm_data_location=os.path.join(LSM_INPUT_DATA_PATH, 'era20cm'), 
        simulation_start_datetime=datetime(1980, 1, 1),
        simulation_end_datetime=datetime(2014, 1, 31),
        ensemble_list=range(10),
        generate_rapid_namelist_file=False,
        run_rapid_simulation=False,
        use_all_processors=True,
    )
    
    for i in range(10):
        #CHECK OUTPUT    
        #m3_riv
        m3_file_name = "m3_riv_bas_era_20cm_t159_3hr_20000129to20000130_{0}.nc".format(i)
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
        generated_m3_file_solution = os.path.join(INFLOW_COMPARE_DATA_PATH, m3_file_name)

        #check other info in netcdf file
        d1 = Dataset(generated_m3_file)
        d2 = Dataset(generated_m3_file_solution)
        assert_almost_equal(d1.variables['m3_riv'][:], d2.variables['m3_riv'][:], decimal=5)
        if 'rivid' in d2.variables.keys():
            ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
        if 'lat' in d2.variables.keys():
            ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
        if 'lon' in d2.variables.keys():
            ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
        d1.close()
        d2.close()

    #cleanup
    rmtree(os.path.join(OUTPUT_DATA_PATH, "input"))
    rmtree(os.path.join(OUTPUT_DATA_PATH, "output"))

def test_generate_erai_t255_inflow():
    """
    Checks generating inflow file from ERA Interim t255 LSM
    """
    print("TEST 4: TEST GENERATE INFLOW FILE FROM ERA Interim t255 DATA")
    rapid_output_path = os.path.join(OUTPUT_DATA_PATH, "output", "x-x")
    #Create testing environment
    try:
        os.mkdir(RAPID_DATA_PATH)
    except OSError:
        pass
    try:
        os.mkdir(os.path.join(OUTPUT_DATA_PATH, "output"))
    except OSError:
        pass
    
    try:
        copytree(os.path.join(COMPARE_DATA_PATH, "gis","x-x"),
                 os.path.join(RAPID_DATA_PATH, "x-x"))    
    except OSError:
        pass

    #run main process
    run_lsm_rapid_process(
        rapid_executable_location=RAPID_EXE_PATH,
        cygwin_bin_location=CYGWIN_BIN_PATH,
        rapid_io_files_location=OUTPUT_DATA_PATH,
        lsm_data_location=os.path.join(LSM_INPUT_DATA_PATH, 'erai3t255'), 
        simulation_start_datetime=datetime(1980, 1, 1),
        simulation_end_datetime=datetime(2014, 12, 31),
        generate_rapid_namelist_file=False,
        run_rapid_simulation=False,
        use_all_processors=True,
    )
    
    #CHECK OUTPUT    
    #m3_riv
    m3_file_name = "m3_riv_bas_erai_t255_3hr_20140820to20140821.nc"
    generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
    generated_m3_file_solution = os.path.join(INFLOW_COMPARE_DATA_PATH, m3_file_name)
    #check other info in netcdf file
    d1 = Dataset(generated_m3_file)
    d2 = Dataset(generated_m3_file_solution)
    assert_almost_equal(d1.variables['m3_riv'][:], d2.variables['m3_riv'][:], decimal=5)
    if 'rivid' in d2.variables.keys():
        ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    if 'lat' in d2.variables.keys():
        ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
    if 'lon' in d2.variables.keys():
        ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
    d1.close()
    d2.close()

    #cleanup
    rmtree(os.path.join(OUTPUT_DATA_PATH, "input"))
    rmtree(os.path.join(OUTPUT_DATA_PATH, "output"))
    
def test_generate_gldas2_inflow():
    """
    Checks generating inflow file from GLDAS V2 LSM
    """
    print("TEST 5: TEST GENERATE INFLOW FILE FROM GLDAS V2 DATA")
    rapid_output_path = os.path.join(OUTPUT_DATA_PATH, "output", "x-x")
    #Create testing environment
    try:
        os.mkdir(RAPID_DATA_PATH)
    except OSError:
        pass
    try:
        os.mkdir(os.path.join(OUTPUT_DATA_PATH, "output"))
    except OSError:
        pass
    
    try:
        copytree(os.path.join(COMPARE_DATA_PATH, "gis","x-x"),
                 os.path.join(RAPID_DATA_PATH, "x-x"))    
    except OSError:
        pass
    
    #run main process    
    run_lsm_rapid_process(
        rapid_executable_location=RAPID_EXE_PATH,
        cygwin_bin_location=CYGWIN_BIN_PATH,
        rapid_io_files_location=OUTPUT_DATA_PATH,
        lsm_data_location=os.path.join(LSM_INPUT_DATA_PATH, 'gldas2'), 
        simulation_start_datetime=datetime(1980, 1, 1),
        simulation_end_datetime=datetime(2014, 1, 31),
        generate_rapid_namelist_file=False,
        run_rapid_simulation=False,
        use_all_processors=True,
    )
    
    #CHECK OUTPUT    
    #m3_riv
    m3_file_name = "m3_riv_bas_nasa_gldas2_3hr_20101231to20101231.nc"
    generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
    generated_m3_file_solution = os.path.join(INFLOW_COMPARE_DATA_PATH, m3_file_name)
    #check other info in netcdf file
    d1 = Dataset(generated_m3_file)
    d2 = Dataset(generated_m3_file_solution)
    assert_almost_equal(d1.variables['m3_riv'][:], d2.variables['m3_riv'][:], decimal=5)
    if 'rivid' in d2.variables.keys():
        ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    if 'lat' in d2.variables.keys():
        ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
    if 'lon' in d2.variables.keys():
        ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
    d1.close()
    d2.close()

    #cleanup
    rmtree(os.path.join(OUTPUT_DATA_PATH, "input"))
    rmtree(os.path.join(OUTPUT_DATA_PATH, "output"))
    
def test_generate_lis_inflow():
    """
    Checks generating inflow file from LIS LSM
    """
    print("TEST 5: TEST GENERATE INFLOW FILE FROM LIS DATA")
    rapid_output_path = os.path.join(OUTPUT_DATA_PATH, "output", "u-k")
    #Create testing environment
    try:
        os.mkdir(RAPID_DATA_PATH)
    except OSError:
        pass
    try:
        os.mkdir(os.path.join(OUTPUT_DATA_PATH, "output"))
    except OSError:
        pass
    
    try:
        copytree(os.path.join(COMPARE_DATA_PATH, "gis","u-k"),
                 os.path.join(RAPID_DATA_PATH, "u-k"))    
    except OSError:
        pass
    
    #run main process    
    run_lsm_rapid_process(
        rapid_executable_location=RAPID_EXE_PATH,
        cygwin_bin_location=CYGWIN_BIN_PATH,
        rapid_io_files_location=OUTPUT_DATA_PATH,
        lsm_data_location=os.path.join(LSM_INPUT_DATA_PATH, 'lis'), 
        simulation_start_datetime=datetime(1980, 1, 1),
        simulation_end_datetime=datetime(2014, 1, 31),
        generate_rapid_namelist_file=False,
        run_rapid_simulation=False,
        use_all_processors=True,
    )
    
    #CHECK OUTPUT
    m3_file_name = "m3_riv_bas_nasa_lis_3hr_20110121to20110121.nc"
    generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
    generated_m3_file_solution = os.path.join(INFLOW_COMPARE_DATA_PATH, m3_file_name)
    #check other info in netcdf file
    d1 = Dataset(generated_m3_file)
    d2 = Dataset(generated_m3_file_solution)
    assert_almost_equal(d1.variables['m3_riv'][:], d2.variables['m3_riv'][:], decimal=5)
    if 'rivid' in d2.variables.keys():
        ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    if 'lat' in d2.variables.keys():
        ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
    if 'lon' in d2.variables.keys():
        ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
    d1.close()
    d2.close()

    #cleanup
    rmtree(os.path.join(OUTPUT_DATA_PATH, "input"))
    rmtree(os.path.join(OUTPUT_DATA_PATH, "output"))

def test_generate_joules_inflow():
    """
    Checks generating inflow file from Joules LSM
    """
    print("TEST 5: TEST GENERATE INFLOW FILE FROM Joules DATA")
    rapid_output_path = os.path.join(OUTPUT_DATA_PATH, "output", "u-k")
    #Create testing environment
    try:
        os.mkdir(RAPID_DATA_PATH)
    except OSError:
        pass
    try:
        os.mkdir(os.path.join(OUTPUT_DATA_PATH, "output"))
    except OSError:
        pass
    
    try:
        copytree(os.path.join(COMPARE_DATA_PATH, "gis","u-k"),
                 os.path.join(RAPID_DATA_PATH, "u-k"))    
    except OSError:
        pass
    
    #run main process    
    run_lsm_rapid_process(
        rapid_executable_location=RAPID_EXE_PATH,
        cygwin_bin_location=CYGWIN_BIN_PATH,
        rapid_io_files_location=OUTPUT_DATA_PATH,
        lsm_data_location=os.path.join(LSM_INPUT_DATA_PATH, 'joules'), 
        simulation_start_datetime=datetime(1980, 1, 1),
        simulation_end_datetime=datetime(2014, 1, 31),
        file_datetime_re_pattern = r'\d{8}_\d{2}',
        file_datetime_pattern = "%Y%m%d_%H",      
        generate_rapid_namelist_file=False,
        run_rapid_simulation=False,
        use_all_processors=True,
    )
    
    #CHECK OUTPUT
    m3_file_name = "m3_riv_bas_met_office_joules_3hr_20080803to20080803.nc"
    generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
    generated_m3_file_solution = os.path.join(INFLOW_COMPARE_DATA_PATH, m3_file_name)
    #check other info in netcdf file
    d1 = Dataset(generated_m3_file)
    d2 = Dataset(generated_m3_file_solution)
    assert_almost_equal(d1.variables['m3_riv'][:], d2.variables['m3_riv'][:], decimal=5)
    if 'rivid' in d2.variables.keys():
        ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
    if 'lat' in d2.variables.keys():
        ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
    if 'lon' in d2.variables.keys():
        ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
    d1.close()
    d2.close()

    #cleanup
    rmtree(os.path.join(OUTPUT_DATA_PATH, "input"))
    rmtree(os.path.join(OUTPUT_DATA_PATH, "output"))

if __name__ == '__main__':
    import nose
    nose.main()
