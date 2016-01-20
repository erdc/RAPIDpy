# -*- coding: utf-8 -*-
##
##  test_rapidpy.py
##  RAPIDpy
##
##  Created by Alan D. Snow 2016.
##  Copyright © 2015 Alan D Snow. All rights reserved.
##

from filecmp import cmp as fcmp
from nose.tools import raises, ok_
import os
from shutil import copy


#local import
import sys
from RAPIDpy.rapid import RAPID


#GLOBAL VARIABLES
MAIN_TESTS_FOLDER = os.path.dirname(os.path.abspath(__file__))
COMPARE_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'compare_data')
INPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'input_data')
OUTPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'tmp_output_files')

#------------------------------------------------------------------------------
# MAIN TEST SCRIPTS
#------------------------------------------------------------------------------
def test_generate_rapid_input_file():
    """
    Checks RAPID input file generation with valid input
    """
    print "TEST 1: GENERATE NAMELIST FILE"
    rapid_manager = RAPID(rapid_executable_location="",
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
                                                 "rapid_namelist-GENERATE-SOLUTION")


    ok_(fcmp(generated_input_file, generated_input_file_solution))
    
    try:
        os.remove(generated_input_file)
    except OSError:
        pass

def test_update_rapid_input_file():
    """
    Checks RAPID input file update with valid input
    """
    print "TEST 2: UPDATE NAMELIST FILE"
    rapid_manager = RAPID(rapid_executable_location="",
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
                                               "rapid_namelist-UPDATE-SOLUTION")


    ok_(fcmp(updated_input_file, updated_input_file_solution))
    
    try:
        os.remove(updated_input_file)
    except OSError:
        pass

def test_update_rapid_invalid_input_file():
    """
    Checks RAPID input file update with valid input
    """
    print "TEST 3: UPDATE WITH INVALID NAMELIST FILE"
    rapid_manager = RAPID(rapid_executable_location="",
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
                                               "rapid_namelist-UPDATE-INVALID-SOLUTION")


    ok_(fcmp(updated_input_file, updated_input_file_solution))
    
    try:
        os.remove(updated_input_file)
    except OSError:
        pass

def test_update_rapid_numbers_input_file():
    """
    Checks RAPID input file update with number validation
    """
    rapid_input_data_path = os.path.join(INPUT_DATA_PATH,
                                         'nfie_texas_gulf_region-huc_2_12')

    print "TEST 4: GENERATE NUMBERS FOR NAMELIST FILE"
    rapid_manager = RAPID(rapid_executable_location="",
                          use_all_processors=True,
                          rapid_connect_file=os.path.join(rapid_input_data_path, 'rapid_connect.csv'),
                          riv_bas_id_file=os.path.join(rapid_input_data_path, 'riv_bas_id.csv'),
                         )
    rapid_manager.update_reach_number_data()
                          
    rapid_manager.update_parameters(rapid_connect_file='rapid_connect.csv',
                                    Vlat_file='m3_riv_bas_erai_t511_3hr_19800101to19800101.nc',
                                    riv_bas_id_file='riv_bas_id.csv',
                                    k_file='k.csv',
                                    x_file='x.csv',
                                    Qout_file='Qout.nc'
                                    )

    generated_input_file = os.path.join(OUTPUT_DATA_PATH, 
                                      "rapid_namelist-GENERATE-NUMBERS")

    rapid_manager.generate_namelist_file(generated_input_file)
                          
    generated_input_file_solution = os.path.join(COMPARE_DATA_PATH, 
                                               "rapid_namelist-GENERATE-NUMBERS-SOLUTION")


    ok_(fcmp(generated_input_file, generated_input_file_solution))
    
    try:
        os.remove(generated_input_file)
    except OSError:
        pass

def test_update_rapid_numbers_input_file():
    """
    Test Running RAPID Simulation
    """
    
    rapid_executable_location = os.path.join(MAIN_TESTS_FOLDER,
                                             "..", "..",
                                             "rapid", "src", "rapid")
                                             
    rapid_input_data_path = os.path.join(INPUT_DATA_PATH,
                                         'nfie_texas_gulf_region-huc_2_12')

    generated_qout_file = os.path.join(OUTPUT_DATA_PATH, 'Qout_erai_t511_3hr_19800101.nc')



    print "TEST 5: TEST RUNNING RAPID SIMULATION"
    rapid_manager = RAPID(rapid_executable_location,
                          use_all_processors=True,
                          rapid_connect_file=os.path.join(rapid_input_data_path, 'rapid_connect.csv'),
                          riv_bas_id_file=os.path.join(rapid_input_data_path, 'riv_bas_id.csv'),
                          Vlat_file=os.path.join(rapid_input_data_path, 'm3_riv_bas_erai_t511_3hr_19800101.nc'),
                          k_file=os.path.join(rapid_input_data_path, 'k.csv'),
                          x_file=os.path.join(rapid_input_data_path, 'x.csv'),
                          Qout_file=generated_qout_file
                         )

    rapid_manager.run()

    generated_qout_file_solution = os.path.join(COMPARE_DATA_PATH,
                                                'Qout_erai_t511_3hr_19800101.nc')

    ok_(fcmp(generated_qout_file, generated_qout_file_solution))
    
    try:
        os.remove(generated_qout_file)
    except OSError:
        pass

if __name__ == '__main__':
    import nose
    nose.main()

