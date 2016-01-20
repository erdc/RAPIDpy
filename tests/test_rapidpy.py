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

def test_generate_rapid_input_file():
    """
    Checks RAPID input file generation with valid input
    """
    main_tests_folder = os.path.dirname(os.path.abspath(__file__))
    
    compare_data_path = os.path.join(main_tests_folder, 'compare_data')
    output_data_path = os.path.join(main_tests_folder, 'tmp_output_files')

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
    generated_input_file = os.path.join(output_data_path, 
                                        "rapid_namelist-GENERATE")
    rapid_manager.generate_namelist_file(generated_input_file)
    generated_input_file_solution = os.path.join(compare_data_path, 
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
    main_tests_folder = os.path.dirname(os.path.abspath(__file__))
    
    compare_data_path = os.path.join(main_tests_folder, 'compare_data')
    input_data_path = os.path.join(main_tests_folder, 'input_data')
    output_data_path = os.path.join(main_tests_folder, 'tmp_output_files')

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

    original_input_file = os.path.join(input_data_path, 
                                      "rapid_namelist_valid")
    
    updated_input_file = os.path.join(output_data_path, 
                                      "rapid_namelist-UPDATE")

    copy(original_input_file, updated_input_file)
    rapid_manager.update_namelist_file(updated_input_file)
    updated_input_file_solution = os.path.join(compare_data_path, 
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
    main_tests_folder = os.path.dirname(os.path.abspath(__file__))
    
    compare_data_path = os.path.join(main_tests_folder, 'compare_data')
    input_data_path = os.path.join(main_tests_folder, 'input_data')
    output_data_path = os.path.join(main_tests_folder, 'tmp_output_files')

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

    original_input_file = os.path.join(input_data_path, 
                                      "rapid_namelist_invalid")
    
    updated_input_file = os.path.join(output_data_path, 
                                      "rapid_namelist-UPDATE-INVALID")

    copy(original_input_file, updated_input_file)
    rapid_manager.update_namelist_file(updated_input_file)
    updated_input_file_solution = os.path.join(compare_data_path, 
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
    main_tests_folder = os.path.dirname(os.path.abspath(__file__))
    
    compare_data_path = os.path.join(main_tests_folder, 'compare_data')
    input_data_path = os.path.join(main_tests_folder, 'input_data',
                                   'nfie_texas_gulf_region-huc_2_12')
    output_data_path = os.path.join(main_tests_folder, 'tmp_output_files')

    print "TEST 2: UPDATE NAMELIST FILE"
    rapid_manager = RAPID(rapid_executable_location="",
                          use_all_processors=True,
                          rapid_connect_file=os.path.join(input_data_path, 'rapid_connect.csv'),
                          riv_bas_id_file=os.path.join(input_data_path, 'riv_bas_id.csv'),
                         )
    rapid_manager.update_reach_number_data()
                          
    rapid_manager.update_parameters(rapid_connect_file='rapid_connect.csv',
                                    Vlat_file='m3_riv_bas_erai_t511_3hr_19800101to19800101.nc',
                                    riv_bas_id_file='riv_bas_id.csv',
                                    k_file='k.csv',
                                    x_file='x.csv',
                                    Qout_file='Qout.nc'
                                    )

    generated_input_file = os.path.join(output_data_path, 
                                      "rapid_namelist-GENERATE-NUMBERS")

    rapid_manager.generate_namelist_file(generated_input_file)
                          
    generated_input_file_solution = os.path.join(compare_data_path, 
                                               "rapid_namelist-GENERATE-NUMBERS-SOLUTION")


    ok_(fcmp(generated_input_file, generated_input_file_solution))
    
    try:
        os.remove(generated_input_file)
    except OSError:
        pass

if __name__ == '__main__':
    import nose
    nose.main()

